// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006  Tiago de Paula Peixoto <tiago@forked.de>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#include <algorithm>
#include <tr1/unordered_set>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/random.hpp>

#include "graph.hh"
#include "histogram.hh"
#include "graph_filtering.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

//==============================================================================
// average_nearest_neighbours_degree
// return generalized average nearest neighbours degree
//==============================================================================

template <class DegreeSelectorOrigin, class DegreeSelectorNeighbours>
struct get_average_nearest_neighbours_degree
{
    get_average_nearest_neighbours_degree(DegreeSelectorOrigin& origin_deg, DegreeSelectorNeighbours& neighbours_deg)
	: _origin_degree(origin_deg), _neighbours_degree(neighbours_deg) {}

    template <class Graph, class AvgDeg>
    void operator()(const Graph& g, AvgDeg& avg_deg) const
    {
	tr1::unordered_map<double,size_t> count;
	tr1::unordered_set<typename graph_traits<Graph>::vertex_descriptor> neighbour_set;

	typename graph_traits<Graph>::vertex_iterator v, v_begin, v_end;
	tie(v_begin,v_end) = vertices(g);
	for(v = v_begin; v != v_end; ++v)
	{
	    typename graph_traits<Graph>::out_edge_iterator e, e_begin, e_end;
	    tie(e_begin,e_end) = out_edges(*v,g);
	    for(e = e_begin; e != e_end; ++e)
	    {		
		if (neighbour_set.find(target(*e,g)) != neighbour_set.end())
		    continue;
		else
		    neighbour_set.insert(target(*e,g));
		typename AvgDeg::value_type::second_type::first_type deg = _neighbours_degree(target(*e,g),g);
		typename AvgDeg::key_type orig_deg = _origin_degree(*v,g);
		avg_deg[orig_deg].first += deg;
		avg_deg[orig_deg].second += deg*deg;
		count[orig_deg]++;
	    }
	    neighbour_set.clear();
	}

	for (typeof(avg_deg.begin()) iter = avg_deg.begin(); iter != avg_deg.end(); ++iter)
	{
	    size_t N = count[iter->first];
	    iter->second.first /= N;
	    if (N > 1)
		iter->second.second = sqrt((iter->second.second - N*iter->second.first*iter->second.first)/(N*(N-1)));
	    else
		iter->second.second = 0.0;
	}
    }
    DegreeSelectorOrigin& _origin_degree;
    DegreeSelectorNeighbours& _neighbours_degree;
};

template <class DegreeSelectors>
struct choose_average_nearest_neighbours_degree
{
    choose_average_nearest_neighbours_degree(const GraphInterface &g, GraphInterface::deg_t origin_deg, GraphInterface::deg_t neighbour_deg, GraphInterface::avg_corr_t &avg_deg)
	: _g(g), _avg_deg(avg_deg) 
    {
	tie(_origin_deg, _origin_deg_name) = get_degree_type(origin_deg);
	tie(_neighbour_deg, _neighbour_deg_name) = get_degree_type(neighbour_deg);
    }

    template <class OriginDegreeSelector>
    struct choose_neighbour_degree
    {
	choose_neighbour_degree(choose_average_nearest_neighbours_degree<DegreeSelectors>& parent):_parent(parent) {}
	template <class DegreeSelector>
	void operator()(DegreeSelector)
	{
	    if ( mpl::at<degree_selector_index, DegreeSelector>::type::value == _parent._neighbour_deg)
	    {
		OriginDegreeSelector origin_deg(_parent._origin_deg_name, _parent._g);
		DegreeSelector deg(_parent._neighbour_deg_name, _parent._g);
		check_filter(_parent._g, bind<void>(get_average_nearest_neighbours_degree<OriginDegreeSelector,DegreeSelector>(origin_deg, deg),
						    _1, var(_parent._avg_deg)),
			     reverse_check(),directed_check()); 
	    }
	}
	choose_average_nearest_neighbours_degree<DegreeSelectors> &_parent;
    };

    template <class DegreeSelector>
    void operator()(DegreeSelector)
    {
	if (mpl::at<degree_selector_index, DegreeSelector>::type::value == _origin_deg)
	    mpl::for_each<DegreeSelectors>(choose_neighbour_degree<DegreeSelector>(*this));
    }

    const GraphInterface &_g;
    GraphInterface::avg_corr_t &_avg_deg;
    GraphInterface::degree_t _origin_deg;
    string _origin_deg_name;
    GraphInterface::degree_t _neighbour_deg;
    string _neighbour_deg_name;
};

//==============================================================================
// GetAverageNearestNeighboursDegree(neigh, orign_deg, neighbours_deg)
//==============================================================================
GraphInterface::avg_corr_t
GraphInterface::GetAverageNearestNeighboursDegree(deg_t origin_deg, deg_t neighbours_deg ) const
{
    GraphInterface::avg_corr_t avg_corr;

    try 
    {
	typedef mpl::vector<in_degreeS,out_degreeS,total_degreeS,scalarS> degrees;
	mpl::for_each<degrees>(choose_average_nearest_neighbours_degree<degrees>(*this, origin_deg, neighbours_deg, avg_corr));
    }
    catch (dynamic_get_failure &e)
    {
	throw GraphException("error getting scalar property: " + string(e.what()));
    }

    return avg_corr;
}
