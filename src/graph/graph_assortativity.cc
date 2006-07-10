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
// GetAssortativityCoefficient(neigh, deg1, deg2)
//==============================================================================

template <class EdgeSelector, class DegreeSelector1, class DegreeSelector2>
struct get_assortativity_coefficient
{
    get_assortativity_coefficient(DegreeSelector1& deg1, DegreeSelector2& deg2): _deg1(deg1), _deg2(deg2) {}

    template <class Graph>
    void operator()(const Graph &g, double &a) const
    {
	tr1::unordered_map<double,int> e_kk;
	size_t n_edges = 0, n_vertices = 0;
	tr1::unordered_map<double,int> n_k;
	
	typename graph_traits<Graph>::vertex_iterator v, v_begin, v_end;
	tie(v_begin,v_end) = vertices(g);
	for(v = v_begin; v != v_end; ++v)
	{
	    double k1, k2;
	    k1 = _deg1(*v,g);
	    n_k[k1]++;
	    n_vertices++;
	    
	    typename EdgeSelector::template iterator<Graph>::type e, e_begin, e_end;
	    tie(e_begin, e_end) = _edges(*v,g);
	    for (e = e_begin; e != e_end; ++e)
	    {
		k2 = _deg2(_edges.target(*e,g),g);
		if (k1 == k2)
		    e_kk[k1]++;
		n_edges++;
	    }
	}
	
	double t1=0.0, t2=0.0;
	
	for (typeof(e_kk.begin()) iter = e_kk.begin(); iter != e_kk.end(); ++iter)
	    t1 += double(iter->second)/n_edges;
	
	double avg_k = GetHistogramMean(n_k);
	
	for (typeof(n_k.begin()) iter = n_k.begin(); iter != n_k.end(); ++iter)
	    t2 += double(iter->second * iter->second * iter->first * iter->first)/(avg_k*avg_k*n_vertices*n_vertices);
	
	a = (t1 == 1.0)?1.0:(t1 - t2)/(1.0 - t2);
    }
    
    EdgeSelector _edges;
    DegreeSelector1& _deg1;
    DegreeSelector2& _deg2;
};

template <class DegreeSelectors>
struct choose_assortativity_coefficient
{
    choose_assortativity_coefficient(const GraphInterface &g, GraphInterface::neighbours_t neigh,
				     GraphInterface::deg_t deg1, GraphInterface::deg_t deg2, double &a)
	: _g(g), _neigh(neigh), _a(a) 
    {
	tie(_deg1, _deg_name1) = get_degree_type(deg1);
	tie(_deg2, _deg_name2) = get_degree_type(deg2);
    }

    template <class NeighbourSelector, class DegreeSelector1>
    struct choose_degree_two
    {
	choose_degree_two(choose_assortativity_coefficient<DegreeSelectors> &parent):_parent(parent) {}
	template <class DegreeSelector2>
	void operator()(DegreeSelector2)
	{
	    if (mpl::at<degree_selector_index, DegreeSelector2>::type::value == _parent._deg2)
	    {
		DegreeSelector1 deg1(_parent._deg_name1, _parent._g);
		DegreeSelector2 deg2(_parent._deg_name2, _parent._g);
		check_filter(_parent._g, bind<void>(get_assortativity_coefficient<NeighbourSelector,DegreeSelector1,DegreeSelector2>(deg1,deg2), 
						    _1, var(_parent._a)), 
			     reverse_check(),directed_check());
	    }
	}
	choose_assortativity_coefficient<DegreeSelectors> &_parent;
    };

    template <class NeighbourSelector>
    struct choose_degree_one
    {
	choose_degree_one(choose_assortativity_coefficient<DegreeSelectors> &parent):_parent(parent) {}
	template <class DegreeSelector1>
	void operator()(DegreeSelector1)
	{
	    if (mpl::at<degree_selector_index, DegreeSelector1>::type::value == _parent._deg1)
		mpl::for_each<DegreeSelectors>(choose_degree_two<NeighbourSelector, DegreeSelector1>(_parent));
	}
	choose_assortativity_coefficient<DegreeSelectors> &_parent;
    };

    template <class NeighbourSelector>
    void operator()(NeighbourSelector)
    {
	if (mpl::at<edges_selector_index, NeighbourSelector>::type::value == _neigh)
	    mpl::for_each<DegreeSelectors>(choose_degree_one<NeighbourSelector>(*this));
    }

    const GraphInterface &_g;
    GraphInterface::neighbours_t _neigh;
    double &_a;
    GraphInterface::degree_t _deg1;
    string _deg_name1;
    GraphInterface::degree_t _deg2;
    string _deg_name2;
};


double 
GraphInterface::GetAssortativityCoefficient(GraphInterface::neighbours_t neigh, GraphInterface::deg_t deg1, GraphInterface::deg_t deg2 ) const
{
    double a;
    try
    {
	typedef mpl::vector<out_edgeS, in_edgeS, any_edgeS> neighbours;
	typedef mpl::vector<in_degreeS, out_degreeS, total_degreeS, scalarS> degrees;
	mpl::for_each<neighbours>(choose_assortativity_coefficient<degrees>(*this, neigh, deg1, deg2, a));
    }
    catch (dynamic_get_failure &e)
    {
	throw GraphException("error getting scalar property: " + string(e.what()));
    }

    return a;
}
