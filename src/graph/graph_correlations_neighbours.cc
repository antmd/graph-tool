// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007  Tiago de Paula Peixoto <tiago@forked.de>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
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
#include "shared_map.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

void operator+=(pair<double, double>&a, const pair<double, double>&b)
{
    a.first += b.first;
    a.second += b.second;
}

//==============================================================================
// average_nearest_neighbours_correlation
// return generalized average nearest neighbours correlation
//==============================================================================

template <class DegreeSelectorOrigin, class DegreeSelectorNeighbours>
struct get_average_nearest_neighbours_correlation
{
    get_average_nearest_neighbours_correlation(DegreeSelectorOrigin& origin_deg, DegreeSelectorNeighbours& neighbours_deg)
        : _origin_degree(origin_deg), _neighbours_degree(neighbours_deg) {}

    template <class Graph, class WeightMap, class AvgDeg>
    void operator()(const Graph& g, WeightMap weight, AvgDeg& avg_deg) const
    {
        tr1::unordered_map<double,double> count;
        SharedMap<tr1::unordered_map<double,double> > s_count(count);
        SharedMap<AvgDeg> s_avg_deg(avg_deg);

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) firstprivate(s_count, s_avg_deg) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            typename AvgDeg::key_type orig_deg = _origin_degree(v,g);

            typename graph_traits<Graph>::out_edge_iterator e, e_begin, e_end;
            tie(e_begin,e_end) = out_edges(v,g);
            for(e = e_begin; e != e_end; ++e)
            {
                typename AvgDeg::value_type::second_type::first_type deg = _neighbours_degree(target(*e,g),g);
                s_avg_deg[orig_deg].first += deg*get(weight, *e);
                s_avg_deg[orig_deg].second += deg*deg;
                s_count[orig_deg] += get(weight,*e);
            }
        }

        s_count.Gather();
        s_avg_deg.Gather();

        for (typeof(avg_deg.begin()) iter = avg_deg.begin(); iter != avg_deg.end(); ++iter)
        {
            double N = count[iter->first];
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

template <class WeightMap, class DegreeSelectors>
struct choose_average_nearest_neighbours_correlation
{
    choose_average_nearest_neighbours_correlation(const GraphInterface &g, WeightMap weight, GraphInterface::deg_t origin_deg, GraphInterface::deg_t neighbour_deg, GraphInterface::avg_corr_t &avg_deg)
        : _g(g), _weight(weight), _avg_deg(avg_deg) 
    {
        tie(_origin_deg, _origin_deg_name) = get_degree_type(origin_deg);
        tie(_neighbour_deg, _neighbour_deg_name) = get_degree_type(neighbour_deg);
    }

    template <class OriginDegreeSelector>
    struct choose_neighbour_degree
    {
        choose_neighbour_degree(choose_average_nearest_neighbours_correlation<WeightMap,DegreeSelectors>& parent):_parent(parent) {}
        template <class DegreeSelector>
        void operator()(DegreeSelector)
        {
            if ( mpl::at<degree_selector_index, DegreeSelector>::type::value == _parent._neighbour_deg)
            {
                OriginDegreeSelector origin_deg(_parent._origin_deg_name, _parent._g);
                DegreeSelector deg(_parent._neighbour_deg_name, _parent._g);
                check_filter(_parent._g, bind<void>(get_average_nearest_neighbours_correlation<OriginDegreeSelector,DegreeSelector>(origin_deg, deg),
                                                    _1, var(_parent._weight), var(_parent._avg_deg)),
                             reverse_check(),directed_check()); 
            }
        }
        choose_average_nearest_neighbours_correlation<WeightMap,DegreeSelectors> &_parent;
    };

    template <class DegreeSelector>
    void operator()(DegreeSelector)
    {
        if (mpl::at<degree_selector_index, DegreeSelector>::type::value == _origin_deg)
            mpl::for_each<DegreeSelectors>(choose_neighbour_degree<DegreeSelector>(*this));
    }

    const GraphInterface &_g;
    WeightMap _weight;
    GraphInterface::avg_corr_t &_avg_deg;
    GraphInterface::degree_t _origin_deg;
    string _origin_deg_name;
    GraphInterface::degree_t _neighbour_deg;
    string _neighbour_deg_name;
};

//==============================================================================
// GetAverageNearestNeighboursCorrelation(neigh, orign_deg, neighbours_deg)
//==============================================================================
GraphInterface::avg_corr_t
GraphInterface::GetAverageNearestNeighboursCorrelation(deg_t origin_deg, deg_t neighbours_deg, std::string weight) const
{
    GraphInterface::avg_corr_t avg_corr;

    try 
    {
        typedef mpl::vector<in_degreeS,out_degreeS,total_degreeS,scalarS> degrees;
        if(weight != "")
        {
            try 
            {
                // FIXME: it would be good also to have a version for a static map (vector_property_map), 
                //        but adding this makes GCC use more than 1 GB of RAM in my system.

                dynamic_property_map& weight_prop = find_property_map(_properties, weight, typeid(graph_traits<multigraph_t>::edge_descriptor));
                typedef DynamicPropertyMapWrap<double,graph_traits<multigraph_t>::edge_descriptor> weight_map_t;
                weight_map_t weight_map(weight_prop);
                mpl::for_each<degrees>(choose_average_nearest_neighbours_correlation<weight_map_t,degrees>(*this, weight_map, origin_deg, neighbours_deg, avg_corr));
            }
            catch (property_not_found& e)
            {
                throw GraphException("error getting scalar property: " + string(e.what()));
            }
        }
        else
        {
            typedef ConstantPropertyMap<double,graph_traits<multigraph_t>::edge_descriptor>  weight_map_t;
            weight_map_t weight_map(1.0); 
            mpl::for_each<degrees>(choose_average_nearest_neighbours_correlation<weight_map_t,degrees>(*this, weight_map, origin_deg, neighbours_deg, avg_corr));
        }
    }
    catch (dynamic_get_failure &e)
    {
        throw GraphException("error getting scalar property: " + string(e.what()));
    }

    return avg_corr;
}
