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
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/lambda/bind.hpp>

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
// bfs_distance_sum_visitor
// This event visitor will record and sum all the distances during a BFS.
//==============================================================================

struct normal_distance 
{
    double operator()(double d) const { return d; }
};

struct harmonic_distance 
{
    double operator()(double d) const { return 1.0/d; }
};

template <class DistanceSelector, class DistanceMap, class Graph> 
class bfs_distance_sum_visitor: public default_bfs_visitor
{
public:    
    bfs_distance_sum_visitor(DistanceMap distance_map, double &sum)
	:_distmap(distance_map), _distsum(sum) { }

    typedef typename graph_traits<Graph>::edge_descriptor edge_t;
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

    void tree_edge(edge_t e, const Graph & g) 
    {
        double d = _distmap[source(e,g)] + 1.0;
        _distmap[target(e,g)] = d; // record distance
        _distsum += _distance(d);
    }

    void initialize_vertex(vertex_t u, const Graph &g)
    {
       _distmap[u] = 0.0;
    }
        
private:
    DistanceMap _distmap;
    double &_distsum;
    DistanceSelector _distance;
};

//==============================================================================
// GetAverageDistance()
// retrieves the average vertex-vertex distance
//==============================================================================
struct no_weightS {};

template <class DistanceSelector>
struct get_average_distance
{

    template <class Graph, class IndexMap, class WeightMap>
    void operator()(const Graph &g, IndexMap index_map, WeightMap weights, double &dist) const
    {	
	typedef DescriptorHash<typename graph_traits<Graph>::vertex_descriptor,IndexMap> hashfc_t;
	typedef tr1::unordered_map<typename graph_traits<Graph>::vertex_descriptor, double, hashfc_t>  map_t;
	hashfc_t hasher(index_map);
	map_t vertex_to_dist(0,hasher);
	typedef associative_property_map<map_t> dist_map_t;
	dist_map_t dist_map(vertex_to_dist);

	double distsum = 0;
	size_t n = 0;
	// select get_sum_vertex_dists based on the existence of weights
	typedef typename mpl::if_<is_same<WeightMap, no_weightS>,
                     	          get_sum_dists_bfs,
                                  get_sum_dists_djk>::type get_sum_vertex_dists;
	get_sum_vertex_dists get_sum_dists;
	typename graph_traits<Graph>::vertex_iterator v, v_begin, v_end;
	tie(v_begin, v_end) = vertices(g);
	for(v = v_begin; v != v_end; ++v)
	{
	    distsum += get_sum_dists(g, *v, index_map, dist_map, weights);
	    n++;
	}
	DistanceSelector distance;
	dist = n<2?0.0:distance(distsum/(n*(n-1)));
    }


    // weighted version. Use dijkstra_shortest_paths()
    struct get_sum_dists_djk
    {
	template <class Graph, class Vertex, class IndexMap, class DistanceMap, class WeightMap>
	double operator()(const Graph& g, Vertex s, IndexMap index_map, DistanceMap dist_map, WeightMap weights) const
	{
	    double distsum = 0.0;
	    dijkstra_shortest_paths(g, s, vertex_index_map(index_map).weight_map(weights).distance_map(dist_map));
	    DistanceSelector get_dist;
	    typename graph_traits<Graph>::vertex_iterator v, v_begin, v_end;
	    tie(v_begin, v_end) = vertices(g);
	    for(v = v_begin; v != v_end; ++v)
		if (dist_map[*v] != std::numeric_limits<double>::max() && *v != s)
		    distsum += get_dist(dist_map[*v]);
	    return distsum;
	}
    };

    // unweighted version. Use BFS.
    struct get_sum_dists_bfs
    {
	template <class Graph, class Vertex, class IndexMap, class DistanceMap>
	double operator()(const Graph& g, Vertex s, IndexMap index_map, DistanceMap dist_map, no_weightS) const
	{
	    double distsum = 0.0;
	    bfs_distance_sum_visitor<DistanceSelector,DistanceMap,Graph> bfs_sum_dists(dist_map, distsum); 
	    breadth_first_search(g, s, visitor(bfs_sum_dists));
	    return distsum;
	}
    };

    
};

double GraphInterface::GetAverageDistance(string weight) const
{
    double avg_dist = 0;
    if (weight == "")
    {
	check_filter(*this, bind<void>(get_average_distance<normal_distance>(), _1, _vertex_index, no_weightS(), var(avg_dist)),
		     reverse_check(), directed_check()); 
    }
    else
    {
	try 
	{
	    dynamic_property_map& weight_prop = find_property_map(_properties, weight, typeid(graph_traits<multigraph_t>::edge_descriptor));
	    try 
	    {
		vector_property_map<double, edge_index_map_t> weight_map;
		weight_map = get_static_property_map<vector_property_map<double, edge_index_map_t> >(weight_prop);
		check_filter(*this, bind<void>(get_average_distance<normal_distance>(), _1, _vertex_index, weight_map, var(avg_dist)),
			     reverse_check(), directed_check()); 
	    }
	    catch (bad_cast)
	    {
		DynamicPropertyMapWrap<double, graph_traits<multigraph_t>::edge_descriptor> weight_map(weight_prop);
		check_filter(*this, bind<void>(get_average_distance<normal_distance>(), _1, _vertex_index, weight_map, var(avg_dist)),
			     reverse_check(), directed_check()); 
	    }
	}
	catch (property_not_found& e)
	{
	    throw GraphException("error getting scalar property: " + string(e.what()));
	}
    }
    return avg_dist;
}

//==============================================================================
// GetAverageHarmonicDistance()
// retrieves the average vertex-vertex harmonic distance
//==============================================================================

double GraphInterface::GetAverageHarmonicDistance(string weight) const
{
    double avg_dist = 0;
    if (weight == "")
    {
	check_filter(*this, bind<void>(get_average_distance<harmonic_distance>(), _1, _vertex_index, no_weightS(), var(avg_dist)),
		     reverse_check(), directed_check()); 
    }
    else
    {
	try 
	{
	    dynamic_property_map& weight_prop = find_property_map(_properties, weight, typeid(graph_traits<multigraph_t>::edge_descriptor));
	    try 
	    {
		vector_property_map<double, edge_index_map_t> weight_map;
		weight_map = get_static_property_map<vector_property_map<double, edge_index_map_t> >(weight_prop);
		check_filter(*this, bind<void>(get_average_distance<harmonic_distance>(), _1, _vertex_index, weight_map, var(avg_dist)),
			     reverse_check(), directed_check()); 
	    }
	    catch (bad_cast)
	    {
		DynamicPropertyMapWrap<double, graph_traits<multigraph_t>::edge_descriptor> weight_map(weight_prop);
		check_filter(*this, bind<void>(get_average_distance<harmonic_distance>(), _1, _vertex_index, weight_map, var(avg_dist)),
			     reverse_check(), directed_check()); 
	    }
	}
	catch (property_not_found& e)
	{
	    throw GraphException("error getting scalar property: " + string(e.what()));
	}
    }
    return avg_dist;
}
