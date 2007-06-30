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
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
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

typedef boost::mt19937 rng_t;

//==============================================================================
// GetSampledDistanceHistogram()
// retrieves the histogram of sampled vertex-vertex distances
//==============================================================================
struct no_weightS {};

struct get_sampled_distances
{

    template <class Graph, class IndexMap, class WeightMap, class Hist>
    void operator()(const Graph &g, IndexMap index_map, WeightMap weights, Hist& hist, size_t samples, size_t seed) const
    {        
        // select get_sum_vertex_dists based on the existence of weights
        typedef typename mpl::if_<is_same<WeightMap, no_weightS>,
                                       get_dists_bfs,
                                  get_dists_djk>::type get_vertex_dists_t;
        get_vertex_dists_t get_vertex_dists;

        tr1::unordered_map<size_t, typename graph_traits<Graph>::vertex_descriptor> descriptors;

        typename graph_traits<Graph>::vertex_iterator v, v_end;
        int i = 0, N = 0;
        for(tie(v, v_end) = vertices(g); v != v_end; ++v,++i)
        {
            descriptors[i] = *v;
            N++;
        }

        rng_t rng(static_cast<rng_t::result_type>(seed));
        uniform_int<size_t> sampler(0,descriptors.size()-1);

        #pragma omp parallel for default(shared) private(i,v,v_end) 
        for(i=0; i < int(samples); ++i)
        {
            typedef HashedDescriptorMap<IndexMap,double> dist_map_t;
            dist_map_t dist_map(index_map);

            for(tie(v, v_end) = vertices(g); v != v_end; ++v)
                dist_map[*v] = numeric_limits<double>::max();
            typename graph_traits<Graph>::vertex_descriptor s,t;

            #pragma omp critical
            {
                s = descriptors[sampler(rng)];
                do
                {
                    t = descriptors[sampler(rng)];
                }
                while (t == s && N != 1);
            }

            dist_map[s] = 0.0;
            get_vertex_dists(g, s, index_map, dist_map, weights);
            if (dist_map[t] != numeric_limits<double>::max() && dist_map[t] != 0.0)
            {
                #pragma omp atomic
                hist[dist_map[t]]++;
            }
        }
        
    }

    // weighted version. Use dijkstra_shortest_paths()
    struct get_dists_djk
    {
        template <class Graph, class Vertex, class IndexMap, class DistanceMap, class WeightMap>
        void operator()(const Graph& g, Vertex s, IndexMap index_map, DistanceMap dist_map, WeightMap weights) const
        {
            dijkstra_shortest_paths(g, s, vertex_index_map(index_map).weight_map(weights).distance_map(dist_map));            
        }
    };

    // unweighted version. Use BFS.
    struct get_dists_bfs
    {
        template <class Graph, class Vertex, class IndexMap, class DistanceMap>
        void operator()(const Graph& g, Vertex s, IndexMap index_map, DistanceMap dist_map, no_weightS) const
        {            
            breadth_first_search(g, s, visitor(make_bfs_visitor(record_distances(dist_map, on_tree_edge()))));
        }
    };
    
};

GraphInterface::hist_t 
GraphInterface::GetSampledDistanceHistogram(string weight, size_t samples, size_t seed) const
{
    hist_t hist;

    if (weight == "")
    {
        check_filter(*this, bind<void>(get_sampled_distances(), _1, _vertex_index, no_weightS(), var(hist), samples, seed),
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
                check_filter(*this, bind<void>(get_sampled_distances(), _1, _vertex_index, weight_map, var(hist), samples, seed),
                             reverse_check(), directed_check()); 
            }
            catch (bad_cast)
            {
                DynamicPropertyMapWrap<double, graph_traits<multigraph_t>::edge_descriptor> weight_map(weight_prop);
                check_filter(*this, bind<void>(get_sampled_distances(), _1, _vertex_index, weight_map, var(hist), samples, seed),
                             reverse_check(), directed_check()); 
            }
        }
        catch (property_not_found& e)
        {
            throw GraphException("error getting scalar property: " + string(e.what()));
        }
    }
    return hist;
}
