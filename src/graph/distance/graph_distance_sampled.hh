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
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#ifndef GRAPH_DISTANCE_SAMPLED_HH
#define GRAPH_DISTANCE_SAMPLED_HH

#include <tr1/unordered_set>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/random.hpp>

namespace graph_tool
{
using namespace std;
using namespace boost;

typedef boost::mt19937 rng_t;

// retrieves the histogram of sampled vertex-vertex distances

struct no_weightS {};

struct get_sampled_distances
{

    template <class Graph, class IndexMap, class WeightMap, class Hist>
    void operator()(const Graph* gp, IndexMap index_map, WeightMap weights,
                    Hist& hist, size_t samples, size_t seed) const
    {
        const Graph& g = *gp;
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

        // select get_sum_vertex_dists based on the existence of weights
        typedef typename mpl::if_<is_same<WeightMap, no_weightS>,
                                       get_dists_bfs,
                                  get_dists_djk>::type get_vertex_dists_t;
        get_vertex_dists_t get_vertex_dists;

        tr1::unordered_map<size_t,vertex_t> descriptors;

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
            vertex_t s,t;

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
            if (dist_map[t] != numeric_limits<double>::max() &&
                dist_map[t] != 0.0)
            {
                #pragma omp atomic
                hist[dist_map[t]]++;
            }
        }

    }

    // weighted version. Use dijkstra_shortest_paths()
    struct get_dists_djk
    {
        template <class Graph, class Vertex, class IndexMap, class DistanceMap,
                  class WeightMap>
        void operator()(const Graph& g, Vertex s, IndexMap index_map,
                        DistanceMap dist_map, WeightMap weights) const
        {
            dijkstra_shortest_paths(g, s, vertex_index_map(index_map).
                                    weight_map(weights).distance_map(dist_map));
        }
    };

    // unweighted version. Use BFS.
    struct get_dists_bfs
    {
        template <class Graph, class Vertex, class IndexMap, class DistanceMap>
        void operator()(const Graph& g, Vertex s, IndexMap index_map,
                        DistanceMap dist_map, no_weightS) const
        {
            breadth_first_search(g, s,
                                 visitor(make_bfs_visitor
                                         (record_distances(dist_map,
                                                           on_tree_edge()))));
        }
    };

};

} // graph_tool namespace

#endif // GRAPH_DISTANCE_SAMPLED_HH
