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

#ifndef GRAPH_DISTANCE_HH
#define GRAPH_DISTANCE_HH

#include <tr1/unordered_set>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

namespace graph_tool
{
using namespace std;
using namespace boost;

// retrieves the vertex-vertex distance histogram

struct no_weightS {};

struct get_distance_histogram
{

    template <class Graph, class IndexMap, class WeightMap, class Hist>
    void operator()(const Graph *gp, IndexMap index_map, WeightMap weights,
                    Hist& hist) const
    {
        const Graph& g = *gp;
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

        // select get_vertex_dists based on the existence of weights
        typedef typename mpl::if_<is_same<WeightMap, no_weightS>,
                                       get_dists_bfs,
                                  get_dists_djk>::type get_vertex_dists_t;

        get_vertex_dists_t get_vertex_dists;
        int i, N = num_vertices(g);

        #pragma omp parallel for default(shared) private(i)
        for (i = 0; i < N; ++i)
        {
            vertex_t v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            typedef tr1::unordered_map<vertex_t,double,
                                       DescriptorHash<IndexMap> > dmap_t;
            dmap_t dmap(0, DescriptorHash<IndexMap>(index_map));
            InitializedPropertyMap<dmap_t>
                dist_map(dmap, numeric_limits<double>::max());

            dist_map[v] = 0.0;
            get_vertex_dists(g, v, index_map, dist_map, weights);

            typename graph_traits<Graph>::vertex_iterator v2, v_end;
            for (tie(v2, v_end) = vertices(g); v2 != v_end; ++v2)
                if (*v2 != v && dist_map[*v2] != numeric_limits<double>::max())
                {
                    double dist = dist_map[*v2];
                    #pragma omp atomic
                    hist[dist]++;
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
            typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
            typedef tr1::unordered_map<vertex_t,default_color_type,
                                       DescriptorHash<IndexMap> > cmap_t;
            cmap_t cmap(0, DescriptorHash<IndexMap>(index_map));
            InitializedPropertyMap<cmap_t>
                color_map(cmap, color_traits<default_color_type>::white());

            breadth_first_visit(g, s,
                                visitor(make_bfs_visitor
                                        (record_distances(dist_map,
                                                          on_tree_edge()))).
                                color_map(color_map));
        }
    };
};

} // boost namespace

#endif // GRAPH_DISTANCE_HH
