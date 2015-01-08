// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2015 Tiago de Paula Peixoto <tiago@skewed.de>
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

#ifndef GRAPH_CLOSENESS_HH
#define GRAPH_CLOSENESS_HH

#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <boost/python/object.hpp>
#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>

#include "histogram.hh"
#include "numpy_bind.hh"

namespace graph_tool
{
using namespace std;

struct no_weightS {};

template <class Map>
struct get_val_type
{
    typedef typename boost::property_traits<Map>::value_type type;
};

template <>
struct get_val_type<no_weightS>
{
    typedef size_t type;
};

struct get_closeness
{
    typedef void result_type;
    template <class Graph, class VertexIndex, class WeightMap, class Closeness>
    void operator()(const Graph& g, VertexIndex vertex_index, WeightMap weights,
                    Closeness closeness, bool harmonic, bool norm)
        const
    {
        using namespace boost;
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

        // select get_vertex_dists based on the existence of weights
        typedef typename mpl::if_<std::is_same<WeightMap, no_weightS>,
                                  get_dists_bfs,
                                  get_dists_djk>::type get_vertex_dists_t;

        // distance type
        typedef typename get_val_type<WeightMap>::type val_type;

        get_vertex_dists_t get_vertex_dists;
        size_t HN = HardNumVertices()(g);
        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) schedule(runtime) if (N > 100)
        for (i = 0; i < N; ++i)
        {
            vertex_t v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            unchecked_vector_property_map<val_type,VertexIndex>
                dist_map(vertex_index, num_vertices(g));

            for (int j = 0; j < N; ++j)
            {
                if (vertex(j, g) != graph_traits<Graph>::null_vertex())
                    dist_map[vertex(j, g)] = numeric_limits<val_type>::max();
            }

            dist_map[v] = 0;

            size_t comp_size = 0;
            get_vertex_dists(g, v, vertex_index, dist_map, weights, comp_size);

            closeness[v] = 0;
            typename graph_traits<Graph>::vertex_iterator v2, v_end;
            for (tie(v2, v_end) = vertices(g); v2 != v_end; ++v2)
            {
                if (*v2 != v && dist_map[*v2] != numeric_limits<val_type>::max())
                {
                    if (!harmonic)
                        closeness[v] += dist_map[*v2];
                    else
                        closeness[v] += 1. / dist_map[*v2];
                }
            }
            if (!harmonic)
                closeness[v] = 1 / closeness[v];
            if (norm)
            {
                if (harmonic)
                    closeness[v] /= HN - 1;
                else
                    closeness[v] *= comp_size - 1;
            }
         }
    }

    class component_djk_visitor: public boost::dijkstra_visitor<>
    {
    public:
        //component_visitor() { }
        component_djk_visitor(size_t& comp_size)
            : _comp_size(comp_size) { }

        template <class Vertex, class Graph>
        void discover_vertex(Vertex u, const Graph&)
        {
            ++_comp_size;
        }

    private:
        size_t& _comp_size;
    };

    // weighted version. Use dijkstra_shortest_paths()
    struct get_dists_djk
    {
        template <class Graph, class Vertex, class VertexIndex,
                  class DistanceMap, class WeightMap>
        void operator()(const Graph& g, Vertex s, VertexIndex vertex_index,
                        DistanceMap dist_map, WeightMap weights,
                        size_t& comp_size) const
        {
            using namespace boost;
            component_djk_visitor vis(comp_size);
            dijkstra_shortest_paths(g, s, vertex_index_map(vertex_index).
                                    weight_map(weights).distance_map(dist_map).visitor(vis));
        }
    };

    template <class DistMap>
    class component_bfs_visitor: public boost::bfs_visitor<>
    {
    public:
        //component_visitor() { }
        component_bfs_visitor(DistMap dist_map, size_t& comp_size)
            : _dist_map(dist_map), _comp_size(comp_size) { }

        template <class Vertex, class Graph>
        void discover_vertex(Vertex, const Graph&)
        {
            ++_comp_size;
        }

        template <class Edge, class Graph>
        void tree_edge(Edge e, const Graph& g)
        {
            _dist_map[target(e, g)] = _dist_map[source(e, g)] + 1;
        }

    private:
        DistMap _dist_map;
        size_t& _comp_size;
    };


    // unweighted version. Use BFS.
    struct get_dists_bfs
    {
        template <class Graph, class Vertex, class VertexIndex,
                  class DistanceMap>
        void operator()(const Graph& g, Vertex s, VertexIndex vertex_index,
                        DistanceMap dist_map, no_weightS, size_t& comp_size) const
        {
            using namespace boost;
            typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
            typedef unordered_map<vertex_t,default_color_type,
                                  DescriptorHash<VertexIndex> > cmap_t;
            cmap_t cmap(0, DescriptorHash<VertexIndex>(vertex_index));
            InitializedPropertyMap<cmap_t>
                color_map(cmap, color_traits<default_color_type>::white());
            component_bfs_visitor<DistanceMap> vis(dist_map, comp_size);
            breadth_first_visit(g, s, visitor(vis).
                                color_map(color_map));
        }
    };
};

} // boost namespace

#endif // GRAPH_CLOSENESS_HH
