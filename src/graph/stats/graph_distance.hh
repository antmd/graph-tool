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

#ifndef GRAPH_DISTANCE_HH
#define GRAPH_DISTANCE_HH

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
using namespace boost;

// retrieves the vertex-vertex distance histogram

struct no_weightS {};

template <class Map>
struct get_val_type
{
    typedef typename property_traits<Map>::value_type type;
};

template <>
struct get_val_type<no_weightS>
{
    typedef size_t type;
};

struct get_distance_histogram
{

    template <class Graph, class VertexIndex, class WeightMap>
    void operator()(const Graph& g, VertexIndex vertex_index, WeightMap weights,
                    const vector<long double>& obins, python::object& phist)
        const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

        // select get_vertex_dists based on the existence of weights
        typedef typename mpl::if_<std::is_same<WeightMap, no_weightS>,
                                  get_dists_bfs,
                                  get_dists_djk>::type get_vertex_dists_t;

        // distance type
        typedef typename get_val_type<WeightMap>::type val_type;
        typedef Histogram<val_type, size_t, 1> hist_t;

        std::array<vector<val_type>,1> bins;
        bins[0].resize(obins.size());
        for (size_t i = 0; i < obins.size(); ++i)
            bins[0][i] = obins[i];

        hist_t hist(bins);
        SharedHistogram<hist_t> s_hist(hist);

        typename hist_t::point_t point;
        get_vertex_dists_t get_vertex_dists;
        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i,point) \
            firstprivate(s_hist) schedule(runtime) if (N > 100)
        for (i = 0; i < N; ++i)
        {
            vertex_t v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            unchecked_vector_property_map<val_type,VertexIndex>
                dist_map(vertex_index, num_vertices(g));

            for (size_t j = 0; j < size_t(N); ++j)
            {
                if (vertex(j,g) != graph_traits<Graph>::null_vertex())
                    dist_map[vertex(j,g)] =  numeric_limits<val_type>::max();
            }

            dist_map[v] = 0;
            get_vertex_dists(g, v, vertex_index, dist_map, weights);

            typename graph_traits<Graph>::vertex_iterator v2, v_end;
            for (tie(v2, v_end) = vertices(g); v2 != v_end; ++v2)
                if (*v2 != v &&
                    dist_map[*v2] != numeric_limits<val_type>::max())
                {
                    point[0] = dist_map[*v2];
                    s_hist.PutValue(point);
                }
        }
        s_hist.Gather();

        python::list ret;
        ret.append(wrap_multi_array_owned<size_t,1>(hist.GetArray()));
        ret.append(wrap_vector_owned<val_type>(hist.GetBins()[0]));
        phist = ret;
    }

    // weighted version. Use dijkstra_shortest_paths()
    struct get_dists_djk
    {
        template <class Graph, class Vertex, class VertexIndex,
                  class DistanceMap, class WeightMap>
        void operator()(const Graph& g, Vertex s, VertexIndex vertex_index,
                        DistanceMap dist_map, WeightMap weights) const
        {
            dijkstra_shortest_paths(g, s, vertex_index_map(vertex_index).
                                    weight_map(weights).distance_map(dist_map));
        }
    };

    // unweighted version. Use BFS.
    struct get_dists_bfs
    {
        template <class Graph, class Vertex, class VertexIndex,
                  class DistanceMap>
        void operator()(const Graph& g, Vertex s, VertexIndex vertex_index,
                        DistanceMap dist_map, no_weightS) const
        {
            typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
            typedef unordered_map<vertex_t,default_color_type,
                                  DescriptorHash<VertexIndex> > cmap_t;
            cmap_t cmap(0, DescriptorHash<VertexIndex>(vertex_index));
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
