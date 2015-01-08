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

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_properties.hh"
#include "graph_selectors.hh"

#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <boost/python.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

template <class DistMap>
class bfs_diam_visitor:
    public boost::bfs_visitor<null_visitor>
{
public:
    bfs_diam_visitor(DistMap dist_map, size_t& v)
        : _dist_map(dist_map), _v(v), _dist(0), _max_dist(0),
          _min_k(numeric_limits<size_t>::max()) {}

    template <class Graph>
    void tree_edge(typename graph_traits<Graph>::edge_descriptor e,
                   Graph& g)
    {
        typename graph_traits<Graph>::vertex_descriptor v = target(e,g);
        size_t dist = _dist_map[source(e,g)] + 1;
        if ((dist > _max_dist) ||
            (dist == _max_dist && total_degreeS()(v, g) <= _min_k))
        {
            _max_dist = dist;
            _min_k = total_degreeS()(v, g);
            _v = v;
        }
        _dist_map[v] = dist;
    }

private:
    DistMap _dist_map;
    size_t& _v;
    size_t _dist;
    size_t _max_dist;
    size_t _min_k;
};

template <class DistMap>
class djk_diam_visitor:
    public boost::dijkstra_visitor<null_visitor>
{
public:
    djk_diam_visitor(DistMap dist_map, size_t& v)
        : _dist_map(dist_map), _v(v), _max_dist(0),
          _min_k(numeric_limits<size_t>::max()) {}

    template <class Graph>
    void examine_vertex(typename graph_traits<Graph>::vertex_descriptor v,
                        Graph& g)
    {
        if ((_dist_map[v] > _max_dist) ||
            (_dist_map[v] == _max_dist && total_degreeS()(v, g) <= _min_k))
        {
            _max_dist = _dist_map[v];
            _min_k = total_degreeS()(v, g);
            _v = v;
        }
    }

private:
    DistMap _dist_map;
    size_t& _v;
    typename property_traits<DistMap>::value_type _max_dist;
    size_t _min_k;


};


struct do_bfs_search
{
    template <class Graph, class VertexIndexMap>
    void operator()(const Graph& g, size_t source, VertexIndexMap vertex_index,
                    size_t& target, long double& max_dist) const
    {
        typedef unchecked_vector_property_map<size_t, VertexIndexMap> dist_map_t;
        dist_map_t dist_map(vertex_index, num_vertices(g));
        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) schedule(runtime) if (N > 100)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            dist_map[v] = numeric_limits<size_t>::max();
        }
        dist_map[vertex(source,g)] = 0;

        unchecked_vector_property_map<boost::default_color_type, VertexIndexMap>
            color_map(vertex_index, num_vertices(g));
        target = source;
        breadth_first_search(g, vertex(source, g),
                             visitor(bfs_diam_visitor<dist_map_t>
                                      (dist_map, std::ref(target))).
                             vertex_index_map(vertex_index).
                             color_map(color_map));
        max_dist = dist_map[vertex(target, g)];
    }
};

struct do_djk_search
{
    template <class Graph, class VertexIndexMap, class WeightMap>
    void operator()(const Graph& g, size_t source, VertexIndexMap vertex_index,
                    WeightMap weight, size_t& target, long double& max_dist) const
    {
        typedef unchecked_vector_property_map<typename property_traits<WeightMap>::value_type,
                                              VertexIndexMap> dist_map_t;
        dist_map_t dist_map(vertex_index, num_vertices(g));
        target = source;
        dijkstra_shortest_paths(g, vertex(source, g),
                                weight_map(weight).
                                distance_map(dist_map).
                                vertex_index_map(vertex_index).
                                visitor(djk_diam_visitor<dist_map_t>
                                        (dist_map, std::ref(target))));
        max_dist = dist_map[vertex(target, g)];
    }
};

python::object get_diam(GraphInterface& gi, size_t source, boost::any weight)
{
    size_t target;
    long double max_dist;
    if (weight.empty())
    {
        run_action<>()
            (gi, std::bind(do_bfs_search(), placeholders::_1, source, gi.GetVertexIndex(),
                           std::ref(target), std::ref(max_dist)))();
    }
    else
    {
        run_action<>()
            (gi, std::bind(do_djk_search(), placeholders::_1, source, gi.GetVertexIndex(),
                           placeholders::_2, std::ref(target), std::ref(max_dist)),
             edge_scalar_properties())(weight);

    }
    return python::make_tuple(target, max_dist);
}

void export_diam()
{
    python::def("get_diam", &get_diam);
};
