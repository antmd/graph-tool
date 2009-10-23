// Copyright (C) 2008  Tiago de Paula Peixoto <tiago@forked.de>
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

struct stop_search {};

template <class DistMap, class PredMap>
class bfs_max_visitor:
    public boost::bfs_visitor<null_visitor>
{
public:
    bfs_max_visitor(DistMap dist_map, PredMap pred, size_t max_dist)
        : _dist_map(dist_map), _pred(pred), _max_dist(max_dist), _dist(0) {}


    template <class Graph>
    void examine_edge(typename graph_traits<Graph>::edge_descriptor e,
                      Graph& g)
    {
        _pred[target(e,g)] = source(e,g);
    }

    template <class Graph>
    void discover_vertex(typename graph_traits<Graph>::vertex_descriptor v,
                         Graph& g)
    {
        if (_pred[v] == v)
            return;
        size_t dist = _dist_map[_pred[v]] + 1;
        if (dist > _max_dist)
            throw stop_search();
        _dist_map[v] = dist;
    }

private:
    DistMap _dist_map;
    PredMap _pred;
    size_t _max_dist;
    size_t _dist;
};

template <class DistMap, class PredMap>
class djk_max_visitor:
    public boost::dijkstra_visitor<null_visitor>
{
public:
    djk_max_visitor(DistMap dist_map, PredMap pred_map,
                    typename property_traits<DistMap>::value_type max_dist)
        : _dist_map(dist_map), _max_dist(max_dist) {}

    template <class Graph>
    void examine_vertex(typename graph_traits<Graph>::vertex_descriptor u,
                        Graph& g)
    {
        if (_dist_map[_pred_map[u]] >= _max_dist)
            throw stop_search();
    }

private:
    DistMap _dist_map;
    PredMap _pred_map;
    typename property_traits<DistMap>::value_type _max_dist;
};


struct do_bfs_search
{
    template <class Graph, class VertexIndexMap, class DistMap>
    void operator()(const Graph& g, size_t source, VertexIndexMap vertex_index,
                    DistMap dist_map, long double max_dist) const
    {
        typedef typename property_traits<DistMap>::value_type dist_t;
        dist_t max_d = (max_dist > 0) ?
            max_dist : numeric_limits<dist_t>::max();

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            dist_map[v] = numeric_limits<dist_t>::max();
        }
        dist_map[vertex(source,g)] = 0;

        typedef unchecked_vector_property_map
            <typename graph_traits<Graph>::vertex_descriptor, VertexIndexMap>
            pred_map_t;
        pred_map_t pred_map(vertex_index, num_vertices(g));
        pred_map[vertex(source, g)] = vertex(source, g);
        unchecked_vector_property_map<boost::default_color_type, VertexIndexMap>
            color_map(vertex_index, num_vertices(g));
        try
        {
            breadth_first_search(g, vertex(source, g),
                                 visitor(bfs_max_visitor<DistMap, pred_map_t>
                                         (dist_map, pred_map, max_d)).
                                 vertex_index_map(vertex_index).
                                 color_map(color_map));
        }
        catch (stop_search&) {}
    }
};

struct do_djk_search
{
    template <class Graph, class VertexIndexMap, class DistMap, class WeightMap>
    void operator()(const Graph& g, size_t source, VertexIndexMap vertex_index,
                    DistMap dist_map, WeightMap weight, long double max_dist)
        const
    {
        typedef typename property_traits<DistMap>::value_type dist_t;
        dist_t max_d = (max_dist > 0) ?
            max_dist : numeric_limits<dist_t>::max();

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            dist_map[v] = numeric_limits<dist_t>::max();
        }
        dist_map[vertex(source,g)] = 0;
        typedef unchecked_vector_property_map
            <typename graph_traits<Graph>::vertex_descriptor, VertexIndexMap>
            pred_map_t;
        pred_map_t pred_map(vertex_index, num_vertices(g));
        unchecked_vector_property_map<boost::default_color_type, VertexIndexMap>
            color_map(vertex_index, num_vertices(g));
        try
        {
            dijkstra_shortest_paths(g, vertex(source, g),
                                    visitor(djk_max_visitor<DistMap, pred_map_t>
                                            (dist_map, pred_map, max_d)).
                                    weight_map(weight).
                                    distance_map(dist_map).
                                    predecessor_map(pred_map).
                                    vertex_index_map(vertex_index).
                                    color_map(color_map));
        }
        catch (stop_search&) {}
    }
};

void get_dists(GraphInterface& gi, size_t source, boost::any dist_map,
               boost::any weight, long double max_dist)
{
    if (weight.empty())
    {
        run_action<>()
            (gi, bind<void>(do_bfs_search(), _1, source, gi.GetVertexIndex(),
                            _2, max_dist),
             writable_vertex_scalar_properties())
            (dist_map);
    }
    else
    {
        run_action<>()
            (gi, bind<void>(do_djk_search(), _1, source, gi.GetVertexIndex(),
                            _2, _3, max_dist),
             writable_vertex_scalar_properties(),
             edge_scalar_properties())
            (dist_map, weight);

    }
}

void export_dists()
{
    python::def("get_dists", &get_dists);
};

