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

//  As a special exception, you have permission to link this program
//  with the CGAL library and distribute executables, as long as you
//  follow the requirements of the GNU GPL in regard to all of the
//  software in the executable aside from CGAL.

#ifndef GRAPH_TRIANGULATION_HH
#define GRAPH_TRIANGULATION_HH

#include <unordered_set>
#include <tuple>

#include <boost/functional/hash.hpp>

#include "graph_util.hh"

namespace graph_tool
{
using namespace std;
using namespace boost;


template <class Graph>
bool is_adjacent(typename graph_traits<Graph>::vertex_descriptor v1,
                 typename graph_traits<Graph>::vertex_descriptor v2,
                 Graph& g)
{
    typename graph_traits<Graph>::out_edge_iterator e, e_end;
    for (tie(e, e_end) = out_edges(v1, g); e != e_end; ++e)
        if (target(*e, g) == v2)
            return true;
    return false;
}

struct hash_point
{
    template <class Vertex>
    std::size_t operator()(const Vertex& v) const
    {
        size_t seed = 42;
        hash_combine(seed, v.point().x());
        hash_combine(seed, v.point().y());
        hash_combine(seed, v.point().z());
        return seed;
    }
};

template <class Triang>
struct get_triangulation
{
    // this will insert edges in the graph
    template <class Graph, class VertexMap>
    class edge_inserter
    {
    public:
        typedef output_iterator_tag iterator_category;
        typedef typename graph_traits<Graph>::vertex_descriptor value_type;
        typedef size_t difference_type;
        typedef typename graph_traits<Graph>::vertex_descriptor* pointer;
        typedef typename graph_traits<Graph>::vertex_descriptor& reference;

        edge_inserter(Graph& g, const typename VertexMap::key_type& v,
                      VertexMap& vertex_map): _g(g), _vertex_map(vertex_map),
                                              _source(vertex_map[v]) {}

        edge_inserter& operator++() { return *this; }
        edge_inserter& operator++(int) { return *this; }
        edge_inserter& operator*() { return *this; }

        template <class Vertex>
        edge_inserter& operator=(const Vertex& v)
        {
            typename VertexMap::iterator iter = _vertex_map.find(*v);
            if (iter != _vertex_map.end())
            {
                typename graph_traits<Graph>::vertex_descriptor target
                    = iter->second;
                if (!is_adjacent(_source, target, _g) && _source != target)
                    add_edge(_source, target, _g);
            }
            return *this;
        }

    private:
        Graph& _g;
        VertexMap _vertex_map;
        typename graph_traits<Graph>::vertex_descriptor _source;
    };

    template <class Graph, class Points, class PosMap>
    void operator()(Graph& g, Points& points, PosMap pos) const
    {
        typedef std::unordered_map
            <typename Triang::Vertex,
             typename graph_traits<Graph>::vertex_descriptor,
             hash_point> vertex_map_t;
        vertex_map_t vertex_map;

        Triang T;
        for (size_t i = 0; i < points.shape()[0]; ++i)
        {
            typename Triang::Point p(points[i][0], points[i][1], points[i][2]);
            typename graph_traits<Graph>::vertex_descriptor v = add_vertex(g);
            vertex_map[*T.insert(p)] = v;
            pos[v].resize(3);
            for (size_t j = 0; j < 3; ++j)
                pos[v][j] = points[i][j];
        }

        typename Triang::Vertex_iterator v_begin, v_end;
        v_begin = T.vertices_begin();
        v_end = T.vertices_end();
        while (v_begin != v_end)
        {
            typename Triang::Vertex_handle v = v_begin++;
            if (vertex_map.find(*v) == vertex_map.end())
                continue;
            edge_inserter<Graph, vertex_map_t> insert(g, *v, vertex_map);
            T.incident_vertices(v, insert);
        }
    }

};

} // namespace graph_tool

#endif // GRAPH_TRIANGULATION_HH
