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

#ifndef GRAPH_WRAP_HH
#define GRAPH_WRAP_HH

#include <utility>
#include "graph_util.hh"
#include "graph_selectors.hh"

// Graph wrapper which takes care of edge index housekeeping

namespace boost
{
using namespace graph_tool;

template <class Graph>
class GraphWrap
{
  public:
    GraphWrap(Graph& g, GraphInterface& gi)
        : _g(g), _gi(gi) {}

    typedef typename Graph::vertex_property_type vertex_property_type;
    typedef typename Graph::edge_property_type edge_property_type;
    typedef typename Graph::graph_tag graph_tag;
    typedef Graph orig_graph_t;

    Graph& _g;
    GraphInterface& _gi;
};

template <class Graph>
GraphWrap<Graph> graph_wrap(Graph& g, GraphInterface& gi)
{
    return GraphWrap<Graph>(g, gi);
}

template <class Graph>
struct graph_traits<GraphWrap<Graph> >: public graph_traits<Graph> {};

template <class Graph>
inline typename graph_traits<GraphWrap<Graph> >::vertex_descriptor
source(typename graph_traits<GraphWrap<Graph> >::edge_descriptor e,
       const GraphWrap<Graph>& g)
{
    return source(e, g._g);
}

template <class Graph>
inline typename graph_traits<GraphWrap<Graph> >::vertex_descriptor
target(typename graph_traits<GraphWrap<Graph> >::edge_descriptor e,
       const GraphWrap<Graph>& g)
{
    return target(e, g._g);
}

template <class Graph>
inline typename graph_traits<GraphWrap<Graph> >::vertex_descriptor
vertex(typename graph_traits<GraphWrap<Graph> >::vertices_size_type n,
       const GraphWrap<Graph>& g)
{
    return vertex(n, g._g);
}

template <class Graph>
inline std::pair<typename graph_traits<GraphWrap<Graph> >::vertex_iterator,
                 typename graph_traits<GraphWrap<Graph> >::vertex_iterator>
vertices(const GraphWrap<Graph>& g)
{
    return vertices(g._g);
}

template <class Graph>
inline std::pair<typename graph_traits<GraphWrap<Graph> >::edge_iterator,
                 typename graph_traits<GraphWrap<Graph> >::edge_iterator>
edges(const GraphWrap<Graph>& g)
{
    return edges(g._g);
}

template <class Graph>
inline std::pair<typename graph_traits<GraphWrap<Graph> >::out_edge_iterator,
                 typename graph_traits<GraphWrap<Graph> >::out_edge_iterator >
out_edges(typename graph_traits<GraphWrap<Graph> >::vertex_descriptor u,
          const GraphWrap<Graph>& g)
{
    return out_edges(u, g._g);
}

template <class Graph>
inline std::pair<typename graph_traits<GraphWrap<Graph> >::in_edge_iterator,
                 typename graph_traits<GraphWrap<Graph> >::in_edge_iterator >
in_edges(typename graph_traits<GraphWrap<Graph> >::vertex_descriptor u,
         const GraphWrap<Graph>& g)
{
    return in_edges(u, g._g);
}

template <class Graph>
inline
std::pair<typename graph_traits<GraphWrap<Graph> >::adjacency_iterator,
          typename graph_traits<GraphWrap<Graph> >::adjacency_iterator>
adjacent_vertices
    (typename graph_traits<GraphWrap<Graph> >::vertex_descriptor u,
     const GraphWrap<Graph>& g)
{
    return adjacent_vertices(u, g._g);
}

template <class Graph>
inline typename graph_traits<GraphWrap<Graph> >::vertices_size_type
num_vertices(const GraphWrap<Graph>& g)
{
    return num_vertices(g._g);
}

template <class Graph>
inline typename graph_traits<GraphWrap<Graph> >::edges_size_type
num_edges(const GraphWrap<Graph>& g)
{
    return g._gi.GetNumberOfEdges();
}

template <class Graph>
inline typename graph_traits<GraphWrap<Graph> >::degree_size_type
out_degree(typename graph_traits<GraphWrap<Graph> >::vertex_descriptor u,
           const GraphWrap<Graph>& g)
{
    return out_degree(u, g._g);
}

template <class Graph>
inline typename graph_traits<GraphWrap<Graph> >::degree_size_type
in_degree(typename graph_traits<GraphWrap<Graph> >::vertex_descriptor u,
          const GraphWrap<Graph>& g)
{
    return in_degree(u, g._g);
}

template <class Graph>
inline typename graph_traits<GraphWrap<Graph> >::degree_size_type
degree(typename graph_traits<GraphWrap<Graph> >::vertex_descriptor u,
       const GraphWrap<Graph>& g)
{
    return degree(u, g._g);
}

template <class Graph>
inline std::pair<typename graph_traits<GraphWrap<Graph> >::edge_descriptor,
                 bool>
add_edge(typename graph_traits<GraphWrap<Graph> >::vertex_descriptor u,
         typename graph_traits<GraphWrap<Graph> >::vertex_descriptor v,
         GraphWrap<Graph>& g)
{
    std::pair<typename graph_traits<GraphWrap<Graph> >::edge_descriptor, bool>
        retval = add_edge(u, v, g._g);
    g._gi.AddEdgeIndex(retval.first);
    return retval;
}

template <class Graph>
inline void remove_edge
(typename graph_traits<GraphWrap<Graph> >::edge_descriptor e,
 GraphWrap<Graph>& g)
{
    g._gi.RemoveEdgeIndex(e);
}

template <class Graph>
inline void remove_edge
(typename graph_traits<GraphWrap<Graph> >::vertex_descriptor u,
 typename graph_traits<GraphWrap<Graph> >::vertex_descriptor v,
 Graph& g)
{
    vector<typename graph_traits<GraphWrap<Graph> >::edge_descriptor>
        removed_edges;

    typename graph_traits<GraphWrap<Graph> >::out_edge_iterator e, e_end;
    for(tie(e, e_end) = out_edges(u, g); e != e_end; ++e)
        if (target(*e, g) == v)
            removed_edges.push_back(*e);
    for (typeof(removed_edges.begin()) iter = removed_edges.begin();
         iter != removed_edges.end(); ++iter)
        remove_edge(*iter, g);
}


template <class Graph>
inline typename graph_traits<GraphWrap<Graph> >::vertex_descriptor
add_vertex(GraphWrap<Graph>& g)
{
    return add_vertex(g._g);
}

template <class Graph>
inline void clear_vertex
(typename graph_traits<GraphWrap<Graph> >::vertex_descriptor u,
 GraphWrap<Graph>& g)
{
    typedef GraphWrap<Graph> graph_t;
    vector<typename graph_traits<graph_t>::edge_descriptor> del_es;
    typename graph_traits<graph_t>::out_edge_iterator e, e_end;
    for (tie(e,e_end) == out_edges(u, g); e != e_end; ++e)
        del_es.push_back(*e);
    if (is_directed::apply<graph_t>::type::value)
    {
        typename in_edge_iteratorS<graph_t>::type e, e_end;
        for (tie(e,e_end) == in_edge_iteratorS<graph_t>::get_edges(u, g);
             e != e_end; ++e)
            del_es.push_back(*e);
    }
    for (size_t i = 0; i < del_es.size(); ++i)
        remove_edge(del_es[i], g);
}

template <class Graph>
inline void remove_vertex
(typename graph_traits<GraphWrap<Graph> >::vertex_descriptor u,
 GraphWrap<Graph>& g)
{
    clear_vertex(u, g);
    remove_vertex(u, g._g);
}

template <class Graph, class Predicate>
inline void
remove_out_edge_if
(typename graph_traits<GraphWrap<Graph> >::vertex_descriptor u,
 Predicate predicate, Graph& g)
{
    vector<typename graph_traits<GraphWrap<Graph> >::edge_descriptor>
        removed_edges;

    typename graph_traits<GraphWrap<Graph> >::out_edge_iterator e, e_end;
    for(tie(e, e_end) = out_edges(u, g); e != e_end; ++e)
        if (predicate(*e))
            removed_edges.push_back(*e);
    for (typeof(removed_edges.begin()) iter = removed_edges.begin();
         iter != removed_edges.end(); ++iter)
        remove_edge(*iter, g);
}

}

#endif // GRAPH_WRAP_HH
