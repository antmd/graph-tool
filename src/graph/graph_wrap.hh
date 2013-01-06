// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2013 Tiago de Paula Peixoto <tiago@skewed.de>
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

#include <boost/graph/filtered_graph.hpp>

#if (GCC_VERSION >= 40400)
#   include <tr1/unordered_map>
#else
#   include <boost/tr1/unordered_map.hpp>
#endif



// Graph wrapper which takes care of edge index housekeeping

namespace boost
{
using namespace graph_tool;

struct graph_wrap_tag {};

template <class Graph>
class GraphWrap
{
  public:
    GraphWrap(Graph& g, GraphInterface& gi)
        : _g(g), _gi(gi) {}

    typedef typename vertex_property_type<Graph>::type vertex_property_type;
    typedef typename edge_property_type<Graph>::type edge_property_type;
    typedef typename Graph::graph_tag orig_wrap_tag;
    typedef graph_wrap_tag graph_tag;
    typedef Graph orig_graph_t;

    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

    Graph& _g;
    GraphInterface& _gi;
};

namespace detail {

  template <typename PM>
  struct graph_wrap_edge_property_map {
    private:
    PM underlying_pm;

    public:
    typedef typename property_traits<PM>::key_type key_type;
    typedef typename property_traits<PM>::value_type value_type;
    typedef typename property_traits<PM>::reference reference;
    typedef typename property_traits<PM>::category category;

    explicit graph_wrap_edge_property_map(const PM& pm): underlying_pm(pm) {}

    friend reference
    get(const graph_wrap_edge_property_map& m,
        const key_type& e) {
      return get(m.underlying_pm, e.underlying_desc);
    }

    friend void
    put(const graph_wrap_edge_property_map& m,
        const key_type& e,
        const value_type& v) {
      put(m.underlying_pm, e, v);
    }

    reference operator[](const key_type& k) {
      return (this->underlying_pm)[k];
    }
  };

  struct graph_wrap_vertex_property_selector {
    template <class ReverseGraph, class Property, class Tag>
    struct bind_ {
      typedef typename ReverseGraph::orig_graph_t Graph;
      typedef property_map<Graph, Tag> PMap;
      typedef typename PMap::type type;
      typedef typename PMap::const_type const_type;
    };
  };

  struct graph_wrap_edge_property_selector {
    template <class ReverseGraph, class Property, class Tag>
    struct bind_ {
      typedef typename ReverseGraph::orig_graph_t Graph;
      typedef property_map<Graph, Tag> PMap;
      typedef graph_wrap_edge_property_map<typename PMap::type> type;
      typedef graph_wrap_edge_property_map<typename PMap::const_type> const_type;
    };
  };

} // namespace detail

template <>
struct vertex_property_selector<graph_wrap_tag> {
  typedef detail::graph_wrap_vertex_property_selector type;
};

template <>
struct edge_property_selector<graph_wrap_tag> {
  typedef detail::graph_wrap_edge_property_selector type;
};


template <class Graph>
GraphWrap<Graph> graph_wrap(Graph& g, GraphInterface& gi)
{
    return GraphWrap<Graph>(g, gi);
}

template <class Graph>
struct graph_traits<GraphWrap<Graph> >: public graph_traits<Graph> {};

template <class Graph>
struct graph_traits<const GraphWrap<Graph> >:
        public graph_traits<const Graph> {};

template <class Graph, class Tag>
struct property_map<GraphWrap<Graph>, Tag>:
        public property_map<Graph, Tag> {};

template <class Graph, class Tag>
typename property_map<Graph, Tag>::type
get(Tag t, const GraphWrap<Graph>& g)
{
    return get(t, g._g);
}

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
         GraphWrap<Graph> g)
{
    std::pair<typename graph_traits<GraphWrap<Graph> >::edge_descriptor, bool>
        retval = add_edge(u, v, g._g);
    g._gi.AddEdgeIndex(retval.first);
    return retval;
}

template <class Graph>
inline void remove_edge
(typename graph_traits<GraphWrap<Graph> >::edge_descriptor e,
 GraphWrap<Graph> g)
{
    g._gi.RemoveEdgeIndex(e);
}

template <class Graph>
inline void remove_edge(typename graph_traits<GraphWrap<Graph> >::vertex_descriptor u,
                        typename graph_traits<GraphWrap<Graph> >::vertex_descriptor v,
                        GraphWrap<Graph> g)
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
add_vertex(GraphWrap<Graph> g)
{
    return add_vertex(g._g);
}

template <class Graph>
inline void clear_vertex(typename graph_traits<GraphWrap<Graph> >::vertex_descriptor u,
                         GraphWrap<Graph> g)
{
    typedef GraphWrap<Graph> graph_t;
    GraphInterface::edge_index_map_t edge_index = g._gi.GetEdgeIndex();

    tr1::unordered_map<size_t, typename graph_traits<graph_t>::edge_descriptor> del_es;
    typename graph_traits<graph_t>::out_edge_iterator e, e_end;
    for (tie(e,e_end) = out_edges(u, g); e != e_end; ++e)
        del_es[edge_index[*e]] = *e;

    typename in_edge_iteratorS<graph_t>::type ie, ie_end;
    for (tie(ie,ie_end) = in_edge_iteratorS<graph_t>::get_edges(u, g);
         ie != ie_end; ++ie)
        del_es[edge_index[*ie]] = *ie;

    for (typeof(del_es.begin()) iter = del_es.begin(); iter != del_es.end();
         ++iter)
        remove_edge(iter->second, g);
}

// filtered graphs lack a remove_vertex function...
template <class Vertex, class Graph, class EdgePredicate, class VertexPredicate>
inline void remove_vertex(Vertex u, const filtered_graph<Graph,EdgePredicate,VertexPredicate>& g)
{
    remove_vertex(u, const_cast<Graph&>(g.m_g));
}

// reverse graphs lack a remove_vertex function...
template <class Vertex, class Graph>
inline void remove_vertex(Vertex u, const reverse_graph<Graph>& g)
{
    remove_vertex(u, const_cast<Graph&>(g.m_g));
}

template <class Graph>
inline void remove_vertex(typename graph_traits<GraphWrap<Graph> >::vertex_descriptor u,
                          GraphWrap<Graph> g)
{
    clear_vertex(u, g);
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    remove_vertex(vertex_t(u), g._g);
}


template <class Graph, class Predicate>
inline void
remove_out_edge_if(typename graph_traits<GraphWrap<Graph> >::vertex_descriptor u,
                   Predicate predicate, GraphWrap<Graph> g)
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
