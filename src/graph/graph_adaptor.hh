// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007-2010 Tiago de Paula Peixoto <tiago@skewed.de>
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

#ifndef GRAPH_ADAPTOR_HH
#define GRAPH_ADAPTOR_HH

#include <list>

#include <boost/config.hpp>
#include <boost/iterator_adaptors.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/filtered_graph.hpp>

namespace boost {

//==============================================================================
// UndirectedAdaptor
// This class encapsulates a directed graph with parallel edges and provides a
// view of the graph as undirected with parallel edges.
// Encapsulated graph must be: VertexListGraph, EdgeListGraph, IncidenceGraph,
//                             AdjacencyGraph, VertexMutableGraph,
//                             EdgeMutableGraph, VertexMutablePropertyGraph,
//                             EdgeMutablePropertyGraph, BidirectionalGraph
// The undirected graph obeys the same concepts.
//==============================================================================
template <class Graph> class UndirectedAdaptor
{
public:
    UndirectedAdaptor(const Graph &g):_g(const_cast<Graph &>(g)){}

    typedef typename vertex_property_type<Graph>::type vertex_property_type;
    typedef typename edge_property_type<Graph>::type edge_property_type;
    typedef typename Graph::graph_tag graph_tag;
    typedef Graph graph_type;

    class EdgeDescriptor;
    typedef Graph original_graph_t;

    typedef typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor
        vertex_descriptor;
    typedef typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor
        edge_descriptor;


#ifndef BOOST_GRAPH_NO_BUNDLED_PROPERTIES
    // Bundled properties support
    template<typename Descriptor>
    typename graph::detail::bundled_result<Graph, Descriptor>::type&
        operator[](Descriptor x) { return this->m_g[x]; }

    template<typename Descriptor>
    typename graph::detail::bundled_result<Graph, Descriptor>::type const&
        operator[](Descriptor x) const { return this->m_g[x]; }
#endif // BOOST_GRAPH_NO_BUNDLED_PROPERTIES


    const Graph& OriginalGraph() const {return _g;}
    Graph& OriginalGraph() {return _g;}

private:
    Graph &_g;
};

#ifndef BOOST_GRAPH_NO_BUNDLED_PROPERTIES
template<typename Graph>
struct vertex_bundle_type<UndirectedAdaptor<Graph> >:
        vertex_bundle_type<Graph> { };

template<typename Graph>
struct edge_bundle_type<UndirectedAdaptor<Graph> >:
        edge_bundle_type<Graph> { };
#endif // BOOST_GRAPH_NO_BUNDLED_PROPERTIES


//==============================================================================
// UndirectedAdaptor::EdgeDescriptor
//==============================================================================
template <class Graph>
class UndirectedAdaptor<Graph>::EdgeDescriptor:
        public graph_traits<Graph>::edge_descriptor
{
public:
    typedef typename graph_traits<Graph>::edge_descriptor original_edge_t;
    EdgeDescriptor(){}
    EdgeDescriptor(original_edge_t e): original_edge_t(e), _inverted(false) {}
    EdgeDescriptor(const original_edge_t &e,  bool inverted):
        original_edge_t(e), _inverted(inverted) {}

    bool IsInverted() const {return _inverted;}

    bool operator==(const EdgeDescriptor& e) const
    {
        return original_edge_t(e) == original_edge_t(*this);
    }

private:
    bool _inverted;
};

//==============================================================================
// UndirectedAdaptorEdgeIterator
//==============================================================================
template <typename Graph>
class UndirectedAdaptorEdgeIterator
    : public iterator<std::bidirectional_iterator_tag,
                      typename UndirectedAdaptor<Graph>::EdgeDescriptor,
                      std::ptrdiff_t,
                      typename UndirectedAdaptor<Graph>::EdgeDescriptor*,
                      typename UndirectedAdaptor<Graph>::EdgeDescriptor>
                      // not a reference!
{
public:
    UndirectedAdaptorEdgeIterator() {}
    explicit UndirectedAdaptorEdgeIterator
        (typename graph_traits<Graph>::edge_iterator &iter):_iter(iter){}
    typename UndirectedAdaptor<Graph>::EdgeDescriptor operator*() const
    {
        return (typename UndirectedAdaptor<Graph>::EdgeDescriptor(*_iter,
                                                                  false));
    }

//    pointer operator->() const {return **this;}

    UndirectedAdaptorEdgeIterator& operator++()
    {
        ++_iter;
        return *this;
    }

    UndirectedAdaptorEdgeIterator operator++(int)
    {
        UndirectedAdaptorEdgeIterator t = *this;
        ++_iter;
        return t;
    }

    UndirectedAdaptorEdgeIterator& operator--()
    {
        --_iter;
        return *this;
    }

    UndirectedAdaptorEdgeIterator operator--(int)
    {
        UndirectedAdaptorEdgeIterator t = *this;
        --_iter;
        return t;
    }

    bool operator==(UndirectedAdaptorEdgeIterator iter) const
    {
        return (_iter == iter._iter);
    }

    bool operator!=(UndirectedAdaptorEdgeIterator iter) const
    {
        return (_iter != iter._iter);
    }


private:
    typename graph_traits<Graph>::edge_iterator _iter;
};

//==============================================================================
// UndirectedAdaptorOutEdgeIterator
// this will iterate through both in_edges and out_edges of the underlying graph
//==============================================================================
template <typename Graph>
class UndirectedAdaptorOutEdgeIterator
    : public iterator<std::bidirectional_iterator_tag,
                      typename UndirectedAdaptor<Graph>::EdgeDescriptor,
                      std::ptrdiff_t,
                      typename UndirectedAdaptor<Graph>::EdgeDescriptor*,
                      typename UndirectedAdaptor<Graph>::EdgeDescriptor>
                      //not a reference
{
public:
    UndirectedAdaptorOutEdgeIterator() {};
    UndirectedAdaptorOutEdgeIterator
        (typename graph_traits<Graph>::out_edge_iterator out_iter,
         typename graph_traits<Graph>::in_edge_iterator in_iter,
         std::pair<typename graph_traits<Graph>::out_edge_iterator,
                   typename graph_traits<Graph>::out_edge_iterator> out_range,
         std::pair<typename graph_traits<Graph>::in_edge_iterator,
                   typename graph_traits<Graph>::in_edge_iterator> in_range)
            : _out_range(out_range), _in_range(in_range),
              _out_iter(out_iter), _in_iter(in_iter) {};

    typename UndirectedAdaptor<Graph>::EdgeDescriptor operator*() const
    {
        if ( _out_iter != _out_range.second )
            return (typename UndirectedAdaptor<Graph>::EdgeDescriptor
                    (*_out_iter,false));
        else
            return (typename UndirectedAdaptor<Graph>::EdgeDescriptor
                    (*_in_iter, true));
    }

//    pointer operator->() const {return **this;}

    UndirectedAdaptorOutEdgeIterator& operator++()
    {
        if (_out_iter != _out_range.second)
            ++_out_iter;
        else
            ++_in_iter;
        return *this;
    }

    UndirectedAdaptorOutEdgeIterator operator++(int)
    {
        UndirectedAdaptorOutEdgeIterator t = *this;
        if (_out_iter != _out_range.second)
            ++_out_iter;
        else
            ++_in_iter;
        return t;
    }

    UndirectedAdaptorOutEdgeIterator& operator--()
    {
        if (_in_iter == _in_range.first)
           --_out_iter;
        else
           --_in_iter;

        return *this;
    }

    UndirectedAdaptorOutEdgeIterator operator--(int)
    {
        UndirectedAdaptorOutEdgeIterator t = *this;
        if (_in_iter == _in_range.first)
           --_out_iter;
        else
           --_in_iter;

        return t;
    }

    bool operator==(UndirectedAdaptorOutEdgeIterator iter) const
    {
        return (_out_iter == iter._out_iter &&  _in_iter == iter._in_iter);
    }

    bool operator!=(UndirectedAdaptorOutEdgeIterator iter) const
    {
        return !(*this == iter);
    }

protected:
    std::pair<typename graph_traits<Graph>::out_edge_iterator,
              typename graph_traits<Graph>::out_edge_iterator> _out_range;
    std::pair<typename graph_traits<Graph>::in_edge_iterator,
              typename graph_traits<Graph>::in_edge_iterator> _in_range;
    typename graph_traits<Graph>::out_edge_iterator _out_iter;
    typename graph_traits<Graph>::in_edge_iterator _in_iter;
};

//==============================================================================
// UndirectedAdaptorAdjacencyIterator
// just keeps an internal reference to out_edge_iterator and calls target() when
// referenced
//==============================================================================
template <typename Graph>
class UndirectedAdaptorAdjacencyIterator
    : public iterator
        <std::bidirectional_iterator_tag,
         typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor,
         std::ptrdiff_t,
         typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor*,
         typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor>
         //not a reference
{
public:
    UndirectedAdaptorAdjacencyIterator(){};
    UndirectedAdaptorAdjacencyIterator
        (UndirectedAdaptorOutEdgeIterator<Graph> iter,
         const UndirectedAdaptor<Graph> &g)
            :_iter(iter), _g(&g) {}

    typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor operator*() const
    {
        return target(*_iter,*_g);
    }

//    pointer operator->() const {return **this;}

    UndirectedAdaptorAdjacencyIterator& operator++()
    {
        ++_iter;
        return *this;
    }

    UndirectedAdaptorAdjacencyIterator operator++(int)
    {
        UndirectedAdaptorAdjacencyIterator t = *this;
        ++_iter;
        return t;
    }

    UndirectedAdaptorAdjacencyIterator& operator--()
    {
        --_iter;
        return *this;
    }

    UndirectedAdaptorAdjacencyIterator operator--(int)
    {
        UndirectedAdaptorAdjacencyIterator t = *this;
        --_iter;
        return t;
    }

    bool operator==(UndirectedAdaptorAdjacencyIterator iter) const
    {
        return (iter._iter == _iter);
    }

    bool operator!=(UndirectedAdaptorAdjacencyIterator iter) const
    {
        return (iter._iter != _iter);
    }

private:
    UndirectedAdaptorOutEdgeIterator<Graph> _iter;
    UndirectedAdaptor<Graph> const * _g;
};

//==============================================================================
// undirected_adaptor_traversal_category
//==============================================================================
struct undirected_adaptor_traversal_category :
    public virtual bidirectional_graph_tag,
    public virtual adjacency_graph_tag,
    public virtual vertex_list_graph_tag,
    public virtual edge_list_graph_tag { };


//==============================================================================
// graph_traits<UndirectedAdaptor>
// this defines all the necessary types associated with UndirectedAdaptor
//==============================================================================
template <class Graph>
struct graph_traits<UndirectedAdaptor<Graph> > {
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef typename UndirectedAdaptor<Graph>::EdgeDescriptor edge_descriptor;

    typedef UndirectedAdaptorAdjacencyIterator<Graph> adjacency_iterator;
    typedef UndirectedAdaptorOutEdgeIterator<Graph> out_edge_iterator;
    typedef typename graph_traits<Graph>::in_edge_iterator in_edge_iterator;
    typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
    typedef UndirectedAdaptorEdgeIterator<Graph> edge_iterator;

    typedef undirected_tag directed_category;
    typedef allow_parallel_edge_tag edge_parallel_category;
    typedef undirected_adaptor_traversal_category traversal_category;
    typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;
    typedef typename graph_traits<Graph>::edges_size_type edges_size_type;
    typedef typename graph_traits<Graph>::degree_size_type degree_size_type;

    static vertex_descriptor null_vertex()
    {
        return graph_traits<Graph>::null_vertex();
    }
};

template <class Graph>
struct graph_traits< const UndirectedAdaptor<Graph> >:
    public graph_traits<UndirectedAdaptor<Graph> > {};

//==============================================================================
// Nonmember functions
// these provide manipulation of the graph
//==============================================================================

//==============================================================================
// source(e,g)
//==============================================================================
template <class Graph>
inline typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor
source(typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor e,
       const UndirectedAdaptor<Graph>& g)
{
    typedef typename graph_traits<Graph>::edge_descriptor original_edge_t;
    if (e.IsInverted())
        return target(original_edge_t(e), g.OriginalGraph());
    else
        return source(original_edge_t(e), g.OriginalGraph());
}

//==============================================================================
// target(e,g)
//==============================================================================
template <class Graph>
inline typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor
target(typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor e,
       const UndirectedAdaptor<Graph>& g)
{
    typedef typename graph_traits<Graph>::edge_descriptor original_edge_t;
    if (e.IsInverted())
        return source(original_edge_t(e), g.OriginalGraph());
    else
        return target(original_edge_t(e), g.OriginalGraph());
}

//==============================================================================
// vertex(n,g)
//==============================================================================
template <class Graph>
inline typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor
vertex(typename graph_traits<UndirectedAdaptor<Graph> >::vertices_size_type n,
       const UndirectedAdaptor<Graph>& g)
{
    return vertex(n, g.OriginalGraph());
}

//==============================================================================
// vertices(g)
//==============================================================================
template <class Graph>
inline
std::pair<typename graph_traits<UndirectedAdaptor<Graph> >::vertex_iterator,
          typename graph_traits<UndirectedAdaptor<Graph> >::vertex_iterator >
vertices(const UndirectedAdaptor<Graph>& g)
{
    return vertices(g.OriginalGraph());
}

//==============================================================================
// edges(g)
//==============================================================================
template <class Graph>
inline
std::pair<typename graph_traits<UndirectedAdaptor<Graph> >::edge_iterator,
          typename graph_traits<UndirectedAdaptor<Graph> >::edge_iterator >
edges(const UndirectedAdaptor<Graph>& g)
{
    std::pair<typename graph_traits<Graph>::edge_iterator,
              typename graph_traits<Graph>::edge_iterator> range;
    range = edges(g.OriginalGraph());
    return make_pair(UndirectedAdaptorEdgeIterator<Graph>(range.first),
                     UndirectedAdaptorEdgeIterator<Graph>(range.second));
}

//==============================================================================
// out_edges(u,g)
//==============================================================================
template <class Graph>
inline
std::pair<typename graph_traits<UndirectedAdaptor<Graph> >::out_edge_iterator,
          typename graph_traits<UndirectedAdaptor<Graph> >::out_edge_iterator >
out_edges(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
          const UndirectedAdaptor<Graph>& g)
{
    std::pair<typename graph_traits<Graph>::out_edge_iterator,
              typename graph_traits<Graph>::out_edge_iterator> out_range;
    std::pair<typename graph_traits<Graph>::in_edge_iterator,
              typename graph_traits<Graph>::in_edge_iterator> in_range;

    out_range = out_edges(u, g.OriginalGraph());
    in_range = in_edges(u, g.OriginalGraph());

    typedef typename graph_traits<UndirectedAdaptor<Graph> >::out_edge_iterator
        OutIter;

    OutIter iter_begin = OutIter(out_range.first, in_range.first, out_range,
                                 in_range);
    OutIter iter_end   = OutIter(out_range.second, in_range.second, out_range,
                                 in_range);

    return std::make_pair(iter_begin, iter_end);
}

//==============================================================================
// adjacent_vertices(u,g)
//==============================================================================
template <class Graph>
inline
std::pair<typename graph_traits<UndirectedAdaptor<Graph> >::adjacency_iterator,
          typename graph_traits<UndirectedAdaptor<Graph> >::adjacency_iterator>
adjacent_vertices
    (typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
     const UndirectedAdaptor<Graph>& g)
{
    typedef typename graph_traits<UndirectedAdaptor<Graph> >::out_edge_iterator
        edge_iter_t;
    typedef typename graph_traits<UndirectedAdaptor<Graph> >::adjacency_iterator
        adj_iter_t;

    std::pair<edge_iter_t, edge_iter_t> edge_range;
    edge_range = out_edges(u,g);

    return std::make_pair(adj_iter_t(edge_range.first, g),
                          adj_iter_t(edge_range.second, g));
}

//==============================================================================
// num_vertices(g)
//==============================================================================
template <class Graph>
inline typename graph_traits<UndirectedAdaptor<Graph> >::vertices_size_type
num_vertices(const UndirectedAdaptor<Graph>& g)
{
    return num_vertices(g.OriginalGraph());
}

//==============================================================================
// num_edges(g)
//==============================================================================
template <class Graph>
inline typename graph_traits<UndirectedAdaptor<Graph> >::edges_size_type
num_edges(const UndirectedAdaptor<Graph>& g)
{
    return num_edges(g.OriginalGraph());
}

//==============================================================================
// out_degree(u,g)
//==============================================================================
template <class Graph>
inline typename graph_traits<UndirectedAdaptor<Graph> >::degree_size_type
out_degree
    (typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
     const UndirectedAdaptor<Graph> &g)
{
    return (out_degree(u, g.OriginalGraph())+in_degree(u,g.OriginalGraph()));
}

//==============================================================================
// degree(u,g)
//==============================================================================
template <class Graph>
inline typename graph_traits<UndirectedAdaptor<Graph> >::degree_size_type
degree(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
       const UndirectedAdaptor<Graph> &g)
{
    return out_degree(u, g);
}


//==============================================================================
// add_vertex(g)
//==============================================================================
template <class Graph>
inline typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor
add_vertex(UndirectedAdaptor<Graph>& g)
{
    return add_vertex(g.OriginalGraph());
}

//==============================================================================
// add_vertex(vp,g)
//==============================================================================
template <class Graph, class VertexProperties>
inline typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor
add_vertex(const VertexProperties &p, UndirectedAdaptor<Graph>& g)
{
    return add_vertex(p, g.OriginalGraph());
}

//==============================================================================
// clear_vertex(u,g)
//==============================================================================
template <class Graph>
inline void clear_vertex
    (typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
     UndirectedAdaptor<Graph>& g)
{
    clear_vertex(u, g.OriginalGraph());
}

//==============================================================================
// remove_vertex(u,g)
//==============================================================================
template <class Graph>
inline void remove_vertex
    (typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
     UndirectedAdaptor<Graph>& g)
{
    remove_vertex(u, g.OriginalGraph());
}

//==============================================================================
// add_edge(u,v,g)
//==============================================================================
template <class Graph>
inline
std::pair<typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor,
          bool>
add_edge(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
         typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor v,
         UndirectedAdaptor<Graph>& g)
{
    std::pair<typename graph_traits<Graph>::edge_descriptor, bool> retval =
        add_edge(u,v,g.OriginalGraph());
    return std::make_pair
        (typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor
         (retval.first,false),
         retval.second);
}

//==============================================================================
// add_edge(u,v,ep,g)
//==============================================================================
template <class Graph, class EdgeProperties>
inline
std::pair<typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor,
          bool>
add_edge(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
         typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor v,
         const EdgeProperties& ep, UndirectedAdaptor<Graph>& g)
{
    std::pair<typename graph_traits<Graph>::edge_descriptor, bool> retval =
        add_edge(u,v,ep,g.OriginalGraph());
    return std::make_pair
        (typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor
         (retval.first,false),
         retval.second);
}

//==============================================================================
// remove_edge(u,v,g)
//==============================================================================
template <class Graph>
inline void remove_edge
    (typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
     typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor v,
     UndirectedAdaptor<Graph>& g)
{
    remove_edge(u,v,g.OriginalGraph());
    remove_edge(v,u,g.OriginalGraph());
}

//==============================================================================
// remove_edge(e,g)
//==============================================================================
template <class Graph>
inline void remove_edge
    (typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor e,
     UndirectedAdaptor<Graph>& g)
{
    remove_edge(typename graph_traits<Graph>::edge_descriptor(e),
                g.OriginalGraph());
}

//==============================================================================
// remove_edge(e_iter,g)
//==============================================================================
template <class Graph>
inline void remove_edge
    (typename graph_traits<UndirectedAdaptor<Graph> >::out_edge_iterator iter,
     UndirectedAdaptor<Graph>& g)
{
    remove_edge(*iter, g);
}

//==============================================================================
// remove_out_edge_if(v,predicate,g)
//==============================================================================
template <class Graph, class Predicate>
inline void remove_out_edge_if
    (typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor v,
     Predicate predicate, UndirectedAdaptor<Graph>& g)
{
    std::list<typename UndirectedAdaptor<Graph>::EdgeDescriptor> removed_edges;
    typedef typename graph_traits<UndirectedAdaptor<Graph> >::out_edge_iterator
        iter_t;
    std::pair<iter_t, iter_t> edge_range;
    edge_range = out_edges(v,g);
    for(iter_t iter = edge_range.first; iter != edge_range.second; ++iter)
        if (predicate(*iter))
            removed_edges.push_front(*iter);

    for(typeof(removed_edges.begin()) iter = removed_edges.begin();
        iter != removed_edges.end(); ++iter)
        remove_edge(*iter,g);
}


//==============================================================================
// Property maps
//==============================================================================

//==============================================================================
// vertex_property<UndirectedAdaptor>
//==============================================================================
template <class Graph>
class vertex_property<UndirectedAdaptor<Graph> >
{
public:
    typedef typename vertex_property<Graph>::type type;
};

//==============================================================================
// vertex_property_type<UndirectedAdaptor>
//==============================================================================
template <class Graph>
class vertex_property_type<UndirectedAdaptor<Graph> >
{
public:
    typedef typename vertex_property_type<Graph>::type type;
};

//==============================================================================
// edge_property<UndirectedAdaptor>
//==============================================================================
template <class Graph>
class edge_property<UndirectedAdaptor<Graph> >
{
public:
    typedef typename edge_property<Graph>::type type;
};

//==============================================================================
// edge_property_type<UndirectedAdaptor>
//==============================================================================
template <class Graph>
class edge_property_type<UndirectedAdaptor<Graph> >
{
public:
    typedef typename edge_property_type<Graph>::type type;
};

//==============================================================================
// property_map<UndirecterdAdaptor, PropertyTag>
//==============================================================================
template <class Graph, class PropertyTag>
class property_map<UndirectedAdaptor<Graph>, PropertyTag>
{
public:
    typedef typename property_map<Graph, PropertyTag>::type type;
    typedef typename property_map<const Graph, PropertyTag>::const_type
        const_type;
    //typedef typename property_map<Graph, PropertyTag>::const_type const type;
};

//==============================================================================
// property_map<UndirectedAdaptor, T Bundle::*>
//==============================================================================
template <typename Graph, typename T, typename Bundle>
class property_map<UndirectedAdaptor<Graph>, T Bundle::*>
{
public:
    typedef typename property_map<Graph, T Bundle::*>::type type;
    typedef typename property_map<Graph, T Bundle::*>::const_type const_type;
};


//==============================================================================
// get(tag,g)
//==============================================================================
template <class PropertyTag, class Graph>
inline typename property_map<UndirectedAdaptor<Graph>, PropertyTag>::type
get(PropertyTag tag, UndirectedAdaptor<Graph> &g)
{
    return get(tag, g.OriginalGraph());
}

//==============================================================================
// const get(tag,g)
//==============================================================================
template <class PropertyTag, class Graph>
inline typename property_map<UndirectedAdaptor<Graph>, PropertyTag>::const_type
get(PropertyTag tag, const UndirectedAdaptor<Graph> &g)
{
    return get(tag, g.OriginalGraph());
}

//==============================================================================
// get(tag,g,v)
//==============================================================================
template <class PropertyTag, class Graph>
inline
typename property_traits
    <typename property_map<UndirectedAdaptor<Graph>,
                           PropertyTag>::const_type >::value_type
get(PropertyTag tag, const UndirectedAdaptor<Graph> &g,
    typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor v)
{
    return get(tag, g.OriginalGraph(), v);
}

//==============================================================================
// get(tag,g,e)
//==============================================================================
template <class PropertyTag, class Graph>
inline
typename property_traits
    <typename property_map<UndirectedAdaptor<Graph>,
                           PropertyTag>::const_type >::value_type
get(PropertyTag tag, const UndirectedAdaptor<Graph> &g,
    typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor e)
{
    return get(tag, g.OriginalGraph(), e.OriginalEdge());
}

//==============================================================================
// put(tag, g, v, value)
//==============================================================================
template <class Graph, class PropertyTag, class Value>
inline
void put(PropertyTag tag, UndirectedAdaptor<Graph> &g,
         typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor v,
         const Value &value)
{
    put(tag, g.OriginalGraph(), v, value);
}

//==============================================================================
// put(tag, g, e, value)
//==============================================================================
template <class Graph, class PropertyTag, class X, class Value>
inline
void put(PropertyTag tag, const UndirectedAdaptor<Graph> &g,
         typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor e,
         const Value &value)
{
    put(tag, g.OriginalGraph(), e.OriginalEdge(), value);
}

//==============================================================================
// get_property(g,tag)
//==============================================================================
template <class Graph, class GraphProperties, class GraphPropertyTag>
inline
typename property_value<GraphProperties, GraphPropertyTag>::type&
get_property(UndirectedAdaptor<Graph> &g, GraphPropertyTag tag)
{
    get_property(g.OriginalGraph(), tag);
}

//==============================================================================
// const get_property(g,tag)
//==============================================================================
template <class Graph, class GraphProperties, class GraphPropertyTag>
inline
const typename property_value<GraphProperties, GraphPropertyTag>::type&
get_property(const UndirectedAdaptor<Graph> &g, GraphPropertyTag tag)
{
    get_property(g.OriginalGraph(), tag);
}

} // namespace boost


#endif // GRAPH_ADAPTOR_HH
