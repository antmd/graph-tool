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

#ifndef GRAPH_ADAPTOR_HH
#define GRAPH_ADAPTOR_HH

#include <list>

#include <boost/config.hpp>
#include <boost/iterator_adaptors.hpp>
#include <boost/range/join.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include "transform_iterator.hh"

namespace boost {

//==============================================================================
// UndirectedAdaptor
// This class encapsulates a directed graph with parallel edges and provides a
// view of the graph as undirected with parallel edges.
// Encapsulated graph can be: VertexListGraph, EdgeListGraph, IncidenceGraph,
//                            AdjacencyGraph, VertexMutableGraph,
//                            EdgeMutableGraph, VertexMutablePropertyGraph,
//                            EdgeMutablePropertyGraph, BidirectionalGraph
// The undirected graph obeys the same concepts.
//==============================================================================
template <class Graph> class UndirectedAdaptor
{
public:
    UndirectedAdaptor(const Graph& g):_g(const_cast<Graph&>(g)){}

    typedef typename vertex_property_type<Graph>::type vertex_property_type;
    typedef typename edge_property_type<Graph>::type edge_property_type;
    typedef typename Graph::graph_tag graph_tag;
    typedef Graph graph_type;

    typedef Graph original_graph_t;

    typedef typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor
        vertex_descriptor;
    typedef typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor
        edge_descriptor;

    typedef undirected_tag directed_category;
    typedef allow_parallel_edge_tag edge_parallel_category;
    typedef typename graph_traits<Graph>::traversal_category traversal_category;

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

    static vertex_descriptor null_vertex() {graph_traits<Graph>::null_vertex();}

private:
    Graph& _g;
};

#ifndef BOOST_GRAPH_NO_BUNDLED_PROPERTIES
template<typename Graph>
struct vertex_bundle_type<UndirectedAdaptor<Graph> >:
        vertex_bundle_type<Graph> { };

template<typename Graph>
struct edge_bundle_type<UndirectedAdaptor<Graph> >:
        edge_bundle_type<Graph> { };
#endif // BOOST_GRAPH_NO_BUNDLED_PROPERTIES

template <class Graph>
struct get_iterator_category
{
    typedef typename graph_traits<Graph>::out_edge_iterator iter_t;
    typedef typename std::iterator_traits<iter_t>::iterator_category type;
};


template <class Graph, class Inverted = mpl::false_>
class joined_edge_iterator
    : public boost::iterator_facade<joined_edge_iterator<Graph, Inverted>,
                                    typename graph_traits<Graph>::edge_descriptor,
                                    typename get_iterator_category<Graph>::type,
                                    typename graph_traits<Graph>::edge_descriptor>
{
 public:
    typedef typename graph_traits<Graph>::in_edge_iterator in_iter_t;
    typedef typename graph_traits<Graph>::out_edge_iterator out_iter_t;

    typedef typename mpl::if_<Inverted, out_iter_t, in_iter_t>::type iter1_t;
    typedef typename mpl::if_<Inverted, in_iter_t, out_iter_t>::type iter2_t;

    joined_edge_iterator() {}
    explicit joined_edge_iterator(const std::pair<iter1_t, iter1_t>& range1,
                                  const std::pair<iter2_t, iter2_t>& range2,
                                  const iter1_t& pos1, const iter2_t& pos2)
        : _range1(range1), _range2(range2), _pos1(pos1), _pos2(pos2)
    { _flip = (_pos1 == _range1.second); }

 private:
    friend class boost::iterator_core_access;
    void increment()
    {
        if (_flip)
        {
            ++_pos2;
        }
        else
        {
            ++_pos1;
            _flip = (_pos1 == _range1.second);
        }
    }

    void decrement()
    {
        if (_flip)
        {
            if (_pos2 == _range2.first)
            {
                _flip = false;
                --_pos1;
            }
            else
            {
                --_pos2;
            }
        }
        else
        {
            --_pos1;
        }
    }

    typedef typename std::iterator_traits<iter1_t>::difference_type diff_t;
    void advance(diff_t n)
    {
        diff_t d1 = _range1.second - _pos1;
        if (n < d1)
        {
            _pos1 += n;
        }
        else
        {
            _pos1 = _range1.second;
            _pos2 += n - d1;
            _flip = true;
        }
    }

    diff_t distance_to(joined_edge_iterator const& other)
    {
        return (other._pos1 - _pos1) + (other._pos2 - _pos2);
    }

    bool equal(joined_edge_iterator const& other) const
    {
        return (_pos2 == other._pos2 && _pos1 == other._pos1);
    }

    typename graph_traits<Graph>::edge_descriptor dereference() const
    {
        typename graph_traits<Graph>::edge_descriptor e;
        if (_flip)
        {
            e = *_pos2;
        }
        else
        {
            e = *_pos1;
            e.inv = true;
        }
        return e;
    }

    std::pair<iter1_t, iter1_t> _range1;
    std::pair<iter2_t, iter2_t> _range2;
    iter1_t _pos1;
    iter2_t _pos2;
    bool _flip;
};



//==============================================================================
// UndirectedAdaptorAdjacencyIterator
// just keeps an internal reference to out_edge_iterator and calls target() when
// referenced
//==============================================================================

template <class Graph>
struct get_undirected_adjacency_iterator
{
    typedef joined_edge_iterator<Graph> out_edge_iter_t;
    typedef typename boost::adjacency_iterator_generator<UndirectedAdaptor<Graph>,
                                                         typename graph_traits<Graph>::vertex_descriptor,
                                                         out_edge_iter_t>::type type;
};


//==============================================================================
// graph_traits<UndirectedAdaptor>
// this defines all the necessary types associated with UndirectedAdaptor
//==============================================================================
template <class Graph>
struct graph_traits<UndirectedAdaptor<Graph> > {
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

    typedef typename get_undirected_adjacency_iterator<Graph>::type adjacency_iterator;
    typedef joined_edge_iterator<Graph, mpl::false_> out_edge_iterator;
    typedef joined_edge_iterator<Graph, mpl::true_> in_edge_iterator;
    typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
    typedef typename graph_traits<Graph>::edge_iterator edge_iterator;


    typedef undirected_tag directed_category;
    typedef allow_parallel_edge_tag edge_parallel_category;
    typedef typename graph_traits<Graph>::traversal_category traversal_category;
    typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;
    typedef typename graph_traits<Graph>::edges_size_type edges_size_type;
    typedef typename graph_traits<Graph>::degree_size_type degree_size_type;

    static vertex_descriptor null_vertex()
    {
        return graph_traits<Graph>::null_vertex();
    }

private:
    typedef is_convertible<typename std::iterator_traits<typename graph_traits<Graph>::out_edge_iterator>::iterator_category,
                           std::random_access_iterator_tag> is_orig_ra;
    typedef is_convertible<typename std::iterator_traits<out_edge_iterator>::iterator_category,
                           std::random_access_iterator_tag> is_ra;
    BOOST_STATIC_ASSERT((!is_orig_ra::value || is_ra::value));
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
source(const typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor& e,
       const UndirectedAdaptor<Graph>& g)
{
    if (e.inv)
        return target(e, g.OriginalGraph());
    else
        return source(e, g.OriginalGraph());
}

//==============================================================================
// target(e,g)
//==============================================================================
template <class Graph>
inline typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor
target(const typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor& e,
       const UndirectedAdaptor<Graph>& g)
{
    if (e.inv)
        return source(e, g.OriginalGraph());
    else
        return target(e, g.OriginalGraph());
}

//==============================================================================
// vertex(n,g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor
vertex(typename graph_traits<UndirectedAdaptor<Graph> >::vertices_size_type n,
       const UndirectedAdaptor<Graph>& g)
{
    return vertex(n, g.OriginalGraph());
}

//==============================================================================
// vertices(g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
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
inline __attribute__((always_inline))
std::pair<typename graph_traits<UndirectedAdaptor<Graph> >::edge_iterator,
          typename graph_traits<UndirectedAdaptor<Graph> >::edge_iterator >
edges(const UndirectedAdaptor<Graph>& g)
{
    return edges(g.OriginalGraph());
}

//==============================================================================
// edge(u, v, g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
std::pair<typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor,
          bool>
edge(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
     typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor v,
     const UndirectedAdaptor<Graph>& g)
{
    auto res = edge(u, v, g.OriginalGraph());

    if (!res.second)
    {
        res = edge(v, u, g.OriginalGraph());
        res.first.inv = true;
    }

    return res;
}

//==============================================================================
// out_edges(u,g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
std::pair<typename graph_traits<UndirectedAdaptor<Graph> >::out_edge_iterator,
          typename graph_traits<UndirectedAdaptor<Graph> >::out_edge_iterator >
out_edges(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
          const UndirectedAdaptor<Graph>& g)
{
    auto range1 = in_edges(u, g.OriginalGraph());
    auto range2 = out_edges(u, g.OriginalGraph());

    joined_edge_iterator<Graph> begin(range1, range2, range1.first, range2.first);
    joined_edge_iterator<Graph> end(range1, range2, range1.second, range2.second);

    return std::make_pair(begin, end);
}

//==============================================================================
// in_edges(u,g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
std::pair<typename graph_traits<UndirectedAdaptor<Graph> >::in_edge_iterator,
          typename graph_traits<UndirectedAdaptor<Graph> >::in_edge_iterator >
in_edges(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
         const UndirectedAdaptor<Graph>& g)
{
    auto range1 = out_edges(u, g.OriginalGraph());
    auto range2 = in_edges(u, g.OriginalGraph());

    joined_edge_iterator<Graph, mpl::true_> begin(range1, range2, range1.first, range2.first);
    joined_edge_iterator<Graph, mpl::true_> end(range1, range2, range1.second, range2.second);

    return std::make_pair(begin, end);
}

//==============================================================================
// adjacent_vertices(u,g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
std::pair<typename graph_traits<UndirectedAdaptor<Graph> >::adjacency_iterator,
          typename graph_traits<UndirectedAdaptor<Graph> >::adjacency_iterator>
adjacent_vertices
    (typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
     const UndirectedAdaptor<Graph>& g)
{
    std::pair<typename graph_traits<UndirectedAdaptor<Graph> >::out_edge_iterator,
              typename graph_traits<UndirectedAdaptor<Graph> >::out_edge_iterator>
        e_range = out_edges(u, g);

    typedef typename graph_traits<UndirectedAdaptor<Graph> >::adjacency_iterator adjacency_iterator;
    return std::make_pair(adjacency_iterator(e_range.first, &g),
                          adjacency_iterator(e_range.second, &g));
}

//==============================================================================
// num_vertices(g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
typename graph_traits<UndirectedAdaptor<Graph> >::vertices_size_type
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
inline __attribute__((always_inline))
typename graph_traits<UndirectedAdaptor<Graph> >::degree_size_type
out_degree(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
           const UndirectedAdaptor<Graph>& g)
{
    return (out_degree(u, g.OriginalGraph()) + in_degree(u, g.OriginalGraph()));
}

//==============================================================================
// in_degree(u,g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
typename graph_traits<UndirectedAdaptor<Graph> >::degree_size_type
in_degree(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
          const UndirectedAdaptor<Graph>& g)
{
    return out_degree(u, g);
}

//==============================================================================
// degree(u,g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
typename graph_traits<UndirectedAdaptor<Graph> >::degree_size_type
degree(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
       const UndirectedAdaptor<Graph>& g)
{
    return out_degree(u, g);
}


//==============================================================================
// add_vertex(g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor
add_vertex(UndirectedAdaptor<Graph>& g)
{
    return add_vertex(g.OriginalGraph());
}

//==============================================================================
// add_vertex(vp,g)
//==============================================================================
template <class Graph, class VertexProperties>
inline typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor
add_vertex(const VertexProperties& p, UndirectedAdaptor<Graph>& g)
{
    return add_vertex(p, g.OriginalGraph());
}

//==============================================================================
// clear_vertex(u,g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
void clear_vertex(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
                  UndirectedAdaptor<Graph>& g)
{
    clear_vertex(u, g.OriginalGraph());
}

//==============================================================================
// remove_vertex(u,g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
void remove_vertex(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
                   UndirectedAdaptor<Graph>& g)
{
    remove_vertex(u, g.OriginalGraph());
}

//==============================================================================
// remove_vertex_fast(u,g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
void remove_vertex_fast(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
                        UndirectedAdaptor<Graph>& g)
{
    remove_vertex_fast(u, g.OriginalGraph());
}

//==============================================================================
// add_edge(u,v,g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
std::pair<typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor,
          bool>
add_edge(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
         typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor v,
         UndirectedAdaptor<Graph>& g)
{
    return add_edge(u, v, g.OriginalGraph());
}

//==============================================================================
// add_edge(u,v,ep,g)
//==============================================================================
template <class Graph, class EdgeProperties>
inline __attribute__((always_inline))
std::pair<typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor,
          bool>
add_edge(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
         typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor v,
         const EdgeProperties& ep, UndirectedAdaptor<Graph>& g)
{
    return add_edge(u, v, ep, g.OriginalGraph());
}

//==============================================================================
// remove_edge(u,v,g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
void remove_edge(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
                 typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor v,
                 UndirectedAdaptor<Graph>& g)
{
    remove_edge(u, v, g.OriginalGraph());
    remove_edge(v, u, g.OriginalGraph());
}

//==============================================================================
// remove_edge(e,g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
void remove_edge (typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor e,
                  UndirectedAdaptor<Graph>& g)
{
    remove_edge(e, g.OriginalGraph());
}

//==============================================================================
// remove_edge(e_iter,g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
void remove_edge(typename graph_traits<UndirectedAdaptor<Graph> >::out_edge_iterator iter,
                 UndirectedAdaptor<Graph>& g)
{
    remove_edge(*iter, g);
}

//==============================================================================
// remove_out_edge_if(v,predicate,g)
//==============================================================================
template <class Graph, class Predicate>
inline __attribute__((always_inline))
void remove_out_edge_if(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor v,
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
// remove_in_edge_if(v,predicate,g)
//==============================================================================
template <class Graph, class Predicate>
inline __attribute__((always_inline))
void remove_in_edge_if(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor v,
                       Predicate predicate, UndirectedAdaptor<Graph>& g)
{
    std::list<typename UndirectedAdaptor<Graph>::EdgeDescriptor> removed_edges;
    typedef typename graph_traits<UndirectedAdaptor<Graph> >::in_edge_iterator
        iter_t;
    std::pair<iter_t, iter_t> edge_range;
    edge_range = in_edges(v,g);
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
    typedef typename property_map<Graph, PropertyTag>::const_type const_type;
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
inline __attribute__((always_inline))
typename property_map<UndirectedAdaptor<Graph>, PropertyTag>::type
get(PropertyTag tag, UndirectedAdaptor<Graph>& g)
{
    return get(tag, g.OriginalGraph());
}

//==============================================================================
// const get(tag,g)
//==============================================================================
template <class PropertyTag, class Graph>
inline __attribute__((always_inline))
typename property_map<UndirectedAdaptor<Graph>, PropertyTag>::const_type
get(PropertyTag tag, const UndirectedAdaptor<Graph>& g)
{
    return get(tag, g.OriginalGraph());
}

//==============================================================================
// get(tag,g,v)
//==============================================================================
template <class PropertyTag, class Graph>
inline __attribute__((always_inline))
typename property_traits
    <typename property_map<UndirectedAdaptor<Graph>,
                           PropertyTag>::const_type >::value_type
get(PropertyTag tag, const UndirectedAdaptor<Graph>& g,
    typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor v)
{
    return get(tag, g.OriginalGraph(), v);
}

//==============================================================================
// get(tag,g,e)
//==============================================================================
template <class PropertyTag, class Graph>
inline __attribute__((always_inline))
typename property_traits
    <typename property_map<UndirectedAdaptor<Graph>,
                           PropertyTag>::const_type >::value_type
get(PropertyTag tag, const UndirectedAdaptor<Graph>& g,
    typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor e)
{
    return get(tag, g.OriginalGraph(), e.OriginalEdge());
}

//==============================================================================
// put(tag, g, v, value)
//==============================================================================
template <class Graph, class PropertyTag, class Value>
inline __attribute__((always_inline))
void put(PropertyTag tag, UndirectedAdaptor<Graph>& g,
         typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor v,
         const Value& value)
{
    put(tag, g.OriginalGraph(), v, value);
}

//==============================================================================
// put(tag, g, e, value)
//==============================================================================
template <class Graph, class PropertyTag, class X, class Value>
inline __attribute__((always_inline))
void put(PropertyTag tag, const UndirectedAdaptor<Graph>& g,
         typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor e,
         const Value &value)
{
    put(tag, g.OriginalGraph(), e.OriginalEdge(), value);
}

//==============================================================================
// get_property(g,tag)
//==============================================================================
template <class Graph, class GraphProperties, class GraphPropertyTag>
inline __attribute__((always_inline))
typename property_value<GraphProperties, GraphPropertyTag>::type&
get_property(UndirectedAdaptor<Graph>& g, GraphPropertyTag tag)
{
    get_property(g.OriginalGraph(), tag);
}

//==============================================================================
// const get_property(g,tag)
//==============================================================================
template <class Graph, class GraphProperties, class GraphPropertyTag>
inline __attribute__((always_inline))
const typename property_value<GraphProperties, GraphPropertyTag>::type&
get_property(const UndirectedAdaptor<Graph>& g, GraphPropertyTag tag)
{
    get_property(g.OriginalGraph(), tag);
}

} // namespace boost


#endif // GRAPH_ADAPTOR_HH
