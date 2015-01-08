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

#ifndef GRAPH_UTIL_HH
#define GRAPH_UTIL_HH

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/filtered_graph.hpp>
#if (BOOST_VERSION / 100 % 1000 >= 48)
    #include <boost/graph/reverse_graph_alt.hpp>
#else
    #include <boost/graph/reverse_graph.hpp>
#endif
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <string>

#include <boost/functional/hash.hpp>

#include <functional>
#include <random>

#include "graph_selectors.hh"

namespace graph_tool
{

//
// Metaprogramming
// ===============

// useful metafunction to determine whether a graph is directed or not

struct is_directed
{
    template <class Graph>
    struct apply
    {
        typedef is_convertible<typename boost::graph_traits<Graph>::directed_category,
                               boost::directed_tag> type;
    };
};


// This will count "by hand" the number of vertices on a graph. Always O(V).
struct HardNumVertices
{
    template <class Graph>
    size_t operator()(Graph& g) const
    {
        size_t n = 0;
        typename boost::graph_traits<Graph>::vertex_iterator v_iter, v_begin, v_end;
        tie(v_begin, v_end) = vertices(g);
        for (v_iter = v_begin; v_iter != v_end; ++v_iter)
            n++;
        return n;
    }
};

// This will return the number of vertices on a graph, as given by
// num_vertices. Can be O(1).
struct SoftNumVertices
{
    template <class Graph>
    size_t operator()(Graph& g) const { return num_vertices(g); }
};

// This will count "by hand" the number of edges on a graph. Always O(E).
struct HardNumEdges
{
    template <class Graph>
    size_t operator()(Graph& g) const
    {
        size_t n = 0;
        typename boost::graph_traits<Graph>::edge_iterator e_iter, e_begin, e_end;
        tie(e_begin, e_end) = edges(g);
        for (e_iter = e_begin; e_iter != e_end; ++e_iter)
            n++;
        return n;
    }
};

// This will return the number of edges on a graph, as given by num_edges. Can
// be O(1).
struct SoftNumEdges
{
    template <class Graph>
    size_t operator()(Graph& g) const { return num_edges(g); }
};

// returns true if vertices u and v are adjacent. This is O(k(u)).
template <class Graph>
bool is_adjacent(typename boost::graph_traits<Graph>::vertex_descriptor u,
                 typename boost::graph_traits<Graph>::vertex_descriptor v,
                 const Graph& g )
{
    typename boost::graph_traits<Graph>::out_edge_iterator e, e_end;
    for (tie(e, e_end) = out_edges(u, g); e != e_end; ++e)
    {
        if (target(*e,g) == v)
            return true;
    }
    return false;
}

// computes the out-degree of a graph, ignoring self-edges
template <class Graph>
inline size_t
out_degree_no_loops(typename boost::graph_traits<Graph>::vertex_descriptor v,
                    const Graph &g)
{
    size_t k = 0;
    typename boost::graph_traits<Graph>::adjacency_iterator a,a_end;
    for (tie(a,a_end) = adjacent_vertices(v,g); a != a_end; ++a)
        if (*a != v)
            k++;
    return k;
}

// computes the out-degree of a graph, ignoring self-edges
template <class Graph, class Weights>
inline typename boost::property_traits<Weights>::value_type
out_degree_no_loops_weighted(typename boost::graph_traits<Graph>::vertex_descriptor v,
                             Weights w, const Graph &g)
{
    typename boost::property_traits<Weights>::value_type k = 0;
    typename boost::graph_traits<Graph>::out_edge_iterator e, e_end;
    for (tie(e, e_end) = out_edges(v, g); e != e_end; ++e)
        if (target(*e, g) != v)
            k += get(w, *e);
    return k;
}


template <class GraphOrig, class GraphTarget>
void graph_copy(const GraphOrig& g, GraphTarget& gt)
{
    typename boost::property_map<GraphOrig, boost::vertex_index_t>::type index = get(boost::vertex_index, g);
    typedef typename boost::graph_traits<GraphTarget>::vertex_descriptor tvertex_t;
    vector<tvertex_t> vmap(num_vertices(g));
    typename boost::graph_traits<GraphOrig>::vertex_iterator v, v_end;
    for (tie(v, v_end) = vertices(g); v != v_end; ++v)
        vmap[index[*v]] = add_vertex(gt);

    typename boost::graph_traits<GraphOrig>::edge_iterator e, e_end;
    for (tie(e, e_end) = edges(g); e != e_end; ++e)
        add_edge(vmap[index[source(*e, g)]], vmap[index[target(*e, g)]], gt);
}

} // namespace graph_tool

// some additional functions for filtered graphs, which don't exist by default
namespace boost
{
//==============================================================================
// vertex(i, filtered_graph<G>)
//==============================================================================
template <class Graph, class EdgePredicate, class VertexPredicate>
inline
typename boost::graph_traits<filtered_graph<Graph,
                                            EdgePredicate,
                                            VertexPredicate>>::vertex_descriptor
vertex(size_t i, const filtered_graph<Graph,EdgePredicate,VertexPredicate>& g)
{
    typename boost::graph_traits<Graph>::vertex_descriptor v = vertex(i, g.m_g);
    if (g.m_vertex_pred(v))
        return v;
    else
        return graph_traits<Graph>::null_vertex();
}

//==============================================================================
// vertex(i, reverse_graph<G>)
//==============================================================================
template <class Graph>
inline
typename boost::graph_traits<reverse_graph<Graph>>::vertex_descriptor
vertex(size_t i, const reverse_graph<Graph>& g)
{
    return vertex(i, g.m_g);
}

//==============================================================================
// add_edge(u, v, filtered_graph<G>)
//==============================================================================
template <class Graph, class EdgePredicate, class VertexPredicate>
inline
std::pair<typename boost::graph_traits
              <filtered_graph<Graph,EdgePredicate,
                              VertexPredicate>>::edge_descriptor, bool>
add_edge(typename boost::graph_traits
              <filtered_graph<Graph,EdgePredicate,
                              VertexPredicate>>::vertex_descriptor u,
         typename boost::graph_traits
              <filtered_graph<Graph,EdgePredicate,
                              VertexPredicate>>::vertex_descriptor v,
         filtered_graph<Graph,EdgePredicate,VertexPredicate>& g)
{
    return add_edge(u,v, const_cast<Graph&>(g.m_g));
}

//==============================================================================
// clear_vertex(v, filtered_graph<G>)
//==============================================================================
template <class Graph, class EdgePredicate, class VertexPredicate>
inline void
clear_vertex(typename boost::graph_traits
             <filtered_graph<Graph,EdgePredicate,VertexPredicate>>::vertex_descriptor v,
             filtered_graph<Graph,EdgePredicate,VertexPredicate>& g)
{
    typedef typename boost::graph_traits<filtered_graph<Graph,EdgePredicate,VertexPredicate>>::edge_descriptor e_t;
    std::vector<e_t> e_list;
    for (auto e : graph_tool::all_edges_range(v, g))
        e_list.push_back(e);
    for (auto& e : e_list)
        remove_edge(e, g);
}

//==============================================================================
// get(vertex_index_t, filtered_graph<G>)
//==============================================================================

// template <class Graph, class EdgePredicate, class VertexPredicate>
// typename property_map<Graph, vertex_index_t>::type
// get(vertex_index_t, const filtered_graph<Graph,EdgePredicate,
//                                          VertexPredicate>& g)
// {
//     return get(vertex_index_t(), g.m_g);
// }


//==============================================================================
// add_edge(u, v, reverse_graph<G>)
//==============================================================================
template <class Graph>
inline
std::pair<typename boost::graph_traits<reverse_graph<Graph>>::edge_descriptor,bool>
add_edge(typename boost::graph_traits<reverse_graph<Graph>>::vertex_descriptor u,
         typename boost::graph_traits<reverse_graph<Graph>>::vertex_descriptor v,
         reverse_graph<Graph>& g)
{
    typedef typename boost::graph_traits<reverse_graph<Graph>>::edge_descriptor e_t;
    std::pair<typename boost::graph_traits<Graph>::edge_descriptor,bool> ret =
        add_edge(v, u, const_cast<Graph&>(g.m_g)); // insert reversed
    return std::make_pair(e_t(ret.first), ret.second);
}

//==============================================================================
//remove_edge(e, filtered_graph<G>)
//==============================================================================
template <class Graph, class EdgePredicate, class VertexPredicate>
inline
void remove_edge(typename boost::graph_traits
                     <filtered_graph<Graph,EdgePredicate,
                                     VertexPredicate>>::edge_descriptor e,
                 filtered_graph<Graph,EdgePredicate,VertexPredicate>& g)
{
    return remove_edge(e,const_cast<Graph&>(g.m_g));
}

//==============================================================================
//remove_vertex(v, filtered_graph<G>)
//==============================================================================
template <class Graph, class EdgePredicate, class VertexPredicate>
inline
void remove_vertex(typename boost::graph_traits
                     <filtered_graph<Graph,EdgePredicate,
                                     VertexPredicate>>::vertex_descriptor v,
                 filtered_graph<Graph,EdgePredicate,VertexPredicate>& g)
{
    return remove_vertex(v,const_cast<Graph&>(g.m_g));
}

//==============================================================================
//remove_vertex_fast(v, filtered_graph<G>)
//==============================================================================
template <class Graph, class EdgePredicate, class VertexPredicate>
inline
void remove_vertex_fast(typename boost::graph_traits
                        <filtered_graph<Graph,EdgePredicate,
                                        VertexPredicate>>::vertex_descriptor v,
                        filtered_graph<Graph,EdgePredicate,VertexPredicate>& g)
{
    return remove_vertex_fast(v,const_cast<Graph&>(g.m_g));
}

//==============================================================================
//remove_edge(e, reverse_graph<G>)
//==============================================================================
template <class Graph>
inline
void remove_edge(typename boost::graph_traits<reverse_graph<Graph>>::edge_descriptor e,
                 reverse_graph<Graph>& g)
{
    return remove_edge(e,const_cast<Graph&>(g.m_g));
}

//==============================================================================
//remove_vertex(v, reverse_graph<G>)
//==============================================================================
template <class Graph>
inline
void remove_vertex(typename boost::graph_traits<reverse_graph<Graph>>::vertex_descriptor v,
                   reverse_graph<Graph>& g)
{
    return remove_vertex(v,const_cast<Graph&>(g.m_g));
}

//==============================================================================
//remove_vertex_fast(v, reverse_graph<G>)
//==============================================================================
template <class Graph>
inline
void remove_vertex_fast(typename boost::graph_traits<reverse_graph<Graph>>::vertex_descriptor v,
                        reverse_graph<Graph>& g)
{
    return remove_vertex_fast(v,const_cast<Graph&>(g.m_g));
}

//==============================================================================
//clear_vertex(v, reverse_graph<G>)
//==============================================================================
template <class Graph>
inline
void clear_vertex(typename boost::graph_traits<reverse_graph<Graph>>::vertex_descriptor v,
                  reverse_graph<Graph>& g)
{
    return clear_vertex(v,const_cast<Graph&>(g.m_g));
}

//==============================================================================
// add_vertex(filtered_graph<G>)
//==============================================================================
template <class Graph, class EdgePredicate, class VertexPredicate>
inline
typename boost::graph_traits
    <filtered_graph<Graph,EdgePredicate,VertexPredicate>>::vertex_descriptor
add_vertex(filtered_graph<Graph,EdgePredicate,VertexPredicate>& g)
{
    return add_vertex(const_cast<Graph&>(g.m_g));
}

//==============================================================================
// add_vertex(reverse_graph<G>)
//==============================================================================
template <class Graph>
inline
typename boost::graph_traits<reverse_graph<Graph>>::vertex_descriptor
add_vertex(reverse_graph<Graph>& g)
{
    return add_vertex(const_cast<Graph&>(g.m_g));
}

} // namespace boost


namespace std
{
// STL omission?
inline bool max(const bool& a, const bool& b) { return a || b; }
}

//
// Data type string representation
// ===============================
//
// String representation of individual data types. We have to take care
// specifically that no information is lost with floating point I/O.
//
// These are implemented in graph_io.cc.

namespace boost
{
using namespace std;

template <>
string lexical_cast<string,uint8_t>(const uint8_t& val);
template <>
uint8_t lexical_cast<uint8_t,string>(const string& val);
template <>
string lexical_cast<string,double>(const double& val);
template <>
double lexical_cast<double,string>(const string& val);
template <>
string lexical_cast<string,long double>(const long double& val);
template <>
long double lexical_cast<long double,string>(const string& val);
}

// std::vector<> stream i/o
namespace std
{
template <class Type>
ostream& operator<<(ostream& out, const vector<Type>& vec)
{
    for (size_t i = 0; i < vec.size(); ++i)
    {
        out << boost::lexical_cast<string>(vec[i]);
        if (i < vec.size() - 1)
            out << ", ";
    }
    return out;
}

template <class Type>
istream& operator>>(istream& in, vector<Type>& vec)
{
    using namespace boost;
    using namespace boost::algorithm;

    vec.clear();
    string data;
    getline(in, data);
    if (data == "")
        return in; // empty strings are OK
    vector<string> split_data;
    split(split_data, data, is_any_of(","));
    for (size_t i = 0; i < split_data.size(); ++i)
    {
        trim(split_data[i]);
        vec.push_back(lexical_cast<Type>(split_data[i]));
    }
    return in;
}

// string vectors need special attention, since separators must be properly
// escaped.
template <>
ostream& operator<<(ostream& out, const vector<string>& vec);

template <>
istream& operator>>(istream& in, vector<string>& vec);

} // std namespace

// This will iterate over a random permutation of a random access sequence, by
// swapping the values of the sequence as it iterates
template <class RandomAccessIterator, class RNG,
          class RandomDist = std::uniform_int_distribution<size_t>>
class random_permutation_iterator : public
    std::iterator<std::input_iterator_tag, typename RandomAccessIterator::value_type>
{
public:
    random_permutation_iterator(RandomAccessIterator begin,
                                RandomAccessIterator end, RNG& rng)
        : _i(begin), _end(end), _rng(&rng)
    {
        if(_i != _end)
        {
            RandomDist random(0,  _end - _i - 1);
            std::iter_swap(_i, _i + random(*_rng));
        }
    }

    typename RandomAccessIterator::value_type operator*()
    {
        return *_i;
    }

    random_permutation_iterator& operator++()
    {
        ++_i;
        if(_i != _end)
        {
            RandomDist random(0,  _end - _i - 1);
            std::iter_swap(_i, _i + random(*_rng));
        }
        return *this;
    }

    bool operator==(const random_permutation_iterator& ri)
    {
        return _i == ri._i;
    }

    bool operator!=(const random_permutation_iterator& ri)
    {
        return _i != ri._i;
    }

    size_t operator-(const random_permutation_iterator& ri)
    {
        return _i - ri._i;
    }

private:
    RandomAccessIterator _i, _end;
    RNG* _rng;
};


//
// Useful hash<> specializations
//

namespace std
{

template <class Val>
void _hash_combine(size_t& seed, const Val& hash)
{
    seed ^= std::hash<Val>()(hash) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

template <size_t pos, class... T>
struct tuple_combine
{
    void operator()(size_t& seed, const std::tuple<T...>& v) const
    {
        std::_hash_combine(seed, std::get<pos-1>(v));
        tuple_combine<pos-1, T...>()(seed, v);
    }
};

template <class... T>
struct tuple_combine<0, T...>
{
    void operator()(size_t& seed, const std::tuple<T...>& v) const {}
};

template <class... T>
struct hash<std::tuple<T...>>
{
    size_t operator()(std::tuple<T...> const& v) const
    {
        std::size_t seed = 0;
        tuple_combine<sizeof...(T), T...>()(seed, v);
        return seed;
    }
};

template <class T1, class T2>
struct hash<std::pair<T1, T2>>
{
    size_t operator()(std::pair<T1, T2> const& v) const
    {
        std::size_t seed = 0;
        std::_hash_combine(seed, v.first);
        std::_hash_combine(seed, v.second);
        return seed;
    }
};

template <class Value>
struct hash<vector<Value>>
{
    size_t operator()(const vector<Value>& v) const
    {
        size_t seed = 0;
        for (const auto& x : v)
            std::_hash_combine(seed, x);
        return seed;
    }
};

}


#endif // GRAPH_UTIL_HH
