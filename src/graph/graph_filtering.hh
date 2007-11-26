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

#ifndef FILTERING_HH
#define FILTERING_HH

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/reverse_graph.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/map.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/python/object.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/extract.hpp>

#include "graph_adaptor.hh"
#include "graph_selectors.hh"

// some additional functions for filtered graphs...
namespace boost
{
//==============================================================================
// vertex(i, filtered_graph<G>)
//==============================================================================
template <class Graph, class EdgePredicate, class VertexPredicate>
inline
typename graph_traits<filtered_graph<Graph,
                                     EdgePredicate,
                                     VertexPredicate> >::vertex_descriptor
vertex(size_t i, const filtered_graph<Graph,EdgePredicate,VertexPredicate>& g)
{
    typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g.m_g);
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
typename graph_traits<reverse_graph<Graph> >::vertex_descriptor
vertex(size_t i, const reverse_graph<Graph>& g)
{
    return vertex(i, g.m_g);
}


//==============================================================================
// add_edge(u, v, filtered_graph<G>)
//==============================================================================
template <class Graph, class EdgePredicate, class VertexPredicate>
inline
std::pair<typename graph_traits
              <filtered_graph<Graph,EdgePredicate,
                              VertexPredicate> >::edge_descriptor, bool>
add_edge(typename graph_traits
              <filtered_graph<Graph,EdgePredicate,
                              VertexPredicate> >::vertex_descriptor u,
         typename graph_traits
              <filtered_graph<Graph,EdgePredicate,
                              VertexPredicate> >::vertex_descriptor v,
         filtered_graph<Graph,EdgePredicate,VertexPredicate>& g)
{
    return add_edge(u,v, const_cast<Graph&>(g.m_g));
}

//==============================================================================
//remove_edge(e, filtered_graph<G>)
//==============================================================================
template <class Graph, class EdgePredicate, class VertexPredicate>
inline
void remove_edge(typename graph_traits
                     <filtered_graph<Graph,EdgePredicate,
                                     VertexPredicate> >::edge_descriptor e,
                 filtered_graph<Graph,EdgePredicate,VertexPredicate>& g)
{
    return remove_edge(e,const_cast<Graph&>(g.m_g));
}


} // namespace boost


namespace graph_tool
{
using namespace boost;

// this will count "by hand" the number of vertices on a graph
struct HardNumVertices
{
    template <class Graph>
    size_t operator()(const Graph &g) const
    {
        size_t n = 0;
        typename graph_traits<Graph>::vertex_iterator v_iter, v_begin, v_end;
        tie(v_begin, v_end) = vertices(g);
        for (v_iter = v_begin; v_iter != v_end; ++v_iter)
            n++;
        return n;
    }
};

// this will return the number of vertices on a graph, as given by num_vertices
struct SoftNumVertices
{
    template <class Graph>
    size_t operator()(const Graph &g) const { return num_vertices(g); }
};

// this will count "by hand" the number of edges on a graph
struct HardNumEdges
{
    template <class Graph>
    size_t operator()(const Graph &g) const
    {
        size_t n = 0;
        typename graph_traits<Graph>::edge_iterator e_iter, e_begin, e_end;
        tie(e_begin, e_end) = edges(g);
        for (e_iter = e_begin; e_iter != e_end; ++e_iter)
            n++;
        return n;
    }
};

// this will return the number of edges on a graph, as given by num_edges
struct SoftNumEdges
{
    template <class Graph>
    size_t operator()(const Graph &g) const { return num_edges(g); }
};

// The following filter predicate will filter out edges or vertices which lay
// outside a given range of values of a specific property. The range is
// specified as a pair of 'double' values

template <class FilteredPropertyMap>
class RangeFilter
{
public:
    RangeFilter(){}
    RangeFilter(const FilteredPropertyMap& filtered_property,
                const std::pair<double, double>& range,
                const std::pair<bool, bool>& include, const bool& invert)
        : _filtered_property(&filtered_property), _range(&range),
          _include(&include), _invert(&invert) {}

    template <class VertexOrEdge>
    bool operator() (VertexOrEdge e) const
    {
        bool retval;
        map_visitor<VertexOrEdge> visitor(e, *_range, *_include, *_invert, retval);
        apply_visitor(visitor, *_filtered_property);
        return retval;
    }

private:
    template <class VertexOrEdge>
    class map_visitor: public static_visitor<void>
    {
    public:
        map_visitor(const VertexOrEdge& descriptor,
                    const std::pair<double, double>& range,
                    const std::pair<bool, bool>& include, bool invert,
                    bool& retval)
            : _descriptor(descriptor), _range(range), _include(include),
              _invert(invert), _retval(retval) {}
        template <class MapType>
        void operator()(MapType& filter_prop)
        {
            // ignore if outside allowed range

            // TODO: This is a critical section. It will be called for every
            //       vertex or edge in the graph, every time they're iterated
            //       through. Therefore, it must be guaranteed this is as
            //       optimized as possible.

            double val = double(get(filter_prop, _descriptor));
            bool lower;
            if (_include.first)
                lower = val < _range.first;
            else
                lower = val <= _range.first;
            bool upper;
            if (_include.second)
                upper = val > _range.second;
            else
                upper = val >= _range.second;
            _retval = !(lower || upper);
            if (_invert)
                _retval = !_retval;
        }

    private:
        const VertexOrEdge& _descriptor;
        const std::pair<double, double>& _range;
        const std::pair<bool, bool>& _include;
        bool _invert;
        bool& _retval;
    };

    const FilteredPropertyMap* _filtered_property;
    const std::pair<double, double>* _range;
    const std::pair<bool, bool>* _include;
    const bool* _invert;
};


// below are some functions which will run a specific algorithm on a properly
// filtered graph, i.e., reversed, undirected and/or range filtered, or a
// combination of these

typedef mpl::vector<mpl::bool_<true>, mpl::bool_<false> > reverse_check;
typedef mpl::vector<mpl::bool_<false> > never_reversed;
typedef mpl::vector<mpl::bool_<true> > always_reversed;
typedef mpl::vector<mpl::bool_<true>, mpl::bool_<false> > directed_check;
typedef mpl::vector<mpl::bool_<true> > always_directed;
typedef mpl::vector<mpl::bool_<false> > always_undirected;

// this will check whether a graph is reversed and run the proper version of the
// algorithm

template <class Graph, class Action>
struct check_reverse
{
    check_reverse(const Graph &g, Action a, bool reverse, bool& found,
                  bool run_all)
        : _g(g), _a(a), _reverse(reverse), _found(found), _run_all(run_all) {}

    template <class Reverse>
    void operator()(Reverse) const
    {
        if (_reverse || _run_all)
        {
            static reverse_graph<Graph> rg(_g);
            _a(rg);
            _found = true;
        }
    }

    void operator()(mpl::bool_<false>) const
    {
        if (!_reverse || _run_all)
        {
            _a(_g);
            _found = true;
        }
    }

    const Graph &_g;
    Action _a;
    bool _reverse;
    bool& _found;
    bool _run_all;
};

// this will check whether a graph is directed and run the proper version of the
// algorithm

template <class Graph, class Action, class ReverseCheck>
struct check_directed
{
    check_directed(const Graph &g, Action a, bool reverse, bool directed,
                   bool& found, bool run_all)
        : _g(g), _a(a), _reverse(reverse), _directed(directed), _found(found),
          _run_all(run_all) {}

    template <class Directed>
    void operator()(Directed)
    {
        if (_directed || _run_all)
            mpl::for_each<ReverseCheck>
                (check_reverse<Graph, Action>(_g, _a, _reverse, _found,
                                              _run_all));
    }

    void operator()(mpl::bool_<false>)
    {
        if (!_directed || _run_all)
        {
            static UndirectedAdaptor<Graph> ug(_g);
            _a(ug);
            _found = true;
        }
    }

    const Graph &_g;
    Action _a;
    bool _reverse;
    bool _directed;
    bool& _found;
    bool _run_all;
};

// this will check whether a graph is range filtered and run the proper version
// of the algorithm

template <class Action, class ReverseCheck, class DirectedCheck>
void check_filter(const GraphInterface &g, Action a, ReverseCheck,
                  DirectedCheck, bool run_all = false)
{
    bool found = false;

    typedef RangeFilter<GraphInterface::vertex_filter_map_t> vertex_filter_t;
    typedef RangeFilter<GraphInterface::edge_filter_map_t> edge_filter_t;

#ifndef NO_RANGE_FILTERING
    if (g._vertex_filter_property != "" || g._edge_filter_property != "" ||
        run_all)
    {
        typedef filtered_graph<GraphInterface::multigraph_t, edge_filter_t,
            vertex_filter_t> fg_t;
        static fg_t fg(g._mg, edge_filter_t(g._edge_filter_map, g._edge_range,
                                            g._edge_range_include,
                                            g._edge_range_invert),
                       vertex_filter_t(g._vertex_filter_map, g._vertex_range,
                                       g._vertex_range_include,
                                       g._vertex_range_invert));
        mpl::for_each<DirectedCheck>
            (check_directed<fg_t,Action,ReverseCheck>(fg, a, g._reversed,
                                                      g._directed, found,
                                                      run_all));
    }
    if (!found || run_all)
    {
        mpl::for_each<DirectedCheck>
            (check_directed<GraphInterface::multigraph_t,Action,
             ReverseCheck>
             (g._mg, a, g._reversed, g._directed, found, run_all));
    }
#else
    mpl::for_each<DirectedCheck>
        (check_directed<GraphInterface::multigraph_t,Action,ReverseCheck>
         (g._mg, a, g._reversed, g._directed, found));
#endif
    if (!found)
        throw GraphException("graph filtering error: filter not found "
                             "(this is a bug in graph-tool)");
}

} //graph_tool namespace

#endif // FILTERING_HH
