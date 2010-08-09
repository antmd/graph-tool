// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007-2010 Tiago de Paula Peixoto <tiago@forked.de>
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
#include <boost/mpl/erase.hpp>
#include <boost/mpl/clear.hpp>
#include <boost/mpl/map.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/mpl/inserter.hpp>
#include <boost/mpl/insert_range.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/plus.hpp>
#include <boost/mpl/divides.hpp>
#include <boost/mpl/arithmetic.hpp>
#include <boost/mpl/greater_equal.hpp>
#include <boost/mpl/comparison.hpp>
#include <boost/mpl/transform_view.hpp>
#include <boost/mpl/quote.hpp>
#include <boost/mpl/range_c.hpp>

#include "graph.hh"
#include "graph_adaptor.hh"
#include "graph_selectors.hh"
#include "graph_util.hh"
#include "graph_wrap.hh"
#include "mpl_nested_loop.hh"

namespace graph_tool
{

using namespace boost;

// Graph filtering
// ---------------
//
// We want to generate versions of a template algorithm for every possible type
// of graph views. The types of graph views are the following:
//
//    - The original directed multigraph
//
//    - Filtered graphs, based on MaskFilter below. This amounts to a
//      filtered_graph for every combination of filtered and unfiltered vertex
//      or edge, i.e., 3.
//
//    - A reversed view of each directed graph (original + filtered)
//
//    - An undirected view of each directed (unreversed) graph (original +
//      filtered)
//
// The total number of graph views is then: 3 * 4 = 12
//
// The specific specialization can be called at run time (and generated at
// compile time) with the run_action() function, which takes as arguments the
// GraphInterface worked on, and the template functor to be specialized, which
// must take as first argument a pointer to a graph type, which itself must be a
// template parameter. Additionally, the function can be called optionally with
// up to 4 extra arguments, which must be MPL sequence instances corresponding
// to the type range of the extra arguments which must be passed to the template
// functor. The run_action() will not run the algorithm itself, but will instead
// return a functor (graph_action) which must be called either with no
// arguments, or with parameters of type boost::any which hold internally the
// type and values of the extra parameters to be passed to the action functor,
// and which will define the specialization to be chosen.
//
// Example:
//
// struct my_algorithm
// {
//     template <class Graph, class ValueType>
//     void operator()(Graph& g, ValueType val) const
//     {
//         ... do something ...
//     }
// };
//
// ...
//
// GraphInterface g;
// typedef mpl::vector<int, double, string> value_types;
// double foo = 42.0;
// run_action(g, my_algorithm(), value_types)(boost::any(foo));
//
// The above line will run my_algorithm::operator() with Graph being the
// appropriate graph view type and ValueType being 'double' and val = 42.0.

// Whenever no implementation is called, the following exception is thrown
class ActionNotFound: public GraphException
{
public:
    ActionNotFound(const boost::any& graph_view, const std::type_info& action,
                   const vector<const std::type_info*>& args);
    virtual const char * what () const throw ();
    virtual ~ActionNotFound() throw () {}
private:
    boost::any _graph_view;
    const std::type_info& _action;
    vector<const std::type_info*> _args;
};

namespace detail
{

// Implementation
// --------------
//
// The class MaskFilter below is the main filter predicate for the filtered
// graph view, based on descriptor property maps.  It filters out edges or
// vertices which are masked according to a property map with bool (actually
// uint8_t) value type.

template <class DescriptorProperty>
class MaskFilter
{
public:
    typedef typename property_traits<DescriptorProperty>::value_type value_t;
    MaskFilter(){}
    MaskFilter(DescriptorProperty filtered_property, bool invert)
        : _filtered_property(filtered_property), _invert(invert) {}

    template <class Descriptor>
    inline bool operator() (Descriptor d) const
    {
        // ignore if masked

        return get(_filtered_property, d) ^ _invert;

        // TODO: This is a critical section. It will be called for every vertex
        //       or edge in the graph, every time they're iterated
        //       through. Therefore, it must be guaranteed this is as optimized
        //       as possible.
    }

private:
    DescriptorProperty _filtered_property;
    bool _invert;
};


// Metaprogramming
// ---------------
//
// We need to generate a type sequence with all the filtered graph views, which
// will be called all_graph_views.

// metafunction class to get the correct filter predicate
template <class Property>
struct get_predicate
{
    typedef MaskFilter<Property> type;
};

template <>
struct get_predicate<keep_all>
{
    typedef keep_all type;
};

// metafunction to get the filtered graph type
struct graph_filter
{
    template <class Graph, class EdgeProperty, class VertexProperty>
    struct apply
    {

        typedef typename get_predicate<EdgeProperty>::type edge_predicate;
        typedef typename get_predicate<VertexProperty>::type vertex_predicate;

        typedef filtered_graph<Graph,
                               edge_predicate,
                               vertex_predicate> filtered_graph_t;

        // If both predicates are keep_all, then return the original graph
        // type. Otherwise return the filtered_graph type.
        typedef typename mpl::if_<
            typename mpl::and_<
                is_same<edge_predicate,
                        keep_all>,
                is_same<vertex_predicate,
                        keep_all>
                >::type,
            Graph,
            filtered_graph_t>::type type;
    };
};

// metafunction to get the undirected graph type
struct graph_undirect
{
    template <class Graph>
    struct apply
    {
        typedef UndirectedAdaptor<Graph> type;
    };
};

// metafunction to get the reversed graph type
struct graph_reverse
{
    template <class Graph>
    struct apply
    {
        typedef reverse_graph<Graph> type;
    };
};

// this wraps an unary metafunction class Bind into a unary metafunction,
// i.e., it is an identity operation. I'm not sure why it's necessary, but
// using pure unary bind expressions below didn't work for me, and this
// fixed it.
template <class Bind>
struct bind_wrap1
{
    template <class T1> struct apply
    { typedef typename Bind::template apply<T1>::type type; };
};

// metafunction which returns a mpl::vector containing all the pair combinations
// of two given type sequences
struct get_all_pairs
{
    struct make_pair
    {
        template <class T1, class T2>
        struct apply
        {
            typedef mpl::pair<T1,T2> type;
        };
    };

    template <class TR1, class TR2>
    struct apply
    {
        struct get_second_types
        {
            template <class T1>
            struct apply
            {
                typedef typename mpl::transform<
                    TR2,
                    bind_wrap1<
                        mpl::bind2<make_pair, T1, mpl::_1>
                    > >::type type;
            };
        };

        typedef typename mpl::transform<
            TR1,
            get_second_types,
            mpl::back_inserter<mpl::vector<> >
            >::type pair_combinations; // nested sequence (vector of vectors) of
                                       // pair combinations

        // joins two sequences
        struct join
        {
            template <class Seq1, class Seq2>
            struct apply
            {
                typedef typename boost::mpl::copy<
                    Seq2,
                    boost::mpl::back_inserter<Seq1>
                    >::type type;
            };
        };

        // flattens a nested sequence
        template<class Sequence>
        struct flatten
        {
            typedef typename boost::mpl::fold<
                Sequence,
                typename boost::mpl::clear<Sequence>::type,
                join
                >::type type;
        };

        // the complete list of combinations
        typedef typename flatten<pair_combinations>::type type;
    };
};

// metafunction class to get the correct property map type
template <class Scalar, class IndexMap>
struct get_property_map_type
{
    typedef typename property_map_type::apply<Scalar, IndexMap>
        ::type::unchecked_t type;
};

template <class IndexMap>
struct get_property_map_type<keep_all, IndexMap>
{
    typedef keep_all type;
};

// this metafunction returns a filtered graph type, given the scalar types to be
// used in the property maps
struct get_graph_filtered
{
    template <class TypePair>
    struct apply
    {
        typedef typename TypePair::first edge_scalar;
        typedef typename TypePair::second vertex_scalar;

        // if the 'scalar' is the index map itself, use simply that, otherwise
        // get the specific property map
        typedef typename mpl::if_<
            is_same<edge_scalar,
                    GraphInterface::edge_index_map_t>,
            GraphInterface::edge_index_map_t,
            typename get_property_map_type<
                edge_scalar,
                GraphInterface::edge_index_map_t>::type
            >::type edge_property_map;

        typedef typename mpl::if_<
            is_same<vertex_scalar,
                    GraphInterface::vertex_index_map_t>,
            GraphInterface::vertex_index_map_t,
            typename get_property_map_type<
                vertex_scalar,
                GraphInterface::vertex_index_map_t>::type
            >::type vertex_property_map;

        typedef typename graph_filter::apply<GraphInterface::multigraph_t,
                                             edge_property_map,
                                             vertex_property_map>::type type;
    };
};

// this metafunction returns all the possible graph views
struct get_all_graph_views
{
    template <class TypePairs,
              class AlwaysDirected = mpl::bool_<false>,
              class NeverDirected = mpl::bool_<false>,
              class AlwaysReversed = mpl::bool_<false>,
              class NeverReversed = mpl::bool_<false>,
              class NeverFiltered = mpl::bool_<false> >
    struct apply
    {
        // filtered graphs
        struct filtered_graphs:
            mpl::if_
            <NeverFiltered,
             mpl::vector<GraphInterface::multigraph_t>,
             typename mpl::transform<TypePairs,
                                     get_graph_filtered>::type>::type {};

        // filtered + reversed graphs
        struct reversed_graphs:
            mpl::if_<AlwaysReversed,
                     typename mpl::transform<filtered_graphs,
                                             graph_reverse>::type,
                     typename mpl::if_<
                         NeverReversed,
                         filtered_graphs,
                         typename mpl::transform<
                             filtered_graphs,
                             graph_reverse,
                             mpl::back_inserter<filtered_graphs>
                             >::type
                         >::type
                >::type {};

        // undirected + filtereed + reversed graphs
        struct undirected_graphs:
            mpl::if_<AlwaysDirected,
                     reversed_graphs,
                     typename mpl::if_<
                         NeverDirected,
                         typename mpl::transform<filtered_graphs,
                                                 graph_undirect>::type,
                         typename mpl::transform<
                             filtered_graphs,
                             graph_undirect,
                             mpl::back_inserter<reversed_graphs>
                             >::type
                         >::type
                     >::type {};
        typedef undirected_graphs type;
    };
};

// useful metafunction to split sequences in half
struct split
{
    template <class Sequence>
    struct get_element
    {
        template <class Index>
        struct apply
        {
            typedef typename mpl::at<Sequence,Index>::type type;
        };
    };

    template <class Sequence>
    struct apply
    {
        typedef typename mpl::size<Sequence>::type size;
        typedef typename mpl::divides<size, mpl::int_<2> >::type half_size;
        typedef typename mpl::transform<mpl::range_c<int, 0, half_size::value>,
                                        get_element<Sequence>,
                                        mpl::back_inserter<mpl::vector<> > >
            ::type first_part;
        typedef typename mpl::transform<mpl::range_c<int, half_size::value,
                                                     size::value>,
                                        get_element<Sequence>,
                                        mpl::back_inserter<mpl::vector<> > >
            ::type second_part;
        typedef typename mpl::pair<first_part,second_part> type;
    };
};



// all scalar types plus edge and vertex index property (we actually only use
// bool)

#ifndef NO_GRAPH_FILTERING
struct edge_scalars:
    mpl::vector<keep_all, uint8_t> {};
struct vertex_scalars:
    mpl::vector<keep_all, uint8_t> {};
#else
struct edge_scalars:
    mpl::vector<keep_all> {};
struct vertex_scalars:
    mpl::vector<keep_all> {};
#endif

// all scalar pairs
struct scalar_pairs: get_all_pairs::apply<edge_scalars,vertex_scalars>::type {};

// finally, this type should hold all the possible graph views
struct all_graph_views:
    get_all_graph_views::apply<scalar_pairs>::type {};

// restricted graph views
struct always_directed:
    get_all_graph_views::apply<scalar_pairs,mpl::bool_<true> >::type {};

struct never_directed:
    get_all_graph_views::apply<scalar_pairs,mpl::bool_<false>,
                               mpl::bool_<true> >::type {};

struct always_reversed:
    get_all_graph_views::apply<scalar_pairs,mpl::bool_<true>,
                               mpl::bool_<false>,mpl::bool_<true> >::type {};

struct never_reversed:
    get_all_graph_views::apply<scalar_pairs,mpl::bool_<false>,
                               mpl::bool_<false>,mpl::bool_<false>,
                               mpl::bool_<true> >::type {};

struct always_directed_never_reversed:
    get_all_graph_views::apply<scalar_pairs,mpl::bool_<true>,
                               mpl::bool_<false>,mpl::bool_<false>,
                               mpl::bool_<true> >::type {};

struct never_filtered:
    get_all_graph_views::apply<scalar_pairs,mpl::bool_<false>,
                               mpl::bool_<false>,mpl::bool_<false>,
                               mpl::bool_<false>,mpl::bool_<true> >::type {};

// sanity check
typedef mpl::size<all_graph_views>::type n_views;
#ifndef NO_GRAPH_FILTERING
BOOST_MPL_ASSERT_RELATION(n_views::value, == , mpl::int_<12>::value);
#else
BOOST_MPL_ASSERT_RELATION(n_views::value, == , mpl::int_<3>::value);
#endif

// run_action() implementation
// ===========================

// wrap action to be called, to deal with property maps, i.e., return version
// with no bounds checking.
template <class Action, class Wrap>
struct action_wrap
{
    action_wrap(Action a, GraphInterface& g, size_t max_v, size_t max_e)
        : _a(a), _g(g), _max_v(max_v), _max_e(max_e) {}

    template <class Type>
    typename checked_vector_property_map
        <Type,GraphInterface::vertex_index_map_t>::unchecked_t
    uncheck(checked_vector_property_map
            <Type,GraphInterface::vertex_index_map_t> a) const
    {
        return a.get_unchecked(_max_v);
    }

    template <class Type>
    typename checked_vector_property_map
        <Type,GraphInterface::edge_index_map_t>::unchecked_t
    uncheck(checked_vector_property_map
            <Type,GraphInterface::edge_index_map_t> a) const
    {
        return a.get_unchecked(_max_e);
    }

    template <class Type>
    scalarS<typename Type::unchecked_t>
    uncheck(scalarS<Type> a) const
    {
        return scalarS<typename Type::unchecked_t>(uncheck(a._pmap));
    }

    //no op
    template <class Type>
    Type uncheck(Type a) const { return a; }

    template <class Graph>
    GraphWrap<Graph> wrap(Graph* g, mpl::true_) const
    {
        return graph_wrap(*g, _g);
    }

    template <class Graph>
    Graph& wrap(Graph* g, mpl::false_) const
    {
        return *g;
    }

    void operator()() const {};
    template <class T1> void operator()(const T1& a1) const { _a(wrap(a1, Wrap())); }
    template <class T1, class T2>
    void operator()(const T1& a1, const T2& a2) const
    { _a(wrap(a1,Wrap()), uncheck(a2)); }
    template <class T1, class T2, class T3>
    void operator()(const T1& a1, const T2& a2, const T3& a3) const
    { _a(wrap(a1,Wrap()), uncheck(a2), uncheck(a3));}
    template <class T1, class T2, class T3, class T4>
    void operator()(const T1& a1, const T2& a2, const T3& a3, const T4& a4)
        const
    { _a(wrap(a1,Wrap()), uncheck(a2), uncheck(a3), uncheck(a4)); }
    template <class T1, class T2, class T3, class T4, class T5>
    void operator()(const T1& a1, const T2& a2, const T3& a3, const T4& a4,
                    const T5& a5) const
    { _a(wrap(a1,Wrap()), uncheck(a2), uncheck(a3), uncheck(a4), uncheck(a5)); }

    Action _a;
    reference_wrapper<GraphInterface> _g;
    size_t _max_v, _max_e;
};

// this functor encapsulates another functor Action, which takes a pointer to a
// graph view as first argument
template <class Action, class GraphViews, class Wrap, class TR1, class TR2,
          class TR3, class TR4 >
struct graph_action
{
    struct graph_view_pointers:
        mpl::transform<GraphViews, mpl::quote1<add_pointer> >::type {};

    graph_action(GraphInterface& g, Action a)
        : _g(g), _a(a, g, num_vertices(g._mg), g._max_edge_index + 1) {}

    void operator()() const
    {
        bool found = false;
        boost::any gview = _g.GetGraphView();
        boost::mpl::for_each<graph_view_pointers>
            (boost::mpl::select_types(_a, found, gview));
        if (!found)
        {
            throw ActionNotFound(gview, typeid(Action),
                                 vector<const std::type_info*>());
        }
    }

    void operator()(boost::any a1) const
    {
        bool found = false;
        boost::any gview = _g.GetGraphView();
        boost::mpl::nested_for_each<graph_view_pointers,TR1>()
            (boost::mpl::select_types(_a, found, gview, a1));
        if (!found)
        {
            vector<const std::type_info*> args;
            args.push_back(&a1.type());
            throw ActionNotFound(gview, typeid(Action), args);
        }
    }

    void operator()(boost::any a1, boost::any a2) const
    {
        bool found = false;
        boost::any gview = _g.GetGraphView();
        boost::mpl::nested_for_each<graph_view_pointers,TR1,TR2>()
            (boost::mpl::select_types(_a, found, gview, a1, a2));
        if (!found)
        {
            vector<const std::type_info*> args;
            args.push_back(&a1.type());
            args.push_back(&a2.type());
            throw ActionNotFound(gview, typeid(Action), args);
        }
    }

    void operator()(boost::any a1, boost::any a2, boost::any a3) const
    {
        bool found = false;
        boost::any gview = _g.GetGraphView();
        boost::mpl::nested_for_each<graph_view_pointers,TR1,TR2,TR3>()
            (boost::mpl::select_types(_a, found, gview, a1, a2, a3));
        if (!found)
        {
            vector<const std::type_info*> args;
            args.push_back(&a1.type());
            args.push_back(&a2.type());
            args.push_back(&a3.type());
            throw ActionNotFound(gview, typeid(Action), args);
        }
    }

    void operator()(boost::any a1, boost::any a2, boost::any a3,
                    boost::any a4) const
    {
        bool found = false;
        boost::any gview = _g.GetGraphView();
        boost::mpl::nested_for_each<graph_view_pointers,TR1,TR2,TR3,TR4>()
            (boost::mpl::select_types(_a, found, gview, a1, a2, a3,a4));
        if (!found)
        {
            vector<const std::type_info*> args;
            args.push_back(&a1.type());
            args.push_back(&a2.type());
            args.push_back(&a3.type());
            args.push_back(&a4.type());
            throw ActionNotFound(gview, typeid(Action), args);
        }
    }

    const GraphInterface &_g;
    action_wrap<Action, Wrap> _a;
};
} // details namespace


// all definitions of run_action with different arity
template <class GraphViews = detail::all_graph_views, class Wrap = mpl::false_>
struct run_action
{
    template <class Action>
    detail::graph_action<Action,GraphViews,Wrap>
    operator()(GraphInterface &g, Action a)
    {
        return detail::graph_action<Action,GraphViews,Wrap>(g, a);
    }

    template <class Action, class TR1>
    detail::graph_action<Action,GraphViews,Wrap,TR1>
    operator()(GraphInterface &g, Action a, TR1)
    {
        return detail::graph_action<Action,GraphViews,Wrap,TR1>(g, a);
    }

    template <class Action, class TR1, class TR2>
    detail::graph_action<Action,GraphViews,Wrap,TR1,TR2>
    operator()(GraphInterface &g, Action a, TR1, TR2)
    {
        return detail::graph_action<Action,GraphViews,Wrap,TR1,TR2>(g, a);
    }

    template <class Action, class TR1, class TR2, class TR3>
    detail::graph_action<Action,GraphViews,Wrap,TR1,TR2,TR3>
    operator()(GraphInterface &g, Action a, TR1, TR2, TR3)
    {
        return detail::graph_action<Action,GraphViews,Wrap,TR1,TR2,TR3>(g, a);
    }

    template <class Action, class TR1, class TR2, class TR3, class TR4>
    detail::graph_action<Action,GraphViews,Wrap,TR1,TR2,TR3,TR4>
    operator()(GraphInterface &g, Action a, TR1, TR2, TR3, TR4)
    {
        return detail::graph_action<Action,GraphViews,Wrap,TR1,TR2,TR3,TR4>(g, a);
    }
};

// returns true if graph filtering was enabled at compile time
bool graph_filtering_enabled();

} //graph_tool namespace

#endif // FILTERING_HH
