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

#ifndef FILTERING_HH
#define FILTERING_HH

#include "graph.hh"
#include <boost/version.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/filtered_graph.hpp>
#if (BOOST_VERSION / 100 % 1000 >= 48)
    #include <boost/graph/reverse_graph_alt.hpp>
#else
    #include <boost/graph/reverse_graph.hpp>
#endif
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
#include <boost/mpl/print.hpp>

#include "graph_adaptor.hh"
#include "graph_selectors.hh"
#include "graph_util.hh"
#include "mpl_nested_loop.hh"

#include <type_traits>

namespace graph_tool
{

// Graph filtering
// ---------------
//
// We want to generate versions of a template algorithm for every possible type
// of graph views. The types of graph views are the following:
//
//    - The original directed multigraph
//
//    - Filtered graphs, based on MaskFilter below
//
//    - A reversed view of each directed graph (original + filtered)
//
//    - An undirected view of each directed (unreversed) graph (original +
//      filtered)
//
// The total number of graph views is then: 1 + 1 + 2 + 2 = 6
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
    typedef typename boost::property_traits<DescriptorProperty>::value_type value_t;
    MaskFilter(){}
    MaskFilter(DescriptorProperty filtered_property, bool invert)
        : _filtered_property(filtered_property), _invert(invert) {}

    template <class Descriptor>
    __attribute__((always_inline)) inline bool operator() (Descriptor&& d) const
    {
        // ignore if masked

        return get(_filtered_property, std::forward<Descriptor>(d)) ^ _invert;

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
struct get_predicate<boost::keep_all>
{
    typedef boost::keep_all type;
};

// metafunction to get the filtered graph type
struct graph_filter
{
    template <class Graph, class EdgeProperty, class VertexProperty>
    struct apply
    {

        typedef typename get_predicate<EdgeProperty>::type edge_predicate;
        typedef typename get_predicate<VertexProperty>::type vertex_predicate;

        typedef boost::filtered_graph<Graph,
                                      edge_predicate,
                                      vertex_predicate> filtered_graph_t;

        // If both predicates are keep_all, then return the original graph
        // type. Otherwise return the filtered_graph type.
        typedef typename boost::mpl::if_<
            typename boost::mpl::and_<
                is_same<edge_predicate,
                        boost::keep_all>,
                is_same<vertex_predicate,
                        boost::keep_all>
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
        typedef boost::UndirectedAdaptor<Graph> type;
    };
};

// metafunction to get the reversed graph type
struct graph_reverse
{
    template <class Graph>
    struct apply
    {
        typedef boost::reverse_graph<Graph> type;
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
struct get_property_map_type<boost::keep_all, IndexMap>
{
    typedef boost::keep_all type;
};

// this metafunction returns a filtered graph type, given the scalar types to be
// used in the property maps
struct get_graph_filtered
{
    template <class Scalar>
    struct apply
    {
        // if the 'scalar' is the index map itself, use simply that, otherwise
        // get the specific property map
        typedef typename get_property_map_type<
            Scalar,
            GraphInterface::edge_index_map_t>::type edge_property_map;

        typedef typename get_property_map_type<
            Scalar,
            GraphInterface::vertex_index_map_t>::type vertex_property_map;

        typedef typename graph_filter::apply<GraphInterface::multigraph_t,
                                             edge_property_map,
                                             vertex_property_map>::type type;
    };
};

// this metafunction returns all the possible graph views
struct get_all_graph_views
{
    template <class FiltType,
              class AlwaysDirected = boost::mpl::bool_<false>,
              class NeverDirected = boost::mpl::bool_<false>,
              class AlwaysReversed = boost::mpl::bool_<false>,
              class NeverReversed = boost::mpl::bool_<false>,
              class NeverFiltered = boost::mpl::bool_<false> >
    struct apply
    {
        // filtered graphs
        struct filtered_graphs:
            boost::mpl::if_<NeverFiltered,
                            boost::mpl::vector<GraphInterface::multigraph_t>,
                            boost::mpl::vector<GraphInterface::multigraph_t,
                                               typename get_graph_filtered::apply<FiltType>::type>
                            >::type {};

        // filtered + reversed graphs
        struct reversed_graphs:
            boost::mpl::if_<AlwaysReversed,
                            typename boost::mpl::transform<filtered_graphs,
                                                           graph_reverse>::type,
                            typename boost::mpl::if_<
                                NeverReversed,
                                filtered_graphs,
                                typename boost::mpl::transform<
                                    filtered_graphs,
                                    graph_reverse,
                                    boost::mpl::back_inserter<filtered_graphs>
                                    >::type
                                >::type
                            >::type {};

        // undirected + filtereed + reversed graphs
        struct undirected_graphs:
            boost::mpl::if_<AlwaysDirected,
                            reversed_graphs,
                            typename boost::mpl::if_<
                                NeverDirected,
                                typename boost::mpl::transform<filtered_graphs,
                                                               graph_undirect>::type,
                                typename boost::mpl::transform<
                                    filtered_graphs,
                                    graph_undirect,
                                    boost::mpl::back_inserter<reversed_graphs>
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
            typedef typename boost::mpl::at<Sequence,Index>::type type;
        };
    };

    template <class Sequence>
    struct apply
    {
        typedef typename boost::mpl::size<Sequence>::type size;
        typedef typename boost::mpl::divides<size, boost::mpl::int_<2> >::type half_size;
        typedef typename boost::mpl::transform<boost::mpl::range_c<int, 0, half_size::value>,
                                               get_element<Sequence>,
                                               boost::mpl::back_inserter<boost::mpl::vector<> > >
            ::type first_part;
        typedef typename boost::mpl::transform<boost::mpl::range_c<int, half_size::value,
                                                                   size::value>,
                                               get_element<Sequence>,
                                               boost::mpl::back_inserter<boost::mpl::vector<> > >
            ::type second_part;
        typedef typename boost::mpl::pair<first_part,second_part> type;
    };
};



typedef uint8_t filt_scalar_type;

// finally, this type should hold all the possible graph views
struct all_graph_views:
    get_all_graph_views::apply<filt_scalar_type>::type {};

// restricted graph views
struct always_directed:
    get_all_graph_views::apply<filt_scalar_type,boost::mpl::bool_<true> >::type {};

struct never_directed:
    get_all_graph_views::apply<filt_scalar_type,boost::mpl::bool_<false>,
                               boost::mpl::bool_<true> >::type {};

struct always_reversed:
    get_all_graph_views::apply<filt_scalar_type,boost::mpl::bool_<true>,
                               boost::mpl::bool_<false>,boost::mpl::bool_<true> >::type {};

struct never_reversed:
    get_all_graph_views::apply<filt_scalar_type,boost::mpl::bool_<false>,
                               boost::mpl::bool_<false>,boost::mpl::bool_<false>,
                               boost::mpl::bool_<true> >::type {};

struct always_directed_never_reversed:
    get_all_graph_views::apply<filt_scalar_type,boost::mpl::bool_<true>,
                               boost::mpl::bool_<false>,boost::mpl::bool_<false>,
                               boost::mpl::bool_<true> >::type {};

struct never_filtered:
    get_all_graph_views::apply<filt_scalar_type,boost::mpl::bool_<false>,
                               boost::mpl::bool_<false>,boost::mpl::bool_<false>,
                               boost::mpl::bool_<false>,boost::mpl::bool_<true> >::type {};

// sanity check
typedef boost::mpl::size<all_graph_views>::type n_views;
#ifndef NO_GRAPH_FILTERING
BOOST_MPL_ASSERT_RELATION(n_views::value, == , boost::mpl::int_<6>::value);
#else
BOOST_MPL_ASSERT_RELATION(n_views::value, == , boost::mpl::int_<3>::value);
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
    boost::checked_vector_property_map<Type,GraphInterface::vertex_index_map_t>&
    uncheck(boost::checked_vector_property_map
            <Type,GraphInterface::vertex_index_map_t>& a, boost::mpl::true_) const
    {
        return a;
    }

    template <class Type>
    boost::unchecked_vector_property_map<Type,GraphInterface::vertex_index_map_t>
    uncheck(boost::checked_vector_property_map
            <Type,GraphInterface::vertex_index_map_t> a, boost::mpl::false_) const
    {
        return a.get_unchecked(_max_v);
    }

    template <class Type>
    boost::checked_vector_property_map<Type,GraphInterface::edge_index_map_t>&
    uncheck(boost::checked_vector_property_map
            <Type,GraphInterface::edge_index_map_t>& a, boost::mpl::true_) const
    {
        return a;
    }

    template <class Type>
    boost::unchecked_vector_property_map<Type,GraphInterface::edge_index_map_t>
    uncheck(boost::checked_vector_property_map
            <Type,GraphInterface::edge_index_map_t> a, boost::mpl::false_) const
    {
        return a.get_unchecked(_max_e);
    }

    template <class Type>
    scalarS<typename Type::unchecked_t>
    uncheck(scalarS<Type> a, boost::mpl::false_) const
    {
        return scalarS<typename Type::unchecked_t>(uncheck(a._pmap,
                                                           boost::mpl::false_()));
    }

    //no op
    template <class Type, class DoWrap>
    Type& uncheck(Type& a, DoWrap) const { return a; }

    void operator()() const {};

    template <class T1> void operator()(T1* a1) const
    { _a(*a1); }

    template <class T1, class... Ts>
    void operator()(T1* a1, Ts&&... as) const
    {
        _a(*a1, uncheck(std::forward<Ts>(as), Wrap())...);

    }

    Action _a;
    reference_wrapper<GraphInterface> _g;
    size_t _max_v, _max_e;
};

// this functor encapsulates another functor Action, which takes a pointer to a
// graph view as first argument
template <class Action, class GraphViews, class Wrap, class... TRS>
struct graph_action
{
    struct graph_view_pointers:
        boost::mpl::transform<GraphViews, boost::mpl::quote1<add_pointer> >::type {};

    graph_action(GraphInterface& g, Action a)
        : _g(g), _a(a, g, num_vertices(*g._mg),
                    max(g._mg->get_last_index(), size_t(1))) {}

    template <class... Args>
    void operator()(Args&&... args) const
    {
        bool found = false;
        boost::any gview = _g.GetGraphView();
        boost::mpl::nested_for_each<graph_view_pointers,TRS...>
            (boost::mpl::select_types(_a, found, gview, std::forward<Args>(args)...));
        if (!found)
        {
            vector<const std::type_info*> args_t = {(&(args).type())...};
            throw ActionNotFound(gview, typeid(Action), args_t);
        }
    }

    const GraphInterface &_g;
    action_wrap<Action, Wrap> _a;
};
} // details namespace


// all definitions of run_action with different arity
template <class GraphViews = detail::all_graph_views, class Wrap = boost::mpl::false_>
struct run_action
{
    template <class Action, class... TRS>
    detail::graph_action<Action,GraphViews,Wrap,TRS...>
    operator()(GraphInterface &g, Action&& a, TRS...)
    {
        return detail::graph_action<Action,GraphViews,Wrap,TRS...>(g, std::forward<Action>(a));
    }
};

// returns true if graph filtering was enabled at compile time
bool graph_filtering_enabled();

} //graph_tool namespace

#endif // FILTERING_HH
