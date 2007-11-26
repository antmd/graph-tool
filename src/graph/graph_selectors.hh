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

#ifndef GRAPH_SELECTORS_HH
#define GRAPH_SELECTORS_HH

#include <utility>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/map.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/find.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/dynamic_property_map.hpp>

#include "graph.hh"

namespace graph_tool
{

// This file contain selector for degree types for different types of
// graphs. Namely, they will return the in degree, out degree, and total degree
// for both directed and undirected graphs. There is also an scalar selector,
// which will return a specific scalar property

// This file also contains selectors for in_edge iterators of graphs, which
// return an empty range for undirected graphs

struct total_degreeS
{
    total_degreeS() {}
    total_degreeS(std::string scalar_property, const GraphInterface& g,
                  bool vertex = true) {}
    template <class Graph, class Vertex>
    size_t operator()(const Vertex& v, const Graph &g) const
    {
        using namespace boost;
        typedef typename is_convertible
            <typename graph_traits<Graph>::directed_category,
             directed_tag>::type is_directed;
        return get_total_degree(v,g,is_directed());
    }

    template <class Graph, class Vertex>
    size_t get_total_degree(const Vertex& v, const Graph &g, boost::true_type)
        const
    {
        return in_degree(v,g)+out_degree(v,g);
    }

    template <class Graph, class Vertex>
    size_t get_total_degree(const Vertex& v, const Graph &g, boost::false_type)
        const
    {
        return out_degree(v,g);
    }

    std::string name() {return "total_degree";}
};

struct in_degreeS
{
    in_degreeS() {}
    in_degreeS(std::string scalar_property, const GraphInterface& g,
               bool vertex = true) {}
    template <class Graph, class Vertex>
    size_t operator()(const Vertex& v, const Graph &g) const
    {
        using namespace boost;
        typedef typename is_convertible
            <typename graph_traits<Graph>::directed_category,
             directed_tag>::type is_directed;
        return get_in_degree(v,g,is_directed());
    }

    template <class Graph, class Vertex>
    size_t get_in_degree(const Vertex& v, const Graph &g, boost::true_type)
        const
    {
        return in_degree(v,g);
    }

    template <class Graph, class Vertex>
    size_t get_in_degree(const Vertex& v, const Graph &g, boost::false_type)
        const
    {
        return 0;
    }

    std::string name() {return "in_degree";}
};

struct out_degreeS
{
    out_degreeS() {}
    out_degreeS(std::string scalar_property, const GraphInterface& g,
                bool vertex = true) {}
    template <class Graph, class Vertex>
    size_t operator()(const Vertex& v, const Graph &g) const
    {
        return out_degree(v,g);
    }
    std::string name() {return "out_degree";}
};

struct scalarS
{
    scalarS(){}
    scalarS(std::string scalar_property, const GraphInterface& g,
            bool is_vertex_prop):
        _name(scalar_property), _g(&g)
    {
        if (is_vertex_prop)
            _scalar_property =
                &find_property_map(g._properties, _name,
                                   typeid(GraphInterface::vertex_t));
        else
            _scalar_property =
                &find_property_map(g._properties, _name,
                                   typeid(GraphInterface::edge_t));
    }

    typedef boost::mpl::vector<double,int,long double,float,long,unsigned long,
                               unsigned int,short,unsigned short,char,
                               unsigned char,bool,std::string> scalar_types;
    template <class Graph>
    double operator()
        (const typename boost::graph_traits<Graph>::vertex_descriptor& v,
         const Graph &g) const
    {
        using namespace boost;

        any value = _scalar_property->get(v);
        double *val = any_cast<double>(&value);
        if (val != 0)
            return *val;
        else
            return get_value<mpl::next<mpl::begin<scalar_types>::type>::type>
                (value);

    }

    template <class Graph>
    double operator()
        (const typename boost::graph_traits<Graph>::edge_descriptor& e,
         const Graph &g) const
    {
        using namespace boost;

        graph_traits<GraphInterface::multigraph_t>::edge_descriptor edge =  e;
        any value = _scalar_property->get(edge);
        double *val = any_cast<double>(&value);
        if (val != 0)
            return *val;
        else
            return get_value<mpl::next<mpl::begin<scalar_types>::type>::type>
                (value);
    }

    template <class ValueIter>
    double get_value(const boost::any& value , ValueIter = ValueIter()) const
    {
        using namespace boost;
        typedef typename mpl::deref<ValueIter>::type val_type;
        const val_type *val = any_cast<val_type>(&value);
        if (val != 0)
            return lexical_cast<double>(*val);
        else
            return get_value(value, typename mpl::next<ValueIter>::type());
    }

    double get_value(const boost::any& value,
                     boost::mpl::end<scalar_types>::type) const
    {
        throw boost::dynamic_get_failure(_name);
    }

    std::string name() {return _name;}
    std::string _name;
    boost::dynamic_property_map* _scalar_property;
    GraphInterface const* _g;
};


typedef boost::mpl::map
    < boost::mpl::pair<in_degreeS,boost::mpl::int_<GraphInterface::IN_DEGREE> >,
      boost::mpl::pair<out_degreeS,
                       boost::mpl::int_<GraphInterface::OUT_DEGREE> >,
      boost::mpl::pair<total_degreeS,
                       boost::mpl::int_<GraphInterface::TOTAL_DEGREE> >,
      boost::mpl::pair<scalarS,
                       boost::mpl::int_<GraphInterface::SCALAR> > >
    degree_selector_index;


// helper types for in_edge_iteratorS

template <class Graph, class IsDirected>
struct get_in_edges
{
    BOOST_MPL_ASSERT((is_same<IsDirected,boost::true_type>));
    BOOST_MPL_ASSERT((is_convertible
                      <typename graph_traits<Graph>::directed_category,
                      boost::directed_tag>));
    typedef typename graph_traits<Graph>::vertex_descriptor
        vertex_descriptor;
    typedef typename graph_traits<Graph>::in_edge_iterator type;
    static std::pair<type,type> in_edges(vertex_descriptor v,
                                         const Graph& g)
    {
        using namespace boost;
        return in_edges(v,g);
    }
};

template <class Graph>
struct get_in_edges<Graph,boost::false_type>
{
    BOOST_MPL_ASSERT((is_convertible
                      <typename graph_traits<Graph>::directed_category,
                      boost::undirected_tag>));
    typedef typename graph_traits<Graph>::vertex_descriptor
        vertex_descriptor;
    typedef typename graph_traits<Graph>::out_edge_iterator type;
    static std::pair<type,type> in_edges(vertex_descriptor v,
                                         const Graph& g)
    {
        return make_pair(type(), type());
    }
};

// this in-edge iterator selector returns the in-edge range for directed graphs
// and an empty out-edge range for undirected graphs. The iterator type is given
// by in_edge_iteratorS<Graph>::type.
template <class Graph>
struct in_edge_iteratorS
{
    typedef typename boost::graph_traits<Graph>::directed_category
        directed_category;
    typedef typename is_convertible<directed_category,
                                    boost::directed_tag>::type is_directed;
    typedef typename get_in_edges<Graph,is_directed>::type type;

    typedef typename boost::graph_traits<Graph>::vertex_descriptor
        vertex_descriptor;
    static std::pair<type,type> in_edges(vertex_descriptor v,
                                         const Graph& g)
    {
        return get_in_edges<Graph,is_directed>::in_edges(v, g);
    }
};

} //namespace graph_tool

#endif
