// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006  Tiago de Paula Peixoto <tiago@forked.de>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

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

#include "graph_adaptor.hh"

namespace graph_tool
{

struct total_degreeS 
{ 
    total_degreeS() {}
    total_degreeS(std::string scalar_property, const GraphInterface& g) {}
    template <class Graph, class Vertex> 
    size_t operator()(const Vertex& v, const Graph &g) const 
    {
	using namespace boost;
	typedef typename is_convertible<typename graph_traits<Graph>::directed_category, directed_tag>::type is_directed;
	return get_total_degree(v,g,is_directed());
    } 

    template <class Graph, class Vertex>
    size_t get_total_degree(const Vertex& v, const Graph &g, boost::true_type) const 
    {
	return in_degree(v,g)+out_degree(v,g);
    }

    template <class Graph, class Vertex>
    size_t get_total_degree(const Vertex& v, const Graph &g, boost::false_type) const 
    {
	return out_degree(v,g);
    }

    std::string name() {return "total_degree";}
};

struct in_degreeS 
{ 
    in_degreeS() {}
    in_degreeS(std::string scalar_property, const GraphInterface& g) {}
    template <class Graph, class Vertex> 
    size_t operator()(const Vertex& v, const Graph &g) const 
    {
	using namespace boost;
	typedef typename is_convertible<typename graph_traits<Graph>::directed_category, directed_tag>::type is_directed;
	return get_in_degree(v,g,is_directed());
    } 

    template <class Graph, class Vertex>
    size_t get_in_degree(const Vertex& v, const Graph &g, boost::true_type) const 
    {
	return in_degree(v,g);
    }

    template <class Graph, class Vertex>
    size_t get_in_degree(const Vertex& v, const Graph &g, boost::false_type) const 
    {
	return 0;
    }

    std::string name() {return "in_degree";}
};

struct out_degreeS 
{ 
    out_degreeS() {}
    out_degreeS(std::string scalar_property, const GraphInterface& g) {}
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
    scalarS(std::string scalar_property, const GraphInterface& g): 
	_scalar_property(scalar_property), _g(&g) {}
    typedef boost::mpl::vector<long double,double,float,long,unsigned long,int,unsigned int,short,unsigned short,char,unsigned char,bool,std::string> scalar_types;
    template <class Graph, class VertexOrEdge> 
    double operator()(const VertexOrEdge& v, const Graph &g) const 
    {
	try 
	{
	    return boost::get(_scalar_property, _g->_properties, v, boost::type<double>());
	}
	catch (boost::bad_any_cast)
	{
	    using namespace boost::mpl;
	    return get_value<next<begin<scalar_types>::type>::type>(v);
	}
    } 

    template <class ValueIter, class VertexOrEdge> 
    double get_value(const VertexOrEdge& v, ValueIter = ValueIter()) const 
    {
	using namespace boost;
	using namespace boost::mpl;
	try 
	{
	    return lexical_cast<double>(get(_scalar_property, _g->_properties, v, type<typename deref<ValueIter>::type>()));
	}
	catch (bad_any_cast)
	{
	    return get_value(v, typename boost::mpl::next<ValueIter>::type());
	}	
    }

    template <class VertexOrEdge> 
    double get_value(const VertexOrEdge& v, boost::mpl::end<scalar_types>::type) const 
    {
	throw boost::dynamic_get_failure(_scalar_property);
    }

    std::string name() {return _scalar_property;}
    std::string _scalar_property;
    GraphInterface const* _g;
};


typedef boost::mpl::map< boost::mpl::pair<in_degreeS, boost::mpl::int_<GraphInterface::IN_DEGREE> >,
			 boost::mpl::pair<out_degreeS, boost::mpl::int_<GraphInterface::OUT_DEGREE> >,
			 boost::mpl::pair<total_degreeS, boost::mpl::int_<GraphInterface::TOTAL_DEGREE> >,
			 boost::mpl::pair<scalarS, boost::mpl::int_<GraphInterface::SCALAR> > > degree_selector_index;

struct out_edgeS
{
    template <class Graph>
    struct iterator
    {
	typedef typename boost::graph_traits<Graph>::out_edge_iterator type;
    };

    template <class Graph>
    struct descriptor
    {
	typedef typename boost::graph_traits<Graph>::edge_descriptor type;
    };

    template <class Graph, class Vertex>
    std::pair<typename iterator<Graph>::type, typename iterator<Graph>::type> 
    operator()(const Vertex& v, const Graph &g) const
    {
	return out_edges(v,g);
    }

    template <class Graph>
    typename boost::graph_traits<Graph>::vertex_descriptor
    target(const typename descriptor<Graph>::type& e, const Graph &g) const
    {
	return boost::target(e,g);
    }
};

struct in_edgeS
{
    template <class Graph, class IsDirected>
    struct get_iterator
    {
	typedef typename boost::graph_traits<Graph>::in_edge_iterator type;
    };

    template <class Graph>
    struct get_iterator<Graph,boost::false_type>
    {
	typedef typename boost::graph_traits<Graph>::out_edge_iterator type;
    };

    template <class Graph>
    struct iterator
    {
	typedef typename boost::is_convertible<typename boost::graph_traits<Graph>::directed_category, boost::directed_tag>::type is_directed;
	typedef typename get_iterator<Graph, is_directed>::type type;
    };


    template <class Graph>
    struct descriptor
    {
	typedef typename boost::graph_traits<Graph>::edge_descriptor type;
    };

    template <class Graph, class Vertex>
    std::pair<typename iterator<Graph>::type, typename iterator<Graph>::type> 
    operator()(const Vertex& v, const Graph &g) const
    {
	typedef typename boost::is_convertible<typename boost::graph_traits<Graph>::directed_category, boost::directed_tag>::type is_directed;
	return get_in_edges(v, g, is_directed());
    }

    template <class Graph, class Vertex>
    std::pair<typename iterator<Graph>::type, typename iterator<Graph>::type> 
    get_in_edges(const Vertex& v, const Graph &g, boost::true_type) const
    {
	return in_edges(v,g);
    }
    
    template <class Graph, class Vertex>
    std::pair<typename iterator<Graph>::type, typename iterator<Graph>::type> 
    get_in_edges(const Vertex& v, const Graph &g, boost::false_type) const
    {
	typename iterator<Graph>::type end = out_edges(v,g).second;
	return make_pair(end,end);
    }

    template <class Graph>
    typename boost::graph_traits<Graph>::vertex_descriptor
    target(const typename descriptor<Graph>::type& e, const Graph &g) const
    {
	return source(e,g);
    }
};

struct any_edgeS
{
    template <class Graph, class IsDirected>
    struct get_iterator
    {
	typedef typename boost::graph_traits<boost::UndirectedAdaptor<Graph> >::out_edge_iterator type;
    };

    template <class Graph>
    struct get_iterator<Graph,boost::false_type>
    {
	typedef typename boost::graph_traits<Graph>::out_edge_iterator type;
    };

    template <class Graph>
    struct iterator
    {
	typedef typename boost::is_convertible<typename boost::graph_traits<Graph>::directed_category, boost::directed_tag>::type is_directed;
	typedef typename get_iterator<Graph, is_directed>::type type;
    };

    template <class Graph, class IsDirected>
    struct get_descriptor
    {
	typedef typename boost::graph_traits<boost::UndirectedAdaptor<Graph> >::edge_descriptor type;
    };

    template <class Graph>
    struct get_descriptor<Graph,boost::false_type>
    {
	typedef typename boost::graph_traits<Graph>::edge_descriptor type;
    };

    template <class Graph>
    struct descriptor
    {
	typedef typename boost::is_convertible<typename boost::graph_traits<Graph>::directed_category, boost::directed_tag>::type is_directed;
	typedef typename get_descriptor<Graph, is_directed>::type type;
    };

    template <class Graph, class Vertex>
    std::pair<typename iterator<Graph>::type, typename iterator<Graph>::type> 
    operator()(const Vertex& v, const Graph &g) const
    {
	typedef typename boost::is_convertible<typename boost::graph_traits<Graph>::directed_category, boost::directed_tag>::type is_directed;
	return get_all_edges(v,g, is_directed());
    }

    template <class Graph, class Vertex>
    std::pair<typename iterator<Graph>::type, typename iterator<Graph>::type> 
    get_all_edges(const Vertex& v, const Graph &g, boost::true_type) const
    {
	boost::UndirectedAdaptor<Graph> ug(g);
	return out_edges(v,ug);
    }
    
    template <class Graph, class Vertex>
    std::pair<typename iterator<Graph>::type, typename iterator<Graph>::type> 
    get_all_edges(const Vertex& v, const Graph &g, boost::false_type) const
    {
	return out_edges(v,g);
    }

    template <class Graph>
    typename boost::graph_traits<Graph>::vertex_descriptor
    target(const typename descriptor<Graph>::type& e, const Graph &g) const
    {
	typedef typename boost::is_convertible<typename boost::graph_traits<Graph>::directed_category, boost::directed_tag>::type is_directed;
	return get_target(e,g,is_directed());
    }

    template <class Graph>
    typename boost::graph_traits<Graph>::vertex_descriptor
    get_target(const typename descriptor<Graph>::type& e, const Graph &g, boost::true_type) const
    {
	boost::UndirectedAdaptor<Graph> ug(g);
	return boost::target(e,ug);
    }

    template <class Graph>
    typename boost::graph_traits<Graph>::vertex_descriptor
    get_target(const typename descriptor<Graph>::type& e, const Graph &g, boost::false_type) const
    {
	return boost::target(e,g);
    }
};


typedef boost::mpl::map< boost::mpl::pair<out_edgeS, boost::mpl::int_<GraphInterface::OUT_NEIGHBOURS> >,
			 boost::mpl::pair<in_edgeS, boost::mpl::int_<GraphInterface::IN_NEIGHBOURS> >,
			 boost::mpl::pair<any_edgeS, boost::mpl::int_<GraphInterface::ALL_NEIGHBOURS> > >::type  edges_selector_index;


} //namespace graph_tool

#endif
