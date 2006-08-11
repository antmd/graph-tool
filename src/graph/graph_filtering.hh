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

namespace graph_tool
{
using namespace boost;

//==============================================================================
// HardNumVertices()
//==============================================================================
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

//==============================================================================
// SoftNumVertices()
//==============================================================================
struct SoftNumVertices
{
    template <class Graph>
    size_t operator()(const Graph &g) const { return num_vertices(g); }
};

//==============================================================================
// HardNumEdges()
//==============================================================================
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

//==============================================================================
// SoftNumEdges()
//==============================================================================
struct SoftNumEdges
{
    template <class Graph> 
    size_t operator()(const Graph &g) const { return num_edges(g); }
};

//==============================================================================
// RangeFilter
//==============================================================================
template <class FilteredPropertyMap>
class RangeFilter
{
public:
    RangeFilter(){}
    typedef typename property_traits<FilteredPropertyMap>::value_type value_type;
    typedef typename property_traits<FilteredPropertyMap>::key_type key_type;
    RangeFilter(FilteredPropertyMap filtered_property, std::pair<value_type, value_type> range)
        : _filtered_property(filtered_property), _range(range) {}
    template <class VertexOrEdge>
    bool operator() (VertexOrEdge e) const 
    {
        // ignore if outside allowed range
        if ( _filtered_property[e] < _range.first || _filtered_property[e] > _range.second)
            return false;
        return true;
    }
private:
    FilteredPropertyMap _filtered_property;
    std::pair<value_type, value_type> _range;
};

//==============================================================================
// PythonFilter
//==============================================================================
template <class Graph, class IndexMap, class HasBase = mpl::bool_<false> >
class PythonFilter
{
public:
    PythonFilter(){}
    PythonFilter(const Graph& g, IndexMap index_map, const dynamic_properties& dp, python::object filter)
        : _g(&g), _index_map(index_map), _dp(&dp), _filter(filter) {}

    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

    struct get_value
    {
	get_value(const std::string& name, dynamic_property_map& map, any key, python::dict& d)
	    : _name(name), _map(map), _key(key), _d(d) {}
	template <class Type>
	void operator()(Type)
	{
	    try 
	    {
		_d[_name] = any_cast<Type>(const_cast<dynamic_property_map&>(_map).get(_key));
	    }
	    catch (bad_any_cast){}
	}
	const std::string& _name;
	const dynamic_property_map& _map;
	any _key;
	python::dict& _d;
    };

    //FIXME: It would be better to do a specialization for vertex_descriptor, and a general
    //       dummy put_degree. But GCC doesn't seem to like nested full template specializations
    template <class Vertex>
    struct put_degree
    {
	put_degree(const Graph& g, Vertex v, python::dict& d, std::string prefix = "")
	    : _g(g), _v(v), _d(d), _prefix(prefix) {}
	template <class Degree>
	void operator()(Degree degree)
	{
	    _d[_prefix+degree.name()] = degree(_v, _g);
	}
	const Graph& _g;
	Vertex _v;
	python::dict& _d;
	std::string _prefix;
    };
    
    template <class Vertex>
    struct put_base_degree: public put_degree<Vertex>
    {
	put_base_degree(const Graph& g, Vertex v, python::dict& d, std::string prefix = ""): put_degree<Vertex>(g,v,d,prefix) {}

	template <class Degree>
	void operator()(Degree degree)
	{
	    this->_d[this->_prefix+degree.name()] = degree(this->_v, this->_g.m_g);
	}
    };

    typedef mpl::vector<in_degreeS,out_degreeS,total_degreeS> degrees;

    inline void put_edge_info(edge_descriptor e, python::dict& d) const
    {
	if (source(e,*_g) == target(e,*_g))
	    d["is_loop"] = true;
	else
	    d["is_loop"] = false;
	for(typeof(_dp->begin()) iter = _dp->begin(); iter != _dp->end(); ++iter)
	{
	    if (iter->second->key() == typeid(vertex_descriptor))
	    {
		typedef mpl::vector<bool,int,long,size_t,float,double,std::string> value_types;
		mpl::for_each<value_types>(get_value("source_"+iter->first, *iter->second, source(e, *_g), d));
		mpl::for_each<value_types>(get_value("target_"+iter->first, *iter->second, target(e, *_g), d));
	    }
	}
	mpl::for_each<degrees>(put_degree<vertex_descriptor>(*_g, source(e,*_g), d, "source_"));
	mpl::for_each<degrees>(put_degree<vertex_descriptor>(*_g, target(e,*_g), d, "target_"));	
    }

    inline void put_edge_info(vertex_descriptor v, python::dict& d) const
    {
    }
    
    template <class VertexOrEdge>
    inline bool operator() (VertexOrEdge e) const 
    {
	BOOST_MPL_ASSERT(( mpl::or_<is_same<VertexOrEdge,vertex_descriptor>,
                                    is_same<VertexOrEdge,edge_descriptor> > ));
	python::dict properties;
	for(typeof(_dp->begin()) iter = _dp->begin(); iter != _dp->end(); ++iter)
	{
	    if (iter->second->key() == typeid(VertexOrEdge))
	    {
		typedef mpl::vector<bool,int,long,size_t,float,double,std::string> value_types;
		mpl::for_each<value_types>(get_value(iter->first, *iter->second, e, properties));
	    }
	}

	typedef typename mpl::if_<is_same<VertexOrEdge,vertex_descriptor>, degrees, mpl::vector<> >::type vertex_degrees;
	mpl::for_each<vertex_degrees>(put_degree<VertexOrEdge>(*_g, e, properties));

	typedef typename mpl::if_<HasBase, vertex_degrees, mpl::vector<> >::type base_degrees;
	mpl::for_each<base_degrees>(put_base_degree<VertexOrEdge>(*_g, e, properties, "orig_"));

	put_edge_info(e, properties);
	
	return python::extract<bool>(_filter(properties));
    }

private:
    Graph const*  _g;
    IndexMap _index_map;
    dynamic_properties const*  _dp;
    python::object _filter;
};


typedef mpl::vector<mpl::bool_<true>, mpl::bool_<false> > reverse_check;
typedef mpl::vector<mpl::bool_<false> > never_reversed;
typedef mpl::vector<mpl::bool_<true> > always_reversed;
typedef mpl::vector<mpl::bool_<true>, mpl::bool_<false> > directed_check;
typedef mpl::vector<mpl::bool_<true> > always_directed;
typedef mpl::vector<mpl::bool_<false> > always_undirected;

template <class Graph, class Action>
struct check_reverse
{
    check_reverse(const Graph &g, Action a, bool reverse, bool& found): _g(g), _a(a), _reverse(reverse), _found(found) {}

    template <class Reverse>
    void operator()(Reverse) const
    {
	if (_reverse)
	{
	    reverse_graph<Graph> rg(_g);
	    _a(rg);
	    _found = true;
	}
    }

    void operator()(mpl::bool_<false>) const
    { 
	if (!_reverse)
	{
	    _a(_g);
	    _found = true;
	}
    }

    const Graph &_g;
    Action _a;
    bool _reverse;
    bool& _found;
};

template <class Graph, class Action, class ReverseCheck>
struct check_directed
{
    check_directed(const Graph &g, Action a, bool reverse, bool directed, bool& found)
	: _g(g), _a(a), _reverse(reverse), _directed(directed), _found(found) {}

    template <class Directed>
    void operator()(Directed)
    {
	if (_directed)
	    mpl::for_each<ReverseCheck>(check_reverse<Graph, Action>(_g, _a, _reverse, _found));
    }

    void operator()(mpl::bool_<false>)
    { 
	if (!_directed)
	{
	    UndirectedAdaptor<Graph> ug(_g);
	    _a(ug);
	    _found = true;
	}
    }

    const Graph &_g;
    Action _a;
    bool _reverse;
    bool _directed;
    bool& _found;
};

template <class Graph, class Action, class ReverseCheck, class DirectedCheck> 
void check_python_filter(const Graph& g, const GraphInterface &gi, Action a, bool& found, ReverseCheck, DirectedCheck)
{
    typedef PythonFilter<Graph,GraphInterface::vertex_index_map_t> vertex_filter_t;
    typedef PythonFilter<Graph,GraphInterface::edge_index_map_t> edge_filter_t;

    if (gi._edge_python_filter != python::object())
    {
	typedef filtered_graph<Graph, edge_filter_t, keep_all> efg_t;
	efg_t efg(g,edge_filter_t(g, gi._edge_index, gi._properties, gi._edge_python_filter), keep_all());

	if (gi._vertex_python_filter != python::object())
	{
	    typedef PythonFilter<efg_t,GraphInterface::vertex_index_map_t, mpl::bool_<true> > vertex_filter_t;
	    typedef filtered_graph<efg_t,keep_all,vertex_filter_t> vefg_t;
	    vefg_t vefg(efg,keep_all(),vertex_filter_t(efg, gi._vertex_index, gi._properties, gi._vertex_python_filter));
	    mpl::for_each<DirectedCheck>(check_directed<vefg_t,Action,ReverseCheck>(vefg, a, gi._reversed, gi._directed, found));
	}
	else
	{
	    mpl::for_each<DirectedCheck>(check_directed<efg_t,Action,ReverseCheck>(efg, a, gi._reversed, gi._directed, found));
	}
    }
    else if (gi._vertex_python_filter != python::object())
    {
	typedef filtered_graph<Graph,keep_all,vertex_filter_t> vfg_t;
	vfg_t vfg(g,keep_all(),vertex_filter_t(g, gi._vertex_index, gi._properties, gi._vertex_python_filter));
	mpl::for_each<DirectedCheck>(check_directed<vfg_t,Action,ReverseCheck>(vfg, a, gi._reversed, gi._directed, found));
    } 
    else
    {
	mpl::for_each<DirectedCheck>(check_directed<Graph,Action,ReverseCheck>(g, a, gi._reversed, gi._directed, found));
    }
}

template <class Action, class ReverseCheck, class DirectedCheck> 
void check_filter(const GraphInterface &g, Action a, ReverseCheck, DirectedCheck)
{
    typedef RangeFilter<GraphInterface::vertex_filter_map_t> vertex_filter_t;
    typedef RangeFilter<GraphInterface::edge_filter_map_t> edge_filter_t;
    
    bool found = false;

    if (g._edge_python_filter == python::object() && g._vertex_python_filter == python::object())
    {
	if (g._vertex_filter_property != "" && g._edge_filter_property != "")
	{	
	    typedef filtered_graph<GraphInterface::multigraph_t, edge_filter_t, vertex_filter_t> fg_t;
	    fg_t fg(g._mg, edge_filter_t(g._edge_filter_map, g._edge_range), vertex_filter_t(g._vertex_filter_map, g._vertex_range));
	    mpl::for_each<DirectedCheck>(check_directed<fg_t,Action,ReverseCheck>(fg, a, g._reversed, g._directed, found));
	}
	else if (g._vertex_filter_property != "")
	{
	    typedef filtered_graph<GraphInterface::multigraph_t, keep_all, vertex_filter_t> fg_t;
	    fg_t fg(g._mg, keep_all(), vertex_filter_t(g._vertex_filter_map, g._vertex_range));
	    mpl::for_each<DirectedCheck>(check_directed<fg_t,Action,ReverseCheck>(fg, a, g._reversed, g._directed, found));
	} 
	else if (g._edge_filter_property != "")
	{
	    typedef filtered_graph<GraphInterface::multigraph_t, edge_filter_t, keep_all> fg_t;
	    fg_t fg(g._mg, edge_filter_t(g._edge_filter_map, g._edge_range), keep_all());
	    mpl::for_each<DirectedCheck>(check_directed<fg_t,Action,ReverseCheck>(fg, a, g._reversed, g._directed, found));
	}
	else
	{
	    mpl::for_each<DirectedCheck>(check_directed<GraphInterface::multigraph_t,Action,ReverseCheck>(g._mg, a, g._reversed, g._directed, found));
	}
    }
    else
    {
	check_python_filter(g._mg, g, a, found, ReverseCheck(), DirectedCheck());
    }

    if (!found)
	throw GraphException("graph filtering error: filter not found");
    
}

template <class Descriptor, class IndexMap>
class DescriptorHash: public std::unary_function<Descriptor, std::size_t> 
{
public:
    DescriptorHash() {}
    DescriptorHash(IndexMap index_map): _index_map(index_map) {}
    std::size_t operator()(Descriptor const& d) const { return boost::hash_value(_index_map[d]); }
private:
    IndexMap _index_map;
};

} //namespace graph_tool

#endif
