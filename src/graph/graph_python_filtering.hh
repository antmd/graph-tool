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

#ifndef PYTHON_FILTERING_HH
#define PYTHON_FILTERING_HH

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
#include <boost/python/class.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/long.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/make_function.hpp>
#include <boost/python/self.hpp>
#include <boost/python/operators.hpp>
#include <boost/python/iterator.hpp>
#include <boost/functional/hash.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/functional.hpp>

namespace graph_tool
{
using namespace boost;

template <class VertexOrEdge>
class get_python_property_value
{
public:
    get_python_property_value(const dynamic_property_map& dmap, const VertexOrEdge& e)
        : _dmap(dmap), _e(e) {}
        
    python::object operator()()
    {
        typedef mpl::vector<bool,int,long,size_t,float,double,std::string> value_types;
        mpl::for_each<value_types>(try_conversion(*this));
        return _retval;
    }

private:
    struct try_conversion
    {
        try_conversion(get_python_property_value& parent): _parent(parent) {}
        
        template <class Type>
        void operator()(Type)
        {
            any any_val = const_cast<dynamic_property_map&>(_parent._dmap).get(_parent._e);
            Type* value = any_cast<Type>(&any_val);
            if (value != 0)
                _parent._retval = python::object(*value);
        }        
        get_python_property_value& _parent;
    };
    
    const dynamic_property_map& _dmap;
    const VertexOrEdge& _e;
    python::object _retval;
};

template <class VertexOrEdge>
class put_python_property_value
{
public:
    put_python_property_value(dynamic_property_map& dmap, const VertexOrEdge& e, python::object val)
        : _dmap(dmap), _e(e), _val(val) {}
        
    void operator()()
    {
        typedef mpl::vector<bool,int,long,size_t,float,double,std::string> value_types;
        mpl::for_each<value_types>(try_conversion(*this));
    }

private:
    struct try_conversion
    {
        try_conversion(put_python_property_value& parent): _parent(parent) {}
        
        template <class Type>
        void operator()(Type)
        {
            if (typeid(Type) == _parent._dmap.value())
            {
                Type val = python::extract<Type>(_parent._val);
                _parent._dmap.put(_parent._e, val);
            }
        }
        
        put_python_property_value& _parent;
    };
    
    dynamic_property_map& _dmap;
    const VertexOrEdge& _e;
    python::object _val;
};


template <class Graph>
class PythonEdge;

//==============================================================================
// PythonVertex
//==============================================================================

template <class Graph, class IsDirected> 
struct in_edge_iterator
{
    BOOST_MPL_ASSERT((is_same<IsDirected, boost::true_type>));
    BOOST_MPL_ASSERT((is_convertible<typename graph_traits<Graph>::directed_category, boost::directed_tag>));
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef typename graph_traits<Graph>::in_edge_iterator type;
    static std::pair<type,type> in_edges(vertex_descriptor v, const Graph& g) { return boost::in_edges(v,g); }
};

template <class Graph>
struct in_edge_iterator<Graph, boost::false_type>
{
    BOOST_MPL_ASSERT((is_convertible<typename graph_traits<Graph>::directed_category, boost::undirected_tag>));
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef typename graph_traits<Graph>::out_edge_iterator type;
    static std::pair<type,type> in_edges(vertex_descriptor v, const Graph& g) { return make_pair(type(), type()); }
};


template <class Graph>
class PythonVertex
{
public:
    PythonVertex(python::object base): _base(base)
    {
        static bool first_time = true;
        if (first_time)
        {
            _python_class.def(python::init<python::object>());
            _python_class.def("in_degree", &PythonVertex::GetInDegree);
            _python_class.def("out_degree", &PythonVertex::GetOutDegree);
            _python_class.def("get_property", &PythonVertex::GetProperty);
            _python_class.def("put_property", &PythonVertex::PutProperty);
            _python_class.def("out_edges", python::range(&PythonVertex::OutEdgesBegin, &PythonVertex::OutEdgesEnd));
            _python_class.def("in_edges", python::range(&PythonVertex::InEdgesBegin, &PythonVertex::InEdgesEnd));
            _python_class.def(python::self == python::self);
            _python_class.def(python::self != python::self);
            _python_class.def("__hash__", &PythonVertex::GetHash);
            first_time = false;
        }
    }
    PythonVertex(): _base(python::object()) { *this = PythonVertex(_base); }

    PythonVertex(const PythonVertex& v)
    {
        if (v._base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            if (get<3>(base))
            {
                vertex_descriptor* new_v = new vertex_descriptor;
                *new_v = get<1>(base);
                base_t new_base(get<0>(base), *new_v, get<2>(base), true);
                _base = python::object(new_base);
            }
        }
    }

    ~PythonVertex()
    {
        if (_base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            if (get<3>(base))
                delete &get<1>(base);
        }
    }

    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
    typedef tuple<const Graph&, const vertex_descriptor&, const dynamic_properties&, bool> base_t;

    python::object GetInDegree()
    {
        if (_base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            return python::object(in_degreeS()(get<1>(base), get<0>(base)));
        }
        else
        {
            return python::object();
        }
    }

    python::object GetOutDegree()
    {
        if (_base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            return python::object(out_degree(get<1>(base), get<0>(base)));
        }
        else
        {
            return python::object();
        }
    }

    python::object GetProperty(std::string prop)
    {
        if (_base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            const dynamic_properties& dp = get<2>(base);
            const vertex_descriptor& v = get<1>(base);
            for(typeof(dp.begin()) iter = dp.lower_bound(prop); iter != dp.end(); ++iter)
            {
                if (iter->second->key() == typeid(vertex_descriptor))
                    return get_python_property_value<vertex_descriptor>(*iter->second, v)();
            }
        }
        return python::object();
    }

    void PutProperty(std::string prop, python::object val)
    {
        if (_base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            const dynamic_properties& dp = get<2>(base);
            const vertex_descriptor& v = get<1>(base);
            for(typeof(dp.begin()) iter = dp.lower_bound(prop); iter != dp.end(); ++iter)
            {
                if (iter->second->key() == typeid(vertex_descriptor))
                    put_python_property_value<vertex_descriptor>(*iter->second, v, val)();
            }
        }
    }


    struct descriptor_wrap
    {
        descriptor_wrap(const Graph &g, const dynamic_properties& dp): _g(g), _dp(dp) {}
    
        python::object operator()(edge_descriptor e)
        {
            edge_descriptor *new_e = new edge_descriptor;
            *new_e = e;
            typename PythonEdge<Graph>::base_t edge_base(_g, *new_e, _dp, true);
            python::object edge = (PythonEdge<Graph>().GetPythonClass())(python::object(edge_base));
            return edge;
        }

        const Graph& _g;
        const dynamic_properties& _dp;
    };

    typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iterator;
    typedef function<python::object (edge_descriptor)> edge_wrap_function;
    typedef transform_iterator<edge_wrap_function, out_edge_iterator> wrapped_out_edge_iterator;

    wrapped_out_edge_iterator OutEdgesBegin()
    {
        if (_base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            const Graph& g = get<0>(base);
            vertex_descriptor v = get<1>(base);
            const dynamic_properties& dp(get<2>(base));
            return wrapped_out_edge_iterator(out_edges(v,g).first, edge_wrap_function(descriptor_wrap(g,dp)));
        }
        return wrapped_out_edge_iterator();
    }

    wrapped_out_edge_iterator OutEdgesEnd()
    {
        if (_base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            const Graph& g = get<0>(base);
            vertex_descriptor v = get<1>(base);
            const dynamic_properties& dp(get<2>(base));
            return wrapped_out_edge_iterator(out_edges(v,g).second, edge_wrap_function(descriptor_wrap(g,dp)));
        }
        return wrapped_out_edge_iterator();
    }

    typedef typename is_convertible<typename graph_traits<Graph>::directed_category,boost::directed_tag>::type is_directed;
    typedef typename in_edge_iterator<Graph,is_directed>::type in_edge_iterator_type;
    typedef transform_iterator<edge_wrap_function, in_edge_iterator_type> wrapped_in_edge_iterator;

    wrapped_in_edge_iterator InEdgesBegin()
    {
        if (_base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            const Graph& g = get<0>(base);
            vertex_descriptor v = get<1>(base);
            const dynamic_properties& dp(get<2>(base));
            return wrapped_in_edge_iterator(in_edge_iterator<Graph,is_directed>::in_edges(v,g).first, edge_wrap_function(descriptor_wrap(g,dp)));
        }
        return wrapped_in_edge_iterator();
    }

    wrapped_in_edge_iterator InEdgesEnd()
    {
        if (_base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            const Graph& g = get<0>(base);
            vertex_descriptor v = get<1>(base);
            const dynamic_properties& dp(get<2>(base));
            return wrapped_in_edge_iterator(in_edge_iterator<Graph,is_directed>::in_edges(v,g).second, edge_wrap_function(descriptor_wrap(g,dp)));
        }
        return wrapped_in_edge_iterator();
    }

    python::object GetHash()
    {
        if (_base != python::object())
        {            
            base_t base = python::extract<base_t>(_base);
            return python::object(hash<vertex_descriptor>()(get<1>(base)));
        }
        else
        {
            return python::object(0);
        }
    }

    bool operator==(const PythonVertex& other)
    {
        if (_base != python::object() && other._base != python::object())
        {            
            const base_t o_base = python::extract<base_t>(other._base);
            base_t base = python::extract<base_t>(_base);
            return python::object(get<1>(base) == get<1>(o_base));
        }
        else
        {            
            return false;
        }
    }

    bool operator!=(const PythonVertex& other)
    {
        return !(*this == other);
    }

    python::class_<PythonVertex> GetPythonClass() { return _python_class; }
        
private:
    python::object _base;
    static python::class_<PythonVertex> _python_class;
};

template <class Graph> python::class_<PythonVertex<Graph> > PythonVertex<Graph>::_python_class =  python::class_<PythonVertex<Graph> >("Vertex", python::no_init);

//==============================================================================
// PythonEdge
//==============================================================================

template <class Graph>
class PythonEdge
{
public:
    PythonEdge(python::object base): _base(base)
    {
        static bool first_time = true;
        if (first_time)
        {
            _python_class.def(python::init<python::object>());
            _python_class.def("source", &PythonEdge::GetSource);
            _python_class.def("target", &PythonEdge::GetTarget);
            _python_class.def("get_property", &PythonEdge::GetProperty);
            _python_class.def("put_property", &PythonEdge::PutProperty);
            _python_class.def("n_parallel", &PythonEdge::GetNParallel);
            _python_class.def(python::self == python::self);
            _python_class.def(python::self != python::self);
            _python_class.def("__hash__", &PythonEdge::GetHash);
            first_time = false;
        }
    }
    PythonEdge(): _base(python::object()) { *this = PythonEdge(_base); }
    PythonEdge(const PythonEdge& e)
    {
        if (e._base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            if (get<3>(base))
            {
                edge_descriptor* new_e = new edge_descriptor;
                base_t new_base(get<0>(base), *new_e, get<2>(base), true);
                _base = python::object(new_base);
            }
        }
    }

    ~PythonEdge()
    {
        if (_base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            if (get<3>(base))
                delete &get<1>(base);
        }
    }

    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
    typedef tuple<const Graph&, const edge_descriptor&, const dynamic_properties&, bool> base_t;

    python::object GetSource()
    {
        if (_base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            vertex_descriptor *v = new vertex_descriptor;
            *v = source(get<1>(base), get<0>(base));
            typename PythonVertex<Graph>::base_t vertex_base(get<0>(base), *v, get<2>(base), true);
            python::object source_vertex = (PythonVertex<Graph>().GetPythonClass())(python::object(vertex_base));
            return source_vertex;
        }
        else
        {
            return python::object();
        }
    }

    python::object GetTarget()
    {
        if (_base != python::object())
        {
            base_t base = python::extract<base_t>(_base);            
            vertex_descriptor *v = new vertex_descriptor;
            *v = target(get<1>(base), get<0>(base));
            typename PythonVertex<Graph>::base_t vertex_base(get<0>(base), *v, get<2>(base), true);
            python::object target_vertex = (PythonVertex<Graph>().GetPythonClass())(python::object(vertex_base));
            return target_vertex;
        }
        else
        {
            return python::object();
        }
    }

    python::object GetProperty(std::string prop)
    {
        if (_base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            const dynamic_properties& dp = get<2>(base);
            const edge_descriptor& e = get<1>(base);
            for(typeof(dp.begin()) iter = dp.lower_bound(prop); iter != dp.end(); ++iter)
            {
                if (iter->second->key() != typeid(vertex_descriptor) && iter->second->key() != typeid(graph_property_tag))
                    return do_get_property(*iter->second, e, typename is_convertible<typename graph_traits<Graph>::directed_category, undirected_tag>::type());
            }
            return python::object();
        }
        else
        {
            return python::object();
        }
    }

    python::object do_get_property(const dynamic_property_map& dmap, const edge_descriptor& e, false_type)
    {
        return get_python_property_value<edge_descriptor>(dmap, e)();
    }

    python::object do_get_property(const dynamic_property_map& dmap, const edge_descriptor& e, true_type)
    {
        typedef typename edge_descriptor::original_edge_t original_edge_t;
        return get_python_property_value<original_edge_t>(dmap, original_edge_t(e))();
    }


    void PutProperty(const std::string& prop, const python::object& val)
    {
        if (_base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            const dynamic_properties& dp = get<2>(base);
            const edge_descriptor& e = get<1>(base);
            for(typeof(dp.begin()) iter = dp.lower_bound(prop); iter != dp.end(); ++iter)
            {
                if (iter->second->key() != typeid(vertex_descriptor) && iter->second->key() != typeid(graph_property_tag))
                    do_put_property(*iter->second, e, val, typename is_convertible<typename graph_traits<Graph>::directed_category, undirected_tag>::type());
            }
        }
    }

    void do_put_property(dynamic_property_map& dmap, const edge_descriptor& e, python::object val, false_type)
    {
        put_python_property_value<edge_descriptor>(dmap, e, val)();
    }

    void do_put_property(dynamic_property_map& dmap, const edge_descriptor& e, python::object val, true_type)
    {
        typedef typename edge_descriptor::original_edge_t original_edge_t;
        const original_edge_t& oe = e;
        put_python_property_value<original_edge_t>(dmap, oe, val)();
    }

    python::object GetNParallel()
    {
        if (_base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            const Graph& g = get<0>(base);
            const edge_descriptor& e = get<1>(base);
            size_t n = 0;
            vertex_descriptor s,t;
            s = source(e, g);
            t = target(e, g);
            typename graph_traits<Graph>::adjacency_iterator a, a_end;
            for(tie(a, a_end) = adjacent_vertices(s, g); a != a_end; ++a)
                if (*a == t)
                    n++;
            return python::object(n-1);
        }
        else
        {
            return python::object();
        }
    };

    python::object GetHash()
    {
        if (_base != python::object())
        {            
            base_t base = python::extract<base_t>(_base);
            const Graph& g = get<0>(base);
            const edge_descriptor& e = get<1>(base);
            vertex_descriptor s,t;
            s = source(e, g);
            t = target(e, g);
            return python::object(hash<std::pair<vertex_descriptor,vertex_descriptor> >()(make_pair(s,t)));
        }
        else
        {
            return python::object(0);
        }
    }

    bool operator==(const PythonEdge& other)
    {
        if (_base != python::object() && other._base != python::object())
        {            
            const base_t o_base = python::extract<base_t>(other._base);
            base_t base = python::extract<base_t>(_base);
            return python::object(get<1>(base) == get<1>(o_base));
        }
        else
        {
            return false;
        }
    }

    bool operator!=(const PythonEdge& other)
    {
        return !(*this == other);
    }

    python::class_<PythonEdge> GetPythonClass() { return _python_class; }
        
private:

    python::object _base;
    static python::class_<PythonEdge> _python_class;
};

template <class Graph> python::class_<PythonEdge<Graph> > PythonEdge<Graph>::_python_class =  python::class_<PythonEdge<Graph> >("Edge", python::no_init);

//==============================================================================
// PythonGraph
//==============================================================================

template <class Graph>
class PythonGraph
{
public:
    PythonGraph(python::object base): _base(base)
    {
        static bool first_time = true;
        if (first_time)
        {
            _python_class.def(python::init<python::object>());
            _python_class.def("num_vertices", &PythonGraph::GetNumVertices);
            _python_class.def("num_edges", &PythonGraph::GetNumEdges);
            _python_class.def("get_vertex", &PythonGraph::GetVertex);
            _python_class.def("get_edge", &PythonGraph::GetEdge);
            _python_class.def("get_property", &PythonGraph::GetProperty);
            _python_class.def("vertices", python::range(&PythonGraph::VertexBegin, &PythonGraph::VertexEnd));
            _python_class.def("edges", python::range(&PythonGraph::EdgeBegin, &PythonGraph::EdgeEnd));
            first_time = false;
        }
    }
    PythonGraph(): _base(python::object()) { *this = PythonGraph(_base); }

    PythonGraph(const PythonGraph& g)
    {
        _base = g._base;
    }

    ~PythonGraph()
    {
    }

    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
    typedef tuple<const Graph&, const edge_descriptor&, const dynamic_properties&, bool> base_t;

    python::object GetNumVertices()
    {
        if (_base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            const Graph& g = get<0>(base);

            return python::object(HardNumVertices()(g));
        }
        else
        {
            return python::object();
        }
    }

    python::object GetNumEdges()
    {
        if (_base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            const Graph& g = get<0>(base);

            return python::object(HardNumEdges()(g));
        }
        else
        {
            return python::object();
        }
    }

    python::object GetVertex(size_t v_index)
    {
        if (_base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            const Graph& g = get<0>(base);
            
            size_t count = 0;
            typename graph_traits<Graph>::vertex_iterator v, v_end;
            for (tie(v,v_end) = vertices(g); v != v_end; ++v)
            {
                if (count == v_index)
                {
                    vertex_descriptor *new_v = new vertex_descriptor;
                    *new_v = *v;
                    typename PythonVertex<Graph>::base_t vertex_base(g, *new_v, get<2>(base), true);
                    python::object vertex = (PythonVertex<Graph>().GetPythonClass())(python::object(vertex_base));
                    return vertex;
                }
                count++;
            }
        }
        return python::object();
    }

    python::object GetEdge(size_t e_index)
    {
        if (_base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            const Graph& g = get<0>(base);
            
            size_t count = 0;
            typename graph_traits<Graph>::edge_iterator e, e_end;
            for (tie(e,e_end) = edges(g); e != e_end; ++e)
            {
                if (count == e_index)
                {
                    edge_descriptor *new_e = new edge_descriptor;
                    *new_e = *e;
                    typename PythonEdge<Graph>::base_t edge_base(g, *new_e, get<2>(base), true);
                    python::object edge = (PythonEdge<Graph>().GetPythonClass())(python::object(edge_base));
                    return edge;
                }
                count++;
            }
        }
        return python::object();
    }

    python::object GetProperty(std::string prop)
    {
        if (_base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            const dynamic_properties& dp = get<2>(base);
            for(typeof(dp.begin()) iter = dp.lower_bound(prop); iter != dp.end(); ++iter)
            {
                if (iter->second->key() == typeid(graph_property_tag))
                    return get_python_property_value<graph_property_tag>(*iter->second, graph_property_tag())();
            }
        }
        return python::object();
    }

    void PutProperty(std::string prop, python::object val)
    {
        if (_base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            const dynamic_properties& dp = get<2>(base);
            for(typeof(dp.begin()) iter = dp.lower_bound(prop); iter != dp.end(); ++iter)
            {
                if (iter->second->key() == typeid(graph_property_tag))
                    put_python_property_value<vertex_descriptor>(*iter->second, graph_property_tag(), val)();
            }
        }
    }

    struct descriptor_wrap
    {
        descriptor_wrap(const Graph &g, const dynamic_properties& dp): _g(g), _dp(dp) {}
        python::object operator()(vertex_descriptor v)
        {
            vertex_descriptor *new_v = new vertex_descriptor;
            *new_v = v;
            typename PythonVertex<Graph>::base_t vertex_base(_g, *new_v, _dp, true);
            python::object vertex = (PythonVertex<Graph>().GetPythonClass())(python::object(vertex_base));
            return vertex;
        }

        python::object operator()(edge_descriptor e)
        {
            edge_descriptor *new_e = new edge_descriptor;
            *new_e = e;
            typename PythonEdge<Graph>::base_t edge_base(_g, *new_e, _dp, true);
            python::object edge = (PythonEdge<Graph>().GetPythonClass())(python::object(edge_base));
            return edge;
        }

        const Graph& _g;
        const dynamic_properties& _dp;
    };

    typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
    typedef function<python::object (vertex_descriptor)> vertex_wrap_function;
    typedef transform_iterator<vertex_wrap_function, vertex_iterator> wrapped_vertex_iterator;

    wrapped_vertex_iterator VertexBegin()
    {
        if (_base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            const Graph& g = get<0>(base);
            const dynamic_properties& dp(get<2>(base));
            return wrapped_vertex_iterator(vertices(g).first, vertex_wrap_function(descriptor_wrap(g,dp)));
        }
        return wrapped_vertex_iterator();
    }

    wrapped_vertex_iterator VertexEnd()
    {
        if (_base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            const Graph& g = get<0>(base);
            const dynamic_properties& dp(get<2>(base));
            return wrapped_vertex_iterator(vertices(g).second, vertex_wrap_function(descriptor_wrap(g,dp)));
        }
        return wrapped_vertex_iterator();
    }

    typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
    typedef function<python::object (edge_descriptor)> edge_wrap_function;
    typedef transform_iterator<edge_wrap_function, edge_iterator> wrapped_edge_iterator;

    wrapped_edge_iterator EdgeBegin()
    {
        if (_base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            const Graph& g = get<0>(base);
            const dynamic_properties& dp(get<2>(base));
            return wrapped_edge_iterator(edges(g).first, edge_wrap_function(descriptor_wrap(g,dp)));
        }
        return wrapped_edge_iterator();
    }

    wrapped_edge_iterator EdgeEnd()
    {
        if (_base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            const Graph& g = get<0>(base);
            const dynamic_properties& dp(get<2>(base));
            return wrapped_edge_iterator(edges(g).second, edge_wrap_function(descriptor_wrap(g,dp)));
        }
        return wrapped_edge_iterator();
    }
    
    python::class_<PythonGraph> GetPythonClass() { return _python_class; }
        
private:

    python::object _base;
    static python::class_<PythonGraph> _python_class;
};

template <class Graph> python::class_<PythonGraph<Graph> > PythonGraph<Graph>::_python_class =  python::class_<PythonGraph<Graph> >("Graph", python::no_init);


//==============================================================================
// populate_python_funcs
//==============================================================================

template <class Graph>
struct init_counter
{
    static bool registered;
};
template <class Graph> bool init_counter<Graph>::registered = false;

template <class Descriptor, class HasBase = mpl::bool_<false> >
class populate_python_funcs
{
public:
    template<class Graph>
    void operator()(const Graph& g, Descriptor& u, const dynamic_properties& dp, python::object& variables)
    {
        if (!init_counter<Graph>::registered)
        {
            base_from_tuple<typename PythonVertex<Graph>::base_t>();
            python::to_python_converter<typename PythonVertex<Graph>::base_t, base_to_tuple<typename PythonVertex<Graph>::base_t> >();
            base_from_tuple<typename PythonEdge<Graph>::base_t>();
            python::to_python_converter<typename PythonEdge<Graph>::base_t, base_to_tuple<typename PythonEdge<Graph>::base_t> >();
            init_counter<Graph>::registered = true;
        }

        variables["Vertex"] = PythonVertex<Graph>().GetPythonClass();
        variables["Edge"] = PythonEdge<Graph>().GetPythonClass();
        variables["Graph"] = PythonGraph<Graph>().GetPythonClass();
        populate_specific(g, u, dp, variables);
    }

private:

    template <class Base>
    struct base_to_tuple
    {
        static PyObject* convert(const Base& b)
        {
            boost::python::tuple t = boost::python::make_tuple(python::long_(size_t(&get<0>(b))), python::long_(size_t(&get<1>(b))), python::long_(size_t(&get<2>(b))), get<3>(b));
            return python::incref(t.ptr());
        }
    };

    template <class Base>
    struct base_from_tuple
    {
        base_from_tuple()
        {
            python::converter::registry::push_back(&convertible, &construct, boost::python::type_id<Base>());
        }

        static void* convertible(PyObject* obj_ptr)
        {
            python::handle<> x(python::borrowed(obj_ptr));
            python::object o(x);
            python::extract<size_t> first(o[0]);
            python::extract<size_t> second(o[1]);
            python::extract<size_t> third(o[2]);
            python::extract<bool> fourth(o[3]);
            if (!first.check() || !second.check() || !third.check() || !fourth.check())
                return 0;
            return obj_ptr;
        }
        
        static void construct(PyObject* obj_ptr, python::converter::rvalue_from_python_stage1_data* data)
        {          
            python::handle<> x(python::borrowed(obj_ptr));
            python::object o(x);
            typedef typename remove_reference<typename tuples::element<0,Base>::type>::type t0;
            typedef typename remove_reference<typename tuples::element<1,Base>::type>::type t1;
            typedef typename remove_reference<typename tuples::element<2,Base>::type>::type t2;
            t0* p0 = static_cast<t0*>((void *) size_t(python::extract<size_t>(o[0])));
            t1* p1 = static_cast<t1*>((void *) size_t(python::extract<size_t>(o[1])));
            t2* p2 = static_cast<t2*>((void *) size_t(python::extract<size_t>(o[2])));
            bool p4 = python::extract<bool>(o[3]);
            Base value(*p0, *p1, *p2, p4);
            
            void* storage = ( (python::converter::rvalue_from_python_storage<Base>*) data)->storage.bytes;
            new (storage) Base(value);
            data->convertible = storage;
        }
    };


    template <class Graph>
    void populate_specific(const Graph& g, typename graph_traits<Graph>::vertex_descriptor& v, const dynamic_properties& dp, python::object& variables)
    {
        typename PythonVertex<Graph>::base_t vertex_base(g, v, dp, false);
        python::object vertex = (PythonVertex<Graph>().GetPythonClass())(python::object(vertex_base));                
        variables["v"] = vertex;
        typename PythonGraph<Graph>::base_t graph_base(g, typename graph_traits<Graph>::edge_descriptor(), dp, false);
        python::object graph = (PythonGraph<Graph>().GetPythonClass())(python::object(graph_base));                
        variables["g"] = graph;
        populate_base(g, v, dp, variables, HasBase());
    }


    template <class Graph>
    void populate_base(const Graph& g, typename graph_traits<Graph>::vertex_descriptor& v, const dynamic_properties& dp, python::object& variables,  mpl::bool_<true>) 
    {
        if (!init_counter<Graph>::registered)
        {
            base_from_tuple<typename PythonVertex<typename Graph::graph_type>::base_t>();
            python::to_python_converter<typename PythonVertex<typename Graph::graph_type>::base_t, base_to_tuple<typename PythonVertex<typename Graph::graph_type>::base_t> >();
            base_from_tuple<typename PythonEdge<typename Graph::graph_type>::base_t>();
            python::to_python_converter<typename PythonEdge<typename Graph::graph_type>::base_t, base_to_tuple<typename PythonEdge<typename Graph::graph_type>::base_t> >();
            init_counter<Graph>::registered = true;
        }

        variables["BaseVertex"] = PythonVertex<typename Graph::graph_type>().GetPythonClass();
        variables["BaseEdge"] = PythonEdge<typename Graph::graph_type>().GetPythonClass();
        typename PythonVertex<typename Graph::graph_type>::base_t vertex_base(g.m_g, v, dp, false);
        python::object vertex = (PythonVertex<typename Graph::graph_type>().GetPythonClass())(python::object(vertex_base));
        variables["base_v"] = vertex;
        typename PythonGraph<typename Graph::graph_type>::base_t graph_base(g.m_g, typename graph_traits<typename Graph::graph_type>::edge_descriptor(), dp, false);
        python::object graph = (PythonGraph<typename Graph::graph_type>().GetPythonClass())(python::object(graph_base));
        variables["base_g"] = graph;
    }

    template <class Graph>
    void populate_base(const Graph& g, typename graph_traits<Graph>::vertex_descriptor& v, const dynamic_properties& dp, python::object& variables,  mpl::bool_<false>)
    {
    }

    template <class Graph>
    void populate_specific(const Graph& g, typename graph_traits<Graph>::edge_descriptor& e, const dynamic_properties& dp, python::object& variables)
    {
        typename PythonEdge<Graph>::base_t edge_base(g, e, dp, false);
        python::object edge = (PythonEdge<Graph>().GetPythonClass())(python::object(edge_base));                
        variables["e"] = edge;
        typename PythonGraph<Graph>::base_t graph_base(g, typename graph_traits<Graph>::edge_descriptor(), dp, false);
        python::object graph = (PythonGraph<Graph>().GetPythonClass())(python::object(graph_base));                
        variables["g"] = graph;
    }

    template <class Graph>
    void populate_specific(const Graph& g, graph_property_tag, const dynamic_properties& dp, python::object& variables)
    {
        typename PythonGraph<Graph>::base_t graph_base(g, typename graph_traits<Graph>::edge_descriptor(), dp, false);
        python::object graph = (PythonGraph<Graph>().GetPythonClass())(python::object(graph_base));                
        variables["g"] = graph;
    }

};


//==============================================================================
// PythonFilter
//==============================================================================
template <class Graph, class Descriptor, class HasBase = mpl::bool_<false> >
class PythonFilter
{
public:
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

    typedef mpl::vector<in_degreeS,out_degreeS,total_degreeS> degrees;

    PythonFilter(){}
    PythonFilter(const Graph& g, const dynamic_properties& dp, python::object filter)
        : _g(&g), _filter(filter[0])
    {
        python::object variables = filter[1];
        populate_python_funcs<Descriptor, HasBase>()(*_g, _u, dp, variables);
    }
    
 
    inline bool operator() (Descriptor u) const
    {              
        _u = u;
        return python::extract<bool>(_filter());
    }

private:
    Graph const*  _g;
    python::object _filter;
    static Descriptor _u;
};

template <class Graph, class Descriptor, class HasBase> 
Descriptor PythonFilter<Graph,Descriptor,HasBase>::_u;

} //graph_tool namespace

#endif
