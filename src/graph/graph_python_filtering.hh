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
#include <boost/functional/hash.hpp>

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

template <class Graph>
class PythonEdge;

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
            _python_class.def("out_edges", &PythonVertex::GetOutEdges);
            _python_class.def("in_edges", &PythonVertex::GetInEdges);
            _python_class.def(python::self == python::self);
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
            return python::object();
        }
        else
        {
            return python::object();
        }
    }

    python::object GetOutEdges()
    {
        if (_base != python::object())
        {
            python::list out_edge_list;
            base_t base = python::extract<base_t>(_base);
            const Graph& g = get<0>(base);
            const vertex_descriptor& v = get<1>(base);
            const dynamic_properties& dp = get<2>(base);
            python::object edge_class = PythonEdge<Graph>().GetPythonClass();
            typename graph_traits<Graph>::out_edge_iterator e, e_end;
            for (tie(e, e_end) = out_edges(v, g); e != e_end; ++e)
            {
                edge_descriptor* new_edge = new edge_descriptor;
                *new_edge = *e;
                typename PythonEdge<Graph>::base_t edge_base(g, *new_edge, dp, true);
                python::object edge = (edge_class)(python::object(edge_base));
                out_edge_list.append(edge);
            }
            return out_edge_list;
        }
        else
        {
            return python::object();
        }
    }

    python::object GetInEdges()
    {
        if (_base != python::object())
        {
            base_t base = python::extract<base_t>(_base);
            const Graph& g = get<0>(base);
            const vertex_descriptor& v = get<1>(base);
            const dynamic_properties& dp = get<2>(base);
            return get_in_edges(g, v, dp, typename is_convertible<typename graph_traits<Graph>::directed_category, undirected_tag>::type());
        }
        else
        {
            return python::object();
        }
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

    python::class_<PythonVertex> GetPythonClass() { return _python_class; }
        
private:
    template<class G>
    python::object get_in_edges(const G& g, const typename graph_traits<G>::vertex_descriptor& v, const dynamic_properties& dp, true_type)
    {
        python::list in_edge_list;
        typename graph_traits<Graph>::in_edge_iterator e, e_end;
        for (tie(e, e_end) = in_edges(v, g); e != e_end; ++e)
        {
            edge_descriptor* new_edge = new edge_descriptor;
            *new_edge = *e;
            typename PythonEdge<Graph>::base_t edge_base(g, *new_edge, dp, true);
            python::object edge = (PythonEdge<Graph>().GetPythonClass())(python::object(edge_base));
            in_edge_list.append(edge);
        }
        return in_edge_list;
    }
        
    template<class G>
    python::object get_in_edges(const G& g, const typename graph_traits<G>::vertex_descriptor& v, const dynamic_properties& dp, false_type)
    {
        python::list props;
        return props;
    }

    python::object _base;
    static python::class_<PythonVertex> _python_class;
};

template <class Graph> python::class_<PythonVertex<Graph> > PythonVertex<Graph>::_python_class =  python::class_<PythonVertex<Graph> >("Vertex", python::no_init);


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
            _python_class.def("n_parallel", &PythonEdge::GetNParallel);
            _python_class.def(python::self == python::self);
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
                if (iter->second->key() != typeid(vertex_descriptor)) // FIXME: graph properties?
                    return get_python_property_value<edge_descriptor>(*iter->second, e)();
            }
            return python::object();
        }
        else
        {
            return python::object();
        }
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

    python::class_<PythonEdge> GetPythonClass() { return _python_class; }
        
private:

    python::object _base;
    static python::class_<PythonEdge> _python_class;
};

template <class Graph> python::class_<PythonEdge<Graph> > PythonEdge<Graph>::_python_class =  python::class_<PythonEdge<Graph> >("Edge", python::no_init);

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
