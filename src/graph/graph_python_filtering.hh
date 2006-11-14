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
#include <boost/python/dict.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/make_function.hpp>

namespace graph_tool
{
using namespace boost;

//==============================================================================
// populate_python_funcs
//==============================================================================

template <class Descriptor, class HasBase = mpl::bool_<false> >
struct populate_python_funcs
{
    template<class Graph>
    void operator()(const Graph& g, Descriptor& u, const dynamic_properties& dp, python::object& variables)
    {
        for(typeof(dp.begin()) iter = dp.begin(); iter != dp.end(); ++iter)
        {
            if (iter->second->key() == typeid(Descriptor))
        	variables[iter->first] = python::make_function(get_value<Descriptor>(*iter->second, u), 
        						       python::default_call_policies(), mpl::vector<python::object>::type());
        }
        populate_specific(g, u, dp, variables);
    }

    typedef mpl::vector<in_degreeS, out_degreeS, total_degreeS> degrees;

    template <class Graph>
    void populate_specific(const Graph& g, typename graph_traits<Graph>::vertex_descriptor& v, const dynamic_properties& dp, python::object& variables)
    {
        mpl::for_each<degrees>(put_degree_function<Graph>(g, v, variables));

        typedef typename mpl::if_<HasBase, degrees, mpl::vector<> >::type base_degrees;
        mpl::for_each<base_degrees>(put_base_degree_function<Graph>(g, v, variables, "orig_"));
    }

    template <class Graph>
    void populate_specific(const Graph& g, typename graph_traits<Graph>::edge_descriptor& e, const dynamic_properties& dp, python::object& variables)
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;

        variables["is_self"] = python::make_function(is_self<Graph>(g, e), python::default_call_policies(), mpl::vector<python::object>::type()); 
        variables["n_parallel"] = python::make_function(n_parallel<Graph>(g, e), python::default_call_policies(), mpl::vector<python::object>::type()); 
        
        for(typeof(dp.begin()) iter = dp.begin(); iter != dp.end(); ++iter)
        {
            if (iter->second->key() == typeid(vertex_descriptor))
            {
        	variables["source_"+iter->first] = python::make_function(get_source_or_target_value<Graph,true>(g, *iter->second, e), 
        								 python::default_call_policies(), mpl::vector<python::object>::type());
        	variables["target_"+iter->first] = python::make_function(get_source_or_target_value<Graph,false>(g, *iter->second, e), 
        								 python::default_call_policies(), mpl::vector<python::object>::type());
            }
        }
        mpl::for_each<degrees>(put_source_or_target_degree_function<Graph,true>(g, e, variables, "source_"));
        mpl::for_each<degrees>(put_source_or_target_degree_function<Graph,false>(g, e, variables, "target_"));
    }


    template <class VertexOrEdge>
    struct get_value
    {
        get_value(const dynamic_property_map& dmap, const VertexOrEdge& e)
            : _dmap(dmap), _e(e) {}
        
        struct try_conversion
        {
            try_conversion(get_value& parent): _parent(parent) {}

            template <class Type>
            void operator()(Type)
            {
        	any any_val = const_cast<dynamic_property_map&>(_parent._dmap).get(_parent._e);
        	Type* value = any_cast<Type>(&any_val);
        	if (value != 0)
        	    _parent._retval = python::object(*value);
            }
            
            get_value& _parent;
        };
        
        python::object operator()()
        {
            typedef mpl::vector<bool,int,long,size_t,float,double,std::string> value_types;
            mpl::for_each<value_types>(try_conversion(*this));
            return _retval;
        }

        const dynamic_property_map& _dmap;
        const VertexOrEdge& _e;
        python::object _retval;
    };

    template <class Graph, bool Source>
    struct get_source_or_target_value
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
        typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

        get_source_or_target_value(const Graph& g, dynamic_property_map& dmap, const edge_descriptor& e)
            : _g(g),_dmap(dmap),_e(e){}

        python::object operator()()
        {
            vertex_descriptor _s;

            if (Source)
        	_s = source(_e, _g);
            else
        	_s = target(_e, _g);

            get_value<vertex_descriptor> get_value(_dmap, _s);
            return get_value();
        }

        const Graph& _g;	
        const dynamic_property_map& _dmap;
        const edge_descriptor& _e;	
    };

    template <class Graph, class Degree>
    struct get_degree
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;

        get_degree(const Graph& g, const vertex_descriptor& v)
            : _g(g), _v(v) {}
        	
        python::object operator()()
        {
            return python::object(_degree(_v, _g));
        }

        const Graph& _g;
        const vertex_descriptor& _v;
        Degree _degree;
    };

    template <class Graph>
    struct put_degree_function
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;

        put_degree_function(const Graph& g, const vertex_descriptor& v, python::object variables, std::string prefix = "")
            : _g(g), _v(v), _variables(variables), _prefix(prefix) {}
        template <class Degree>
        void operator()(Degree degree)
        {
            _variables[_prefix+degree.name()] =  python::make_function(get_degree<Graph,Degree>(_g, _v),
        							       python::default_call_policies(), mpl::vector<python::object>::type());
        }
        const Graph& _g;
        const vertex_descriptor& _v;
        python::object& _variables;
        std::string _prefix;
    };

    template <class Graph>    
    struct put_base_degree_function: public put_degree_function<Graph>
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;

        put_base_degree_function(const Graph& g, const vertex_descriptor& v, python::object& variables, std::string prefix = "")
            : put_degree_function<Graph>(g, v, variables, prefix) {}

        template <class Degree>
        void operator()(Degree degree)
        {
            this->_variables[this->_prefix+degree.name()] = python::make_function(get_degree<typename Graph::graph_type,Degree>(this->_g.m_g, this->_v),
        									  python::default_call_policies(), mpl::vector<python::object>::type());
        }
    };

    template <class Graph, bool Source>
    struct put_source_or_target_degree_function
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
        typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

        put_source_or_target_degree_function(const Graph& g, const edge_descriptor& e, python::object& variables, std::string prefix = "")
            : _g(g), _e(e), _variables(variables), _prefix(prefix) {}
        template <class Degree>
        void operator()(Degree d)
        {
            vertex_descriptor v;	    
            if (Source)
        	v = source(_e, this->_g);
            else
        	v = target(_e, this->_g);
            put_degree_function<Graph>(_g, v, _variables, _prefix)(d);
        }

        const Graph& _g;
        const edge_descriptor& _e;
        python::object& _variables;
        std::string _prefix;
    };

    template<class Graph>
    struct is_self
    {
        typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

        is_self(const Graph& g, const edge_descriptor& e): _g(g), _e(e) {}
        python::object operator()()
        {
            return python::object(source(_e, _g) == target(_e, _g));
        }
        const Graph& _g;
        const edge_descriptor& _e;
    };

    template<class Graph>
    struct n_parallel
    {
        typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

        n_parallel(const Graph& g, const edge_descriptor& e): _g(g), _e(e) {}
        python::object operator()()
        {
            size_t n = 0;
            typename graph_traits<Graph>::vertex_descriptor s,t;
            s = source(_e, _g);
            t = target(_e, _g);
            typename graph_traits<Graph>::adjacency_iterator a, a_end;
            for(tie(a, a_end) = adjacent_vertices(s, _g); a != a_end; ++a)
        	if (*a == t)
        	    n++;
            return python::object(n-1);
        }
        const Graph& _g;
        const edge_descriptor& _e;
    };
    

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
