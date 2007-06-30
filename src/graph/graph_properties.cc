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
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#include <boost/lambda/bind.hpp>

#include "graph.hh"
#include "histogram.hh"
#include "graph_filtering.hh"
#include "graph_python_filtering.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"

#include <boost/mpl/for_each.hpp>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

namespace graph_tool
{
// global property types
const char* type_names[] = {"boolean", "int", "long", "long", "float", "double", "string"};

// scalar types
const char* scalar_names[] = {"boolean", "int", "long", "long", "float", "double"};
}

//==============================================================================
// find_property_map(dp,name,key_type)
//==============================================================================
dynamic_property_map& 
graph_tool::find_property_map(const dynamic_properties& dp, string name, const type_info& key_type)
{
    for(typeof(dp.begin()) iter = dp.begin(); iter != dp.end(); ++iter)
        if (iter->first == name && iter->second->key() == key_type)
            return *iter->second;

    throw property_not_found(name);
}

//==============================================================================
// RemoveVertexProperty(property)
//==============================================================================
void GraphInterface::RemoveVertexProperty(string property)
{
    dynamic_properties_copy dp;
    try
    {
        dynamic_property_map& prop_map = find_property_map(_properties, property, typeid(graph_traits<multigraph_t>::vertex_descriptor));
        for (typeof(_properties.begin()) iter = _properties.begin(); iter != _properties.end(); ++iter)
        {
            if (iter->second != &prop_map)
                dp.insert(iter->first, auto_ptr<dynamic_property_map>(iter->second));
        }
    }
    catch (property_not_found)
    {
        throw GraphException("property '" + property + "' not found");
    }
    _properties = dp;
}

//==============================================================================
// RemoveEdgeProperty(property)
//==============================================================================
void GraphInterface::RemoveEdgeProperty(string property)
{
    dynamic_properties_copy dp;
    try
    {
        dynamic_property_map& prop_map = find_property_map(_properties, property, typeid(graph_traits<multigraph_t>::edge_descriptor));
        for (typeof(_properties.begin()) iter = _properties.begin(); iter != _properties.end(); ++iter)
        {
            if (iter->second != &prop_map)
                dp.insert(iter->first, auto_ptr<dynamic_property_map>(iter->second));
        }
    }
    catch (property_not_found)
    {
        throw GraphException("property '" + property + "' not found");
    }
    _properties = dp;
}

//==============================================================================
// RemoveGraphProperty(property)
//==============================================================================
void GraphInterface::RemoveGraphProperty(string property)
{
    dynamic_properties_copy dp;
    try
    {
        dynamic_property_map& prop_map = find_property_map(_properties, property, typeid(graph_property_tag));
        for (typeof(_properties.begin()) iter = _properties.begin(); iter != _properties.end(); ++iter)
        {
            if (iter->second != &prop_map)
                dp.insert(iter->first, auto_ptr<dynamic_property_map>(iter->second));
        }
    }
    catch (property_not_found)
    {
        throw GraphException("property '" + property + "' not found");
    }
    _properties = dp;
}


//==============================================================================
// edit_property
//==============================================================================
template <class Descriptor>
struct edit_property
{
    template <class Graph>
    void operator()(const Graph& g, const dynamic_properties& dp, dynamic_property_map* prop_map, python::object& op) const
    {
        typedef mpl::vector<in_degreeS,out_degreeS,total_degreeS> degrees;

        python::object operation = op[0], variables = op[1];

        typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
        typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

        typedef typename mpl::if_< is_same<Descriptor,vertex_descriptor>,
                                   vertex_descriptor, 
                                   typename mpl::if_<is_same<Descriptor,graph_property_tag>,
                                                     graph_property_tag, 
                                                     edge_descriptor>::type >::type descriptor_t;
        descriptor_t u;
        
        populate_python_funcs<descriptor_t>()(g, u, dp, variables);

        put_properties(g, u, *prop_map, operation);
    }
    
    template<class Graph>
    void put_properties(const Graph& g, typename graph_traits<Graph>::vertex_descriptor& v, 
                        dynamic_property_map& prop_map, python::object& operation) const
    {
        typename graph_traits<Graph>::vertex_iterator vi,v_end;
        for (tie(vi, v_end) = vertices(g); vi != v_end; ++vi)
        {
            v = *vi;
            python::object val = operation();
            prop_map.put(*vi, val);
        }
    }

    template<class Graph>
    void put_properties(const Graph& g, typename graph_traits<Graph>::edge_descriptor& e, 
                        dynamic_property_map& prop_map, python::object& operation) const
    {
        typename graph_traits<Graph>::edge_iterator ei,e_end;
        for (tie(ei, e_end) = edges(g); ei != e_end; ++ei)
        {
            e = *ei;
            Descriptor& ec = e;
            python::object val = operation();
            prop_map.put(ec, val);
        }
    }

    template<class Graph>
    void put_properties(const Graph& g, graph_property_tag, 
                        dynamic_property_map& prop_map, python::object& operation) const
    {
        python::object val = operation();
        prop_map.put(graph_property_tag(), val);
    }

};

//==============================================================================
// get_property_map()
//==============================================================================

template <class ValueTypes, class Descriptor, class IndexMap>
class get_property_map
{
public:
    get_property_map(GraphInterface& gi, dynamic_properties& dp, IndexMap index_map, string property, string type, const char* types[], python::object op, dynamic_property_map*& pmap)
        : _gi(gi), _dp(dp), _index_map(index_map), _property(property), _type(type), _types(types), _op(op), _pmap(pmap) {}
    
    template <class ValueType>
    class python_dynamic_property_map: public dynamic_property_map
    {
    public:
        python_dynamic_property_map(dynamic_property_map& dmap): _dmap(dmap) {}

        virtual void put(const any& key, const any& val) 
        {
            const python::object& o = any_cast<python::object>(val);
            ValueType value = python::extract<ValueType>(o); 
            _dmap.put(key, value); 
        }

        virtual any get(const any& key) { return _dmap.get(key); }
        virtual string get_string(const any& key) { return _dmap.get_string(key); }
        virtual const std::type_info& key() const { return _dmap.key(); } 
        virtual const std::type_info& value() const { return _dmap.value(); }
    private:
        dynamic_property_map& _dmap;
    };

    template <class ValueType>
    void operator()(ValueType)
    {
        if (_type == _types[mpl::find<ValueTypes,ValueType>::type::pos::value])
        {
            try
            {
                dynamic_property_map& pmap = find_property_map(_dp, _property, typeid(Descriptor));
                if (pmap.value() != typeid(ValueType))
                    throw GraphException("property \""+ _property + "\" already exists with a type other than " + _type + 
                                         ". Remove it first, or use the same type when editing.");
                _pmap = new python_dynamic_property_map<ValueType>(pmap);
            }
            catch (property_not_found)
            {
                typedef vector_property_map<ValueType, IndexMap> prop_map_t;
                prop_map_t prop_map(_index_map);
                _dp.property(_property, prop_map);
                _pmap = new python_dynamic_property_map<ValueType>(find_property_map(_dp, _property, typeid(typename IndexMap::key_type)));
            }
        }
    }

private:
    GraphInterface& _gi;
    dynamic_properties& _dp;
    IndexMap _index_map;
    string _property;
    string _type;
    const char** _types;
    python::object _op;
    dynamic_property_map*& _pmap;
};


//==============================================================================
// EditVertexProperty()
//==============================================================================

void GraphInterface::EditVertexProperty(string property, string type, python::object op)
{
    bool valid = false;
    for(int i = 0; i < mpl::size<value_types>::type::value; ++i)
        if (type == type_names[i])
            valid = true;
    if (!valid)
        throw GraphException("invalid type: " + type);

    dynamic_property_map* pmap;
    typedef graph_traits<multigraph_t>::vertex_descriptor vertex_descriptor;
    mpl::for_each<value_types>(get_property_map<value_types,vertex_descriptor,vertex_index_map_t>(*this, _properties, _vertex_index, property, type, type_names, op, pmap));
    check_filter(*this, lambda::bind<void>(edit_property<vertex_descriptor>(), lambda::_1, var(_properties), var(pmap), var(op)),
                 reverse_check(), directed_check());
    delete pmap;
}

//==============================================================================
// EditEdgeProperty()
//==============================================================================

void GraphInterface::EditEdgeProperty(string property, string type, python::object op)
{
    bool valid = false;
    for(int i = 0; i < mpl::size<value_types>::type::value; ++i)
        if (type == type_names[i])
            valid = true;
    if (!valid)
        throw GraphException("invalid type: " + type);

    dynamic_property_map* pmap;
    typedef graph_traits<multigraph_t>::edge_descriptor edge_descriptor;
    mpl::for_each<value_types>(get_property_map<value_types,edge_descriptor,edge_index_map_t>(*this, _properties, _edge_index, property, type, type_names, op, pmap));
    check_filter(*this, lambda::bind<void>(edit_property<edge_descriptor>(), lambda::_1, var(_properties), var(pmap), var(op)),
                 reverse_check(), directed_check());
    delete pmap;
}

//==============================================================================
// EditGraphProperty()
//==============================================================================

void GraphInterface::EditGraphProperty(string property, string type, python::object op)
{
    bool valid = false;
    for(int i = 0; i < mpl::size<value_types>::type::value; ++i)
        if (type == type_names[i])
            valid = true;
    if (!valid)
        throw GraphException("invalid type: " + type);

    dynamic_property_map* pmap;
    ConstantPropertyMap<size_t,graph_property_tag> graph_index(0);
    mpl::for_each<value_types>(get_property_map<value_types,graph_property_tag,ConstantPropertyMap<size_t,graph_property_tag> >(*this, _properties, graph_index, property, type, type_names, op, pmap));
    check_filter(*this, lambda::bind<void>(edit_property<graph_property_tag>(), lambda::_1, var(_properties), var(pmap), var(op)),
                 reverse_check(), directed_check());
    delete pmap;
}


//==============================================================================
// ListProperties()
//==============================================================================

template <class ValueTypes>
class print_name
{
public:
    print_name(const type_info& type, const char* types[]): _type(type), _types(types) {}

    template <class ValueType>
    void operator()(ValueType)
    {
        if (_type == typeid(ValueType))
            cout << _types[mpl::find<ValueTypes,ValueType>::type::pos::value];
    }
private:
    const type_info& _type;
    const char** _types;
};

void GraphInterface::ListProperties() const
{
    for (typeof(_properties.begin()) p = _properties.begin(); p != _properties.end(); ++p)
    {
        cout << setw(15) << left << p->first << " " << setw(8) << left;
        if (p->second->key() == typeid(graph_traits<multigraph_t>::vertex_descriptor))
            cout << "(vertex)";
        else 
            if (p->second->key() == typeid(graph_traits<multigraph_t>::edge_descriptor))
                cout << "(edge)";
            else
                cout << "(graph)";
        cout << "  type: ";
        mpl::for_each<value_types>(print_name<value_types>(p->second->value(), type_names));
        cout << endl;
    }
}

