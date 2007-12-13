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


#include <boost/lambda/bind.hpp>

#include "graph.hh"
#include "histogram.hh"
#include "graph_filtering.hh"
#include "graph_python_interface.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"

#include <boost/mpl/for_each.hpp>
#include <iostream>
#include <iomanip>
#include <boost/ref.hpp>

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

namespace graph_tool
{

// pos_t i/o
std::ostream& operator<<(std::ostream &o, const pos_t &p )
{
    o << p.x << "," << p.y; return o;
}
std::istream& operator>>(std::istream &o, pos_t &p )
{
    char c; o >> p.x >> c >> p.y; return o;
}

// global property types
const char* type_names[] =
    {"boolean", "int", "long", "size_t", "float", "double", "string", "pos_t"};

// scalar types
const char* scalar_names[] =
    {"boolean", "int", "long", "size_t", "float", "double"};
}

// this function gets the dynamic property map inside dp which matches the given
// name and key type
dynamic_property_map&
graph_tool::find_property_map(const dynamic_properties& dp, string name,
                              const type_info& key_type)
{
    for(typeof(dp.begin()) iter = dp.begin(); iter != dp.end(); ++iter)
        if (iter->first == name && iter->second->key() == key_type)
            return *iter->second;

    throw property_not_found(name);
}

void GraphInterface::RemoveVertexProperty(string property)
{
    dynamic_properties_copy dp;
    try
    {
        dynamic_property_map& prop_map =
            find_property_map(_properties, property, typeid(vertex_t));
        for (typeof(_properties.begin()) iter = _properties.begin();
             iter != _properties.end(); ++iter)
        {
            if (iter->second != &prop_map)
                dp.insert(iter->first,
                          auto_ptr<dynamic_property_map>(iter->second));
        }
    }
    catch (property_not_found)
    {
        throw GraphException("property '" + property + "' not found");
    }
    _properties = dp;
}

void GraphInterface::RemoveEdgeProperty(string property)
{
    dynamic_properties_copy dp;
    try
    {
        dynamic_property_map& prop_map =
            find_property_map(_properties, property, typeid(edge_t));
        for (typeof(_properties.begin()) iter = _properties.begin();
             iter != _properties.end(); ++iter)
        {
            if (iter->second != &prop_map)
                dp.insert(iter->first,
                          auto_ptr<dynamic_property_map>(iter->second));
        }
    }
    catch (property_not_found)
    {
        throw GraphException("property '" + property + "' not found");
    }
    _properties = dp;
}

void GraphInterface::RemoveGraphProperty(string property)
{
    dynamic_properties_copy dp;
    try
    {
        dynamic_property_map& prop_map =
            find_property_map(_properties, property,
                              typeid(graph_property_tag));
        for (typeof(_properties.begin()) iter = _properties.begin();
             iter != _properties.end(); ++iter)
        {
            if (iter->second != &prop_map)
                dp.insert(iter->first,
                          auto_ptr<dynamic_property_map>(iter->second));
        }
    }
    catch (property_not_found)
    {
        throw GraphException("property '" + property + "' not found");
    }
    _properties = dp;
}



// the following will edit a specific property map, given its python edit
// function

template <class Descriptor>
struct edit_property
{
    template <class Graph>
    void operator()(const Graph& g, const dynamic_properties& dp,
                    dynamic_property_map* prop_map, python::object& op,
                    python::object& pg) const
    {
        put_properties(g, Descriptor(), *prop_map, op, pg);
    }

    template<class Graph>
    void put_properties(const Graph& g,
                        typename graph_traits<Graph>::vertex_descriptor,
                        dynamic_property_map& prop_map,
                        python::object& operation, python::object& pg) const
    {
        if (operation == python::object()) // don't set properties if op == None
            return;

        typename graph_traits<Graph>::vertex_iterator v,v_end;
        for (tie(v, v_end) = vertices(g); v != v_end; ++v)
        {
            python::object val =
                operation(python::object(PythonVertex<Graph>(g, *v)), pg);
            prop_map.put(*v, val);
        }
    }

    template<class Graph>
    void put_properties(const Graph& g,
                        typename graph_traits<Graph>::edge_descriptor,
                        dynamic_property_map& prop_map,
                        python::object& operation, python::object& pg) const
    {
        if (operation == python::object()) // don't set properties if op == None
            return;

        typename graph_traits<Graph>::edge_iterator e, e_end;
        for (tie(e, e_end) = edges(g); e != e_end; ++e)
        {
            python::object val =
                operation(python::object(PythonEdge<Graph>(g, *e)), pg);
            prop_map.put(*e, val);
        }
    }
};


// the function below will get a dynamic_property_map (or create a new one if it
// doesn't exist, based on a vector property map) with the type given by the
// "type" string argument. It should be used in conjunction with mpl::for_each,
// as is done below
template <class ValueTypes, class Descriptor, class IndexMap>
class get_property_map
{
public:
    get_property_map(GraphInterface& gi, dynamic_properties& dp,
                     IndexMap index_map, string property, string type,
                     const char* types[], dynamic_property_map*& pmap)
        : _gi(gi), _dp(dp), _index_map(index_map), _property(property),
          _type(type), _types(types), _pmap(pmap) {}

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
        virtual string get_string(const any& key)
        {
            return _dmap.get_string(key);
        }
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
                dynamic_property_map& pmap =
                    find_property_map(_dp, _property, typeid(Descriptor));
                if (pmap.value() != typeid(ValueType))
                    throw GraphException("property \""+ _property +
                                         "\" already exists with a type "
                                         "other than " + _type +
                                         ". Remove it first, or use the same "
                                         "type when editing.");
                _pmap = new python_dynamic_property_map<ValueType>(pmap);
            }
            catch (property_not_found)
            {
                typedef vector_property_map<ValueType, IndexMap> prop_map_t;
                prop_map_t prop_map(_index_map);
                _dp.property(_property, prop_map);
                _pmap = new python_dynamic_property_map<ValueType>
                    (find_property_map(_dp, _property,
                                       typeid(typename IndexMap::key_type)));
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
    dynamic_property_map*& _pmap;
};


void GraphInterface::EditVertexProperty(string property, string type, 
                                        python::object op, python::object g)
{
    bool valid = false;
    for(int i = 0; i < mpl::size<value_types>::type::value; ++i)
        if (type == type_names[i])
            valid = true;
    if (!valid)
        throw GraphException("invalid type: " + type);

    dynamic_property_map* pmap;
    typedef vertex_t vertex_descriptor;
    mpl::for_each<value_types>
        (get_property_map<value_types,vertex_descriptor,vertex_index_map_t>
         (*this, _properties, _vertex_index, property, type, type_names, pmap));
    run_action(*this, lambda::bind<void>(edit_property<vertex_descriptor>(),
                                         lambda::_1, var(_properties),
                                         var(pmap), var(op), var(g)));
    delete pmap;
}


void GraphInterface::EditEdgeProperty(string property, string type,
                                      python::object op, python::object g)
{
    bool valid = false;
    for(int i = 0; i < mpl::size<value_types>::type::value; ++i)
        if (type == type_names[i])
            valid = true;
    if (!valid)
        throw GraphException("invalid type: " + type);

    dynamic_property_map* pmap;
    typedef edge_t edge_descriptor;
    mpl::for_each<value_types>
        (get_property_map<value_types,edge_descriptor,edge_index_map_t>
         (*this, _properties, _edge_index, property, type, type_names,pmap));
    run_action(*this, lambda::bind<void>(edit_property<edge_descriptor>(),
                                         lambda::_1, var(_properties), 
                                         var(pmap), var(op), var(g)));
    delete pmap;
}

void GraphInterface::EditGraphProperty(string property, string type,
                                       python::object op, python::object g)
{
    bool valid = false;
    for(int i = 0; i < mpl::size<value_types>::type::value; ++i)
        if (type == type_names[i])
            valid = true;
    if (!valid)
        throw GraphException("invalid type: " + type);

    dynamic_property_map* pmap;
    ConstantPropertyMap<size_t,graph_property_tag> graph_index(0);
    mpl::for_each<value_types>
        (get_property_map
         <value_types,graph_property_tag,
          ConstantPropertyMap<size_t,graph_property_tag> >
         (*this, _properties, graph_index, property, type, type_names, pmap));

    if (op != python::object()) // don't set property if op == None
    {
        python::object val = op(g);
        pmap->put(graph_property_tag(), val);
    }
    delete pmap;
}

// this will return the name of a given type
template <class ValueTypes>
class get_type_name
{
public:
    get_type_name(const type_info& type, const char* types[], string& name)
        : _type(type), _types(types), _name(name) {}

    template <class ValueType>
    void operator()(ValueType)
    {
        if (_type == typeid(ValueType))
            _name = _types[mpl::find<ValueTypes,ValueType>::type::pos::value];
    }
private:
    const type_info& _type;
    const char** _types;
    string& _name;
};

