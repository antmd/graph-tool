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
// edit_property
//==============================================================================
template <class Descriptor>
struct edit_property
{
    template <class Graph, class PropertyMap>
    void operator()(const Graph& g, const dynamic_properties& dp, PropertyMap prop_map, python::object& op) const
    {
        typedef mpl::vector<in_degreeS,out_degreeS,total_degreeS> degrees;

        python::object operation = op[0], variables = op[1];

        typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
        typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

        typedef typename mpl::if_<is_same<Descriptor,vertex_descriptor>,vertex_descriptor,edge_descriptor>::type descriptor_t;
        descriptor_t u;
        
        populate_python_funcs<descriptor_t>()(g, u, dp, variables);

        put_properties(g, u, prop_map, operation);
    }
    
    template<class Graph, class PropertyMap>
    void put_properties(const Graph& g, typename graph_traits<Graph>::vertex_descriptor& v, 
                        PropertyMap prop_map, python::object& operation) const
    {
        typename graph_traits<Graph>::vertex_iterator vi,v_end;
        for (tie(vi, v_end) = vertices(g); vi != v_end; ++vi)
        {
            v = *vi;
            typename property_traits<PropertyMap>::value_type val = python::extract<typename property_traits<PropertyMap>::value_type>(operation());
            put(prop_map, v, val);
        }
    }

    template<class Graph, class PropertyMap>
    void put_properties(const Graph& g, typename graph_traits<Graph>::edge_descriptor& e, 
                        PropertyMap prop_map, python::object& operation) const
    {
        typename graph_traits<Graph>::edge_iterator ei,e_end;
        for (tie(ei, e_end) = edges(g); ei != e_end; ++ei)
        {
            e = *ei;
            typename property_traits<PropertyMap>::value_type val = python::extract<typename property_traits<PropertyMap>::value_type>(operation());
            put(prop_map, e, val);
        }
    }
};

//==============================================================================
// update_property_map()
//==============================================================================

template <class ValueTypes, class Descriptor, class IndexMap>
class update_property_map
{
public:
    update_property_map(GraphInterface& gi, dynamic_properties& dp, IndexMap index_map, string property, string type, char* types[], python::object op)
        : _gi(gi), _dp(dp), _index_map(index_map), _property(property), _type(type), _types(types), _op(op) {}

    template <class ValueType>
    void operator()(ValueType)
    {
        if (_type == _types[mpl::find<ValueTypes,ValueType>::type::pos::value])
        {
            try
            {
                dynamic_property_map& dpmap = find_property_map(_dp, _property, typeid(Descriptor));
                typedef DynamicPropertyMapWrap<ValueType,Descriptor> prop_map_t;
                prop_map_t prop_map(dpmap);
                check_filter(_gi, bind<void>(edit_property<Descriptor>(), _1, var(_dp), prop_map, var(_op)),
                             reverse_check(), directed_check());
            }
            catch (property_not_found)
            {
                typedef vector_property_map<ValueType, IndexMap> prop_map_t;
                prop_map_t prop_map(_index_map);
                check_filter(_gi, bind<void>(edit_property<Descriptor>(), _1, var(_dp), prop_map, var(_op)),
                             reverse_check(), directed_check());
                _dp.property(_property, prop_map);
            }    

        }
    }

private:
    GraphInterface& _gi;
    dynamic_properties& _dp;
    IndexMap _index_map;
    string _property;
    string _type;
    char** _types;
    python::object _op;
};


//==============================================================================
// EditVertexProperty()
//==============================================================================

void GraphInterface::EditVertexProperty(string property, string type, python::object op)
{
    typedef mpl::vector<bool, int, long, float, double, std::string> value_types;
    char* type_names[] = {"boolean", "int", "long", "float", "double", "string"};

    bool valid = false;
    for(int i = 0; i < mpl::size<value_types>::type::value; ++i)
        if (type == type_names[i])
            valid = true;
    if (!valid)
        throw GraphException("invalid type: " + type);

    typedef graph_traits<multigraph_t>::vertex_descriptor vertex_descriptor;
    mpl::for_each<value_types>(update_property_map<value_types,vertex_descriptor,vertex_index_map_t>(*this, _properties, _vertex_index, property, type, type_names, op));
}

//==============================================================================
// EditEdgeProperty()
//==============================================================================

void GraphInterface::EditEdgeProperty(string property, string type, python::object op)
{
    typedef mpl::vector<bool, int, long, float, double, std::string> value_types;
    char* type_names[] = {"boolean", "int", "long", "float", "double", "string"};

    bool valid = false;
    for(int i = 0; i < mpl::size<value_types>::type::value; ++i)
        if (type == type_names[i])
            valid = true;
    if (!valid)
        throw GraphException("invalid type: " + type);

    typedef graph_traits<multigraph_t>::edge_descriptor edge_descriptor;
    mpl::for_each<value_types>(update_property_map<value_types,edge_descriptor,edge_index_map_t>(*this, _properties, _edge_index, property, type, type_names, op));
}

//==============================================================================
// ListProperties()
//==============================================================================

template <class ValueTypes>
class print_name
{
public:
    print_name(const type_info& type, char* types[]): _type(type), _types(types) {}

    template <class ValueType>
    void operator()(ValueType)
    {
        if (_type == typeid(ValueType))
            cout << _types[mpl::find<ValueTypes,ValueType>::type::pos::value];
    }
private:
    const type_info& _type;
    char** _types;
};

void GraphInterface::ListProperties() const
{
    typedef mpl::vector<bool, int, long, float, double, std::string> value_types;
    char* type_names[] = {"boolean", "int", "long", "float", "double", "string"};

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

