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

#include "graph_properties.hh"
#include "graph_util.hh"
#include "graph.hh"

#include <boost/lambda/bind.hpp>
#include <boost/mpl/for_each.hpp>

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

namespace graph_tool
{

// global property types
const char* type_names[] =
    {"bool", "int32_t", "int64_t", "double", "long double",
     "string", "vector<bool>","vector<int32_t>", "vector<int64_t>",
     "vector<double>", "vector<long double>"};

// this function gets the dynamic property map inside dp which matches the given
// name and key type
dynamic_property_map&
find_property_map(const dynamic_properties& dp, string name,
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


struct add_property_map
{
    template <class ValueType, class IndexMap>
    void operator()(ValueType, const string& type_name, IndexMap index,
                    const string& property_name, dynamic_properties& dp,
                    bool& found) const
    {
        size_t i = mpl::find<value_types,ValueType>::type::pos::value;
        if (type_name == type_names[i])
        {
            vector_property_map<ValueType, IndexMap> prop(index);
            dp.property(property_name, prop);
            found = true;
        }
    }
};

void GraphInterface::AddVertexProperty(string property, string type)
{
    try
    {
        find_property_map(_properties, property, typeid(vertex_t));
    }
    catch (property_not_found)
    {
        throw GraphException("A vertex property named " + property +
                             " already exists");
    }
    bool found = false;
    mpl::for_each<value_types>(lambda::bind<void>(add_property_map(),
                                                  lambda::_1, lambda::var(type),
                                                  _vertex_index,
                                                  lambda::var(property),
                                                  lambda::var(_properties), 
                                                  lambda::var(found)));
    if (!found)
        throw GraphException("Invalid property type " + type);
}

void GraphInterface::AddEdgeProperty(string property, string type)
{
    try
    {
        find_property_map(_properties, property, typeid(edge_t));
    }
    catch (property_not_found)
    {
        throw GraphException("An edge property named " + property +
                             " already exists");
    }
    bool found = false;
    mpl::for_each<value_types>(lambda::bind<void>(add_property_map(), 
                                                  lambda::_1, lambda::var(type),
                                                  _edge_index, 
                                                  lambda::var(property),
                                                  lambda::var(_properties),
                                                  lambda::var(found)));
    if (!found)
        throw GraphException("Invalid property type " + type);
}

void GraphInterface::AddGraphProperty(string property, string type)
{
    ConstantPropertyMap<size_t,graph_property_tag> graph_index(0);
    try
    {
        find_property_map(_properties, property, typeid(graph_property_tag));
    }
    catch (property_not_found)
    {
        throw GraphException("A graph property named " + property +
                             " already exists");
    }
    bool found = false;
    mpl::for_each<value_types>(lambda::bind<void>(add_property_map(), 
                                                  lambda::_1, lambda::var(type),
                                                  graph_index, 
                                                  lambda::var(property),
                                                  lambda::var(_properties),
                                                  lambda::var(found)));
    if (!found)
        throw GraphException("Invalid property type " + type);
}


} // graph_tool namespace

