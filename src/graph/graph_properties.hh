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

#ifndef GRAPH_PROPERTIES_HH
#define GRAPH_PROPERTIES_HH

#include <string>
#include <tr1/unordered_map>
#include <boost/property_map.hpp>
#include <boost/dynamic_property_map.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/functional/hash.hpp>
#include <boost/shared_ptr.hpp>

namespace graph_tool 
{

// global property types
typedef boost::mpl::vector<bool, int, long, size_t, float, double, std::string> value_types;
extern const char* type_names[];

// scalar types
typedef boost::mpl::vector<bool, int, long, size_t, float, double> scalar_types;
extern const char* scalar_names[];



//==============================================================================
// Property Map Utility Functions
//==============================================================================

//==============================================================================
// get_static_property_map(map)
// gets the "static" property map behind the dynamic property map, or throws
// bad_cast if it doesn't match
//==============================================================================

template <class PropertyMap>
PropertyMap& get_static_property_map(boost::dynamic_property_map& map)
{
    return dynamic_cast<boost::detail::dynamic_property_map_adaptor<PropertyMap>&>(map).base();
}

template <class PropertyMap>
PropertyMap* get_static_property_map(boost::dynamic_property_map* map)
{
    boost::detail::dynamic_property_map_adaptor<PropertyMap>* adaptor = 
        dynamic_cast<boost::detail::dynamic_property_map_adaptor<PropertyMap>*>(map);
    if (adaptor)
        return &adaptor->base();
    else
        return 0;
}


//==============================================================================
// find_property_map(dp, name, key_type)
// gets the dynamic property map inside dp which matches the given name and 
// key type
//==============================================================================

boost::dynamic_property_map& find_property_map(const boost::dynamic_properties& dp, std::string name, const std::type_info& key_type);

//==============================================================================
// dynamic_properties_copy
// contains a copy of a property map, which does not delete its members when it
// deconstructs
//==============================================================================

struct dynamic_properties_copy: public boost::dynamic_properties
{
    dynamic_properties_copy() {}
    dynamic_properties_copy(boost::dynamic_properties& dp): boost::dynamic_properties(dp) {}
    dynamic_properties_copy(const boost::function<std::auto_ptr<boost::dynamic_property_map>(const std::string&, const boost::any&, const boost::any&)>& fn)
        : boost::dynamic_properties(fn) {}
    ~dynamic_properties_copy()
    {        
        for (typeof(this->begin()) iter = this->begin(); iter != this->end(); ++iter)
            iter->second = 0; // will be deleted when original dp deconstructs
    }
};

//==============================================================================
// template<ConvertedType,Key>
// get_converted_scalar_value(map, key)
// gets the value in the map corresponding to key, converted to ConvertedType,
// or throws bad_lexical_cast. 
//==============================================================================

template <class T>
struct AttemptAnyConversion; // forward declaration

template <class ConvertedType, class Key>
ConvertedType get_converted_scalar_value(boost::dynamic_property_map& dmap, const Key& key)
{
    typedef typename boost::mpl::vector<long double, double, float, unsigned long long, long long, 
                                        unsigned long, long, unsigned int, int, unsigned short, short, 
                                        unsigned char, char, std::string>::type scalar_types;
    ConvertedType target;
    const boost::any& source = dmap.get(key);
    bool success;
    if (dmap.value() == typeid(ConvertedType))
    {
        target = boost::any_cast<ConvertedType>(dmap.get(key));
        success = true;
    }
    else
    {
        boost::mpl::for_each<scalar_types>(AttemptAnyConversion<ConvertedType>(target, source, success));
    }
    if (!success)
        throw boost::bad_lexical_cast();
    return target;
}

template <class T>
struct AttemptAnyConversion
{
    AttemptAnyConversion(T& value, const boost::any& source, bool& success)
        :_value(value), _source(source), _success(success) 
    { 
        _success = false;
    }

    template <class Source>
    void operator()(Source)
    {
        try
        {
            _value = boost::lexical_cast<T>(boost::any_cast<Source>(_source));
            _success = true;
        }
        catch (boost::bad_any_cast){}
        catch (boost::bad_lexical_cast){}
    }

    T& _value;
    const boost::any& _source;
    bool& _success;
};

//==============================================================================
// DynamicPropertyMapWrap
// wraps a dynamic_property_map, so it can be used as a regular property map
//==============================================================================

template <class Value, class Key>
class DynamicPropertyMapWrap
{
public:
    typedef Value value_type;
    typedef Value reference;
    typedef Key key_type;
    typedef boost::read_write_property_map_tag category;

    DynamicPropertyMapWrap(boost::dynamic_property_map& dmap):_dmap(&dmap) {}
    DynamicPropertyMapWrap() {}

    Value get(const Key& k) const
    {
        return get_converted_scalar_value<Value>(*_dmap, k);
    }

    void put(const Key& k, const Value& val)
    {
        _dmap->put(k, val);
    }

private:
    boost::dynamic_property_map* _dmap;
};

} // graph_tool namespace

namespace boost 
{

template <class Value, class Key>
Value get(const graph_tool::DynamicPropertyMapWrap<Value,Key>& pmap, typename property_traits<graph_tool::DynamicPropertyMapWrap<Value,Key> >::key_type k)
{
    return pmap.get(k);
}

template <class Value, class Key>
void put(graph_tool::DynamicPropertyMapWrap<Value,Key> pmap, typename property_traits<graph_tool::DynamicPropertyMapWrap<Value,Key> >::key_type k, 
         typename property_traits<graph_tool::DynamicPropertyMapWrap<Value,Key> >::value_type val)
{
    pmap.put(k,val);
}

} // boost namespace 

namespace graph_tool
{

//==============================================================================
// DescriptorHash
// property map based on a hashed container
//==============================================================================

template <class IndexMap>
class DescriptorHash: public std::unary_function<typename IndexMap::key_type, std::size_t> 
{
public:
    DescriptorHash() {}
    DescriptorHash(IndexMap index_map): _index_map(index_map) {}
    std::size_t operator()(typename IndexMap::key_type const& d) const { return boost::hash_value(_index_map[d]); }
private:
    IndexMap _index_map;
};

template <class IndexMap, class Value>
class HashedDescriptorMap
    : public boost::put_get_helper<Value&, HashedDescriptorMap<IndexMap,Value> >
{
public:
    typedef DescriptorHash<IndexMap> hashfc_t;
    typedef std::tr1::unordered_map<typename IndexMap::key_type,Value,hashfc_t> map_t;
    typedef boost::associative_property_map<map_t> prop_map_t;

    typedef typename boost::property_traits<prop_map_t>::value_type value_type;
    typedef typename boost::property_traits<prop_map_t>::reference reference;
    typedef typename boost::property_traits<prop_map_t>::key_type key_type;
    typedef typename boost::property_traits<prop_map_t>::category category;

    HashedDescriptorMap(IndexMap index_map): _base_map(new map_t(0, hashfc_t(index_map))), _prop_map(*_base_map) {}
    HashedDescriptorMap(){}
    
    reference operator[](const key_type& k) { return _prop_map[k]; }
    const reference operator[](const key_type& k) const { return _prop_map[k]; }

private:
    boost::shared_ptr<map_t> _base_map;
    prop_map_t _prop_map;
};


//==============================================================================
// InitializedPropertyMap
// this wraps a container as a property map which is automatically initialized
// with a given default value
//==============================================================================

template <class Container>
class InitializedPropertyMap
    : public boost::put_get_helper<typename Container::value_type::second_type&, InitializedPropertyMap<Container> >
{
public:
    typedef typename Container::value_type::second_type value_type;
    typedef value_type& reference;
    typedef typename Container::key_type key_type;
    typedef boost::read_write_property_map_tag category;

    InitializedPropertyMap(Container& base_map, value_type def)
        : _base_map(&base_map), _default(def) {}
    InitializedPropertyMap(){}

    reference operator[](const key_type& k)
    {
        return get(k);
    }

    const reference operator[](const key_type& k) const
    {
        return get(k);
    }

    const reference get(const key_type& k) const
    {
        typename Container::iterator val;
        val = _base_map->find(k);
        if (val == _base_map->end())
            val = _base_map->insert(make_pair(k, _default)).first;
        return val->second;
    }

private:
    Container* _base_map;
    value_type _default;
};

//==============================================================================
// ConstantPropertyMap
// a property map which returns a constant value
//==============================================================================

template <class Value, class Key>
class ConstantPropertyMap
    : public boost::put_get_helper<Value, ConstantPropertyMap<Value,Key> >
{
public:
    typedef Value value_type;
    typedef value_type& reference;
    typedef Key key_type;
    typedef boost::readable_property_map_tag category;

    ConstantPropertyMap(value_type c): _c(c) {}
    ConstantPropertyMap(){}

    const value_type& operator[](const key_type& k) const
    {
        return _c;
    }

private:
    value_type _c;
};


} // graph_tool namespace


#endif
