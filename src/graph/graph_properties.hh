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

#ifndef GRAPH_PROPERTIES_HH
#define GRAPH_PROPERTIES_HH

#include <string>
#include <vector>
#include <tr1/unordered_map>

#include <boost/property_map.hpp>
#include <boost/dynamic_property_map.hpp>
#include <boost/vector_property_map.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/find.hpp>
#include <boost/lambda/bind.hpp>

// this file provides general functions for manipulating graph properties

namespace graph_tool
{
using namespace std;
using namespace boost;

// global property types. only these types are allowed in property maps
typedef mpl::vector<bool, int32_t, int64_t, double, long double, string,
                    vector<bool>, vector<int32_t>, vector<int64_t>,
                    vector<double>, vector<long double> >
    value_types;

extern const char* type_names[]; // respective type names

// scalar types: types contained in value_types which are scalar
typedef mpl::vector<bool, int32_t, int64_t, double, long double> scalar_types;

// integer_types: scalar types which are integer
typedef mpl::vector<bool, int32_t, int64_t> integer_types;

// floating_types: scalar types which are floating point
typedef mpl::vector<double, long double> floating_types;

struct make_vector
{
    template <class ValueType> struct apply { typedef vector<ValueType> type; };
};

// scalar_vector_types: vector types with floating point values
typedef mpl::transform<scalar_types,make_vector>::type scalar_vector_types;

// integer_vector_types: vector types with floating point values
typedef mpl::transform<integer_types,make_vector>::type integer_vector_types;

// floating_vector_types: vector types with floating point values
typedef mpl::transform<floating_types,make_vector>::type floating_vector_types;


// metafunction to get the associated property map types to ValueTypes and
// IndexMap

struct property_map_type
{
    template <class ValueType, class IndexMap>
    struct apply
    {
        typedef vector_property_map<ValueType,IndexMap> type;
    };                           
};

struct property_map_types
{
   // this wraps an unary metafunction class Bind into a unary metafunction,
   // i.e., it is an identity operation. I'm not sure why it's necessary, but
   // using pure unary bind expressions below didn't work for me, and this fixed
   // it.
    template <class Bind>
    struct bind_wrap1
    {
        template <class T1> struct apply 
        { typedef typename Bind::template apply<T1>::type type; };
    };

    template <class ValueTypes, class IndexMap, 
              class IncludeIndexMap = mpl::bool_<false> >
    struct apply
    {
        typedef typename mpl::transform<
            ValueTypes,
            bind_wrap1<mpl::bind2<property_map_type, 
                                  mpl::_1, 
                                  IndexMap> >
            >::type scalar_properties;

        // put index map itself        
        typedef typename mpl::if_<
            IncludeIndexMap,
            typename mpl::push_back<scalar_properties,IndexMap>::type,
            scalar_properties
            >::type type;
    };                           
};


// this function gets the "static" property map behind the dynamic property map,
// or throws bad_cast if it doesn't match the correct type
template <class PropertyMap>
PropertyMap& get_static_property_map(dynamic_property_map& map)
{
    return dynamic_cast
        <detail::dynamic_property_map_adaptor<PropertyMap>&>(map).base();
}

// same as above, but returns a pointer and does not throw an exception, and
// returns 0 in case of failure
template <class PropertyMap>
PropertyMap* get_static_property_map(dynamic_property_map* map)
{
    detail::dynamic_property_map_adaptor<PropertyMap>* adaptor =
        dynamic_cast
        <detail::dynamic_property_map_adaptor<PropertyMap>*>(map);
    if (adaptor)
        return &adaptor->base();
    else
        return 0;
}

// this function gets the dynamic property map inside dp which matches the given
// name and key type
dynamic_property_map&
find_property_map(const dynamic_properties& dp, string name,
                  const type_info& key_type);

// convenience function which finds and returns the appropriate static property
// from the dynamic properties
template <class PropertyMap>
PropertyMap& find_static_property_map(const dynamic_properties& dp, 
                                      string name)
{
    typedef typename property_traits<PropertyMap>::key_type key_type;
    dynamic_property_map& dmap = find_property_map(dp, name, 
                                                   typeid(key_type));
    return get_static_property_map<PropertyMap>(dmap);
}

// this class contains a copy of a dynamic_properties, which does not delete its
// members when it deconstructs
struct dynamic_properties_copy: public dynamic_properties
{
    dynamic_properties_copy() {}
    dynamic_properties_copy(dynamic_properties& dp)
        : dynamic_properties(dp) {}
    dynamic_properties_copy
        (const function<auto_ptr<dynamic_property_map>
         (const string&, const boost::any&, const boost::any&)>& fn)
            : dynamic_properties(fn) {}
    ~dynamic_properties_copy()
{
        for (typeof(this->begin()) iter = this->begin(); iter != this->end();
             ++iter)
            iter->second = 0; // will be deleted when original dp deconstructs
    }
};

// this function gets the value in the map corresponding to the given key,
// converted to ConvertedType, or throws bad_lexical_cast.

template <class T>
struct AttemptAnyConversion; // forward declaration

template <class ConvertedType, class Key>
ConvertedType get_converted_scalar_value(dynamic_property_map& dmap,
                                         const Key& key)
{
    ConvertedType target = ConvertedType();
    const boost::any& source = dmap.get(key);
    bool success;
    if (dmap.value() == typeid(ConvertedType))
    {
        target = any_cast<ConvertedType>(dmap.get(key));
        success = true;
    }
    else
    {
        mpl::for_each<scalar_types>
            (AttemptAnyConversion<ConvertedType>(target, source, success));
    }
    if (!success)
        throw bad_lexical_cast();
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
            _value = lexical_cast<T>(any_cast<Source>(_source));
            _success = true;
        }
        catch (bad_any_cast){}
        catch (bad_lexical_cast){}
    }

    T& _value;
    const boost::any& _source;
    bool& _success;
};

// the following class wraps a dynamic_property_map, so it can be used as a
// regular property map. The values are converted to the desired Value type,
// which may cause a performance impact. Should be used only when property map
// access time is not crucial
template <class Value, class Key>
class DynamicPropertyMapWrap
{
public:
    typedef Value value_type;
    typedef Value reference;
    typedef Key key_type;
    typedef read_write_property_map_tag category;

    DynamicPropertyMapWrap(dynamic_property_map& dmap):_dmap(&dmap) {}
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
    dynamic_property_map* _dmap;
};

template <class Value, class Key>
Value get(const graph_tool::DynamicPropertyMapWrap<Value,Key>& pmap,
          const Key& k)
{
    return pmap.get(k);
}

template <class Value, class Key>
void put(graph_tool::DynamicPropertyMapWrap<Value,Key> pmap,
         const Key& k, const Value& val)
{
    pmap.put(k, val);
}

// the following is hash functor which, will hash a vertex or edge descriptor
// based on its index
template <class IndexMap>
class DescriptorHash
    : public unary_function<typename IndexMap::key_type, size_t>
{
public:
    DescriptorHash() {}
    DescriptorHash(IndexMap index_map): _index_map(index_map) {}
    size_t operator()(typename IndexMap::key_type const& d) const
    {
        return hash_value(_index_map[d]);
    }
private:
    IndexMap _index_map;
};

// the following is a property map based on a hashed container, which uses the
// above hash function for vertex or edge descriptors
template <class IndexMap, class Value>
class HashedDescriptorMap
    : public put_get_helper<Value&, HashedDescriptorMap<IndexMap,Value> >
{
public:
    typedef DescriptorHash<IndexMap> hashfc_t;
    typedef tr1::unordered_map<typename IndexMap::key_type,Value,hashfc_t>
        map_t;
    typedef associative_property_map<map_t> prop_map_t;

    typedef typename property_traits<prop_map_t>::value_type value_type;
    typedef typename property_traits<prop_map_t>::reference reference;
    typedef typename property_traits<prop_map_t>::key_type key_type;
    typedef typename property_traits<prop_map_t>::category category;

    HashedDescriptorMap(IndexMap index_map)
        : _base_map(new map_t(0, hashfc_t(index_map))), _prop_map(*_base_map) {}
    HashedDescriptorMap(){}

    reference operator[](const key_type& k) { return _prop_map[k]; }
    const reference operator[](const key_type& k) const { return _prop_map[k]; }

private:
    shared_ptr<map_t> _base_map;
    prop_map_t _prop_map;
};


// this wraps a container as a property map which is automatically initialized
// with a given default value
template <class Container>
class InitializedPropertyMap
    : public put_get_helper<typename Container::value_type::second_type&,
                            InitializedPropertyMap<Container> >
{
public:
    typedef typename Container::value_type::second_type value_type;
    typedef value_type& reference;
    typedef typename Container::key_type key_type;
    typedef read_write_property_map_tag category;

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

// the following is a property map which always returns a constant value
template <class Value, class Key>
class ConstantPropertyMap
    : public put_get_helper<Value, ConstantPropertyMap<Value,Key> >
{
public:
    typedef Value value_type;
    typedef value_type& reference;
    typedef Key key_type;
    typedef readable_property_map_tag category;

    ConstantPropertyMap(value_type c): _c(c) {}
    ConstantPropertyMap(){}

    const value_type& operator[](const key_type& k) const
    {
        return _c;
    }

private:
    value_type _c;
};

// this wraps an existing property map, but always converts its values to a
// given type
template <class PropertyMap, class Type>
class ConvertedPropertyMap
{
public:
    typedef Type value_type;
    typedef typename property_traits<PropertyMap>::value_type orig_type;
    typedef value_type reference;
    typedef typename property_traits<PropertyMap>::key_type key_type;
    typedef read_write_property_map_tag category;

    ConvertedPropertyMap(PropertyMap base_map)
        : _base_map(base_map) {}
    ConvertedPropertyMap(){}

    value_type get(const key_type& k) const
    {
        return boost::get(_base_map, k);
    }

    void put(const key_type& k, const Type& v)
    {
        put(_base_map, k, orig_type(v));
    }
private:
    PropertyMap _base_map;
};

template <class PropertyMap, class Type>
Type get(const graph_tool::ConvertedPropertyMap<PropertyMap,Type>& pmap,
         const typename property_traits<PropertyMap>::key_type& k)
{
    return pmap.get(k);
}

template <class PropertyMap, class Type>
void put(graph_tool::ConvertedPropertyMap<PropertyMap,Type> pmap,
         const typename property_traits<PropertyMap>::key_type& k,
         const typename property_traits<PropertyMap>::value_type& val)
{
    pmap.put(k, val);
}

struct get_static_prop
{
    get_static_prop(dynamic_property_map& dmap, boost::any& smap)
        : _dmap(dmap), _smap(smap) {}
    
    template <class PropertyMap>
    void operator()(PropertyMap) const
    {
        PropertyMap* smap = get_static_property_map<PropertyMap>(&_dmap);
        if (smap != 0)
            _smap = boost::any(*smap);
    }

    dynamic_property_map& _dmap;
    boost::any& _smap;
};

template <class IndexMap>
boost::any prop(const string& name, IndexMap, 
                const dynamic_properties& dp)
{
    typedef typename property_traits<IndexMap>::key_type key_t;
    dynamic_property_map& dmap = 
        find_property_map(dp, name, typeid(key_t));
    boost::any prop;
    typedef typename property_map_types::apply<value_types,IndexMap>::type 
        properties_t;
    mpl::for_each<properties_t>(get_static_prop(dmap, prop));
    return prop;
}

// this functor tests whether or not a given boost::any object holds a type
// contained in a given type Sequence
template <class Sequence>
struct belongs
{
    struct get_type
    {
        get_type(const boost::any& val, bool& found)
            : _val(val), _found(found) {}
        
        template <class Type>
        void operator()(Type) const
        {
            const Type* ptr = any_cast<Type>(&_val);
            if (ptr != 0)
                _found = true;
        }

        const boost::any& _val;
        bool& _found;
    };

    bool operator()(const boost::any& prop)
    {
        bool found = false;
        mpl::for_each<Sequence>(get_type(prop, found));
        return found;
    }
};

// this will return the name of a given type
template <class TypeSequence = value_types>
class get_type_name
{
public:
    get_type_name(const char* type_names[])
        : _type_names(type_names) {}

    get_type_name()
        : _type_names(type_names) {}

    string operator()(const type_info& type) const
    {
        using namespace lambda;
        string name;
        mpl::for_each<TypeSequence>
            (lambda::bind<void>(find_name(), lambda::_1, var(type), 
                                var(_type_names), var(name)));
        return name;
    }

private:
    struct find_name
    {
        template <class Type>
        void operator()(Type, const type_info& type, 
                        const char** type_names,
                        string& name) const
        {
            size_t index = mpl::find<TypeSequence,Type>::type::pos::value;
            if (type == typeid(Type))
                name = type_names[index];
        }
    };

    const char** _type_names;
};

} // graph_tool namespace

#endif
