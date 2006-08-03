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

#ifndef GRAPH_PROPERTIES_HH
#define GRAPH_PROPERTIES_HH

#include <string>
#include <boost/dynamic_property_map.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/for_each.hpp>

namespace graph_tool 
{

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


} // namespace graph_tool


namespace boost {

using namespace graph_tool;

template <class Value, class Key>
Value get(const DynamicPropertyMapWrap<Value,Key>& pmap, typename property_traits<DynamicPropertyMapWrap<Value,Key> >::key_type k)
{
    return pmap.get(k);
}

template <class Value, class Key>
void put(DynamicPropertyMapWrap<Value,Key> pmap, typename property_traits<DynamicPropertyMapWrap<Value,Key> >::key_type k, 
	 typename property_traits<DynamicPropertyMapWrap<Value,Key> >::value_type val)
{
    pmap.put(k,val);
}

} // namespace boost

#endif
