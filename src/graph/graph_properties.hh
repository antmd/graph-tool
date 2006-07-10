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

namespace graph_tool 
{

//==============================================================================
// Property Map Utility Functions
//==============================================================================

template <class PropertyMap>
PropertyMap& get_static_property_map(boost::dynamic_property_map& map)
{
    return dynamic_cast<boost::detail::dynamic_property_map_adaptor<PropertyMap>&>(map).base();
}

boost::dynamic_property_map& find_property_map(boost::dynamic_properties& dp, std::string name, const std::type_info& key_type);

struct dynamic_properties_copy: public boost::dynamic_properties
{
    dynamic_properties_copy() {}
    dynamic_properties_copy(boost::dynamic_properties& dp): boost::dynamic_properties(dp) {}
    dynamic_properties_copy(const boost::function<std::auto_ptr<boost::dynamic_property_map>(const std::string&, const boost::any&, const boost::any&)>& fn)
	: boost::dynamic_properties(fn) {}
    ~dynamic_properties_copy()
    {	
	for (typeof(this->begin()) iter = this->begin(); iter != this->end(); ++iter)
	    iter->second = 0; // will be deleted when dp deconstructs
    }
};

} // namespace graph_tool

#endif
