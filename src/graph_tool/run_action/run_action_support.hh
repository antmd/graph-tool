// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2015 Tiago de Paula Peixoto <tiago@skewed.de>
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

#include <map>
#include <set>
#include <list>
#include <tr1/unordered_set>
#include <tr1/unordered_map>
#include <tr1/tuple>
#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_properties.hh"
#include "histogram.hh"
#include <boost/bind.hpp>
#include <boost/python.hpp>

using namespace boost;
using namespace std;
using namespace graph_tool;

namespace graph_tool
{

// metafunction to get the correct property map
template <class IndexMap>
struct prop_bind_t
{
    template <class Value>
    struct as
    {
        typedef typename mpl::if_<std::is_same<Value,bool>,
                                  uint8_t, Value>::type val_t;
        typedef typename property_map_type::apply<val_t,IndexMap>::type type;
    };
};

// utility template function to extract the correct property map
template <class PropertyMap>
PropertyMap get_prop(py::object& map)
{
    try
    {
        python::object pmap(python::handle<>
                            (python::borrowed((PyObject*)(map))));

        python::object opmap = pmap.attr("_PropertyMap__map").attr("get_map")();
        return any_cast<PropertyMap>(python::extract<boost::any>(opmap)());
    }
    catch (bad_any_cast&)
    {
        throw GraphException("could not convert property map as requested");
    }
}
}
