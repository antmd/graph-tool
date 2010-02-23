// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007-2010 Tiago de Paula Peixoto <tiago@forked.de>
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

#include "graph.hh"
#include "graph_properties.hh"
#include "graph_filtering.hh"
#include "graph_selectors.hh"
#include "graph_util.hh"

#include <boost/mpl/for_each.hpp>

#include <boost/python/extract.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

namespace graph_tool
{

// global property types' names
const char* type_names[] =
    {"bool", "int32_t", "int64_t", "double", "long double",
     "string", "vector<bool>","vector<int32_t>", "vector<int64_t>",
     "vector<double>", "vector<long double>", "vector<string>",
     "python::object"};


struct shift_vertex_property
{
    template <class PropertyMap>
    void operator()(PropertyMap, const GraphInterface::multigraph_t& g,
                    boost::any map, size_t vi, bool& found) const
    {
        try
        {
            PropertyMap pmap = any_cast<PropertyMap>(map);
            for (size_t i = vi; i < num_vertices(g)-1; ++i)
                pmap[vertex(i,g)] = pmap[vertex(i+1,g)];
            found = true;
        }
        catch (bad_any_cast&) {}
    }
};

// this function will shift all the properties when a vertex is to be deleted
void GraphInterface::ShiftVertexProperty(boost::any prop, size_t index) const
{
    bool found = false;
    mpl::for_each<writable_vertex_properties>
        (bind<void>(shift_vertex_property(), _1, ref(_mg),
                    prop, index, ref(found)));
    if (!found)
        throw GraphException("invalid writable property map");
}

} // graph_tool namespace
