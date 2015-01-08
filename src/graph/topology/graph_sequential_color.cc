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

#include "graph_filtering.hh"
#include "graph.hh"
#include "graph_properties.hh"

#include <boost/graph/sequential_vertex_coloring.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

struct get_coloring
{
    template <class Graph, class OrderMap, class ColorMap>
    void operator()(Graph& g, OrderMap order, ColorMap color, size_t& nc) const
    {
        nc = sequential_vertex_coloring(g, order, color);
    }
};

typedef property_map_types::apply<mpl::vector<int32_t, int64_t>,
                                  GraphInterface::vertex_index_map_t,
                                  mpl::bool_<false> >::type
    int_properties;

size_t sequential_coloring(GraphInterface& gi, boost::any order,
                           boost::any color)
{
    size_t nc = 0;
    run_action<>()
        (gi, std::bind(get_coloring(), placeholders::_1, placeholders::_2,
                       placeholders::_3, std::ref(nc)),
         vertex_integer_properties(), int_properties())(order, color);
    return nc;
}
