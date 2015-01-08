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

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_properties.hh"

#include "graph_average.hh"

#include <boost/python.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

// this will return the vertex average of degrees or scalar properties
python::object
get_vertex_average(GraphInterface& gi, GraphInterface::deg_t deg)
{
    long double a, dev;
    run_action<>()(gi, get_average<VertexAverageTraverse>(a,dev),
                   scalar_selectors())(degree_selector(deg));
    return python::make_tuple(a,dev);
}

// this will return the edge average of scalar properties
python::object
get_edge_average(GraphInterface& gi, boost::any prop)
{
    long double a, dev;
    bool directed = gi.GetDirected();
    gi.SetDirected(true);
    run_action<graph_tool::detail::always_directed>()
        (gi, get_average<EdgeAverageTraverse>(a, dev),
         edge_scalar_properties())(prop);
    gi.SetDirected(directed);

    return python::make_tuple(a, dev);
}

using namespace boost::python;

void export_average()
{
    def("get_vertex_average", &get_vertex_average);
    def("get_edge_average", &get_edge_average);
}
