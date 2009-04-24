// Copyright (C) 2008  Tiago de Paula Peixoto <tiago@forked.de>
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

#include "graph_search.hh"

#include <boost/python.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

// find vertices which match a certain (inclusive) property range
python::list
find_vertex_range(GraphInterface& gi, GraphInterface::deg_t deg,
                  python::tuple range)
{
    python::list ret;

    run_action<>()(gi, lambda::bind<void>(find_vertices(), lambda::_1,
                                          lambda::var(gi), lambda::_2,
                                          lambda::var(range), lambda::var(ret)),
                   all_selectors())(degree_selector(deg));
    return ret;
}

// find edges which match a certain (inclusive) property range
python::list
find_edge_range(GraphInterface& gi, boost::any eprop,
                python::tuple range)
{
    python::list ret;
    typedef property_map_types::apply<value_types,
                                      GraphInterface::edge_index_map_t,
                                      mpl::bool_<true> >::type
        all_edge_props;

    GraphInterface::edge_index_map_t eindex =
        any_cast<GraphInterface::edge_index_map_t>(gi.GetEdgeIndex());
    run_action<>()(gi, lambda::bind<void>(find_edges(), lambda::_1,
                                          lambda::var(gi), eindex,
                                          lambda::_2, lambda::var(range),
                                          lambda::var(ret)),
                   all_edge_props())(eprop);
    return ret;
}


using namespace boost::python;

void export_search()
{
    def("find_vertex_range", &find_vertex_range);
    def("find_edge_range", &find_edge_range);
}
