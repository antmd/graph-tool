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

#include <boost/mpl/joint_view.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

// this will return the vertex average of degrees or scalar properties
python::object
get_vertex_average(GraphInterface& gi, GraphInterface::deg_t deg)
{
    python::object a, dev;
    size_t count;

    typedef mpl::joint_view<vertex_scalar_properties,
                            vertex_scalar_vector_properties>
        vertex_numeric_properties_c;

    typedef property_map_types::apply<mpl::vector<python::object>,
                                      GraphInterface::vertex_index_map_t>::type
        python_properties;

    typedef mpl::joint_view<vertex_numeric_properties_c,
                            python_properties>
        vertex_numeric_properties;

    typedef mpl::transform<vertex_numeric_properties,
                           scalar_selector_type,
                           mpl::back_inserter<degree_selectors> >::type
        numeric_selectors;

    run_action<>()(gi, get_average<VertexAverageTraverse>(a,dev,count),
                   numeric_selectors())(degree_selector(deg));
    return python::make_tuple(a, dev, count);
}

// this will return the edge average of scalar properties
python::object
get_edge_average(GraphInterface& gi, boost::any prop)
{
    typedef mpl::joint_view<edge_scalar_properties,
                            edge_scalar_vector_properties>
        edge_numeric_properties_c;

    typedef property_map_types::apply<mpl::vector<python::object>,
                                      GraphInterface::edge_index_map_t>::type
        python_properties;

    typedef mpl::joint_view<edge_numeric_properties_c,
                            python_properties>
        edge_numeric_properties;

    python::object a, dev;
    size_t count;
    run_action<graph_tool::detail::always_directed>()
        (gi, get_average<EdgeAverageTraverse>(a, dev, count),
         edge_numeric_properties())(prop);
    return python::make_tuple(a, dev, count);
}

using namespace boost::python;

void export_average()
{
    def("get_vertex_average", &get_vertex_average);
    def("get_edge_average", &get_edge_average);
}
