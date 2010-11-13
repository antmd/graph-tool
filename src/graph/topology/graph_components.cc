// Copyright (C) 2008  Tiago de Paula Peixoto <tiago@skewed.de>
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
#include "graph_selectors.hh"

#include "graph_components.hh"

#include <boost/python.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

void  do_label_components(GraphInterface& gi, boost::any prop)
{
    run_action<>()(gi, label_components(),
                   writable_vertex_scalar_properties())(prop);
}

size_t do_label_biconnected_components(GraphInterface& gi, boost::any comp,
                                       boost::any art)
{
    size_t nc;
    run_action<graph_tool::detail::never_directed>()
        (gi, bind<void>(label_biconnected_components(), _1, _2, _3,
                        ref(nc)),
         writable_edge_scalar_properties(),
         writable_vertex_scalar_properties())
        (comp, art);
    return nc;
}

void export_components()
{
    python::def("label_components", &do_label_components);
    python::def("label_biconnected_components",
                &do_label_biconnected_components);
};
