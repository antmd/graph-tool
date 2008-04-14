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

#include "graph_parallel.hh"

#include <boost/python.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

void do_label_parallel_edges(GraphInterface& gi, const string& property)
{
    run_action<>()(gi, bind<void>(label_parallel_edges(), _1,
                                  gi.GetEdgeIndex(), _2),
                   edge_scalar_properties())
            (edge_prop(property, gi));
}

void do_label_self_loops(GraphInterface& gi, const string& property)
{
    run_action<>()(gi, bind<void>(label_self_loops(), _1,
                                  gi.GetEdgeIndex(), _2),
                   edge_scalar_properties())
            (edge_prop(property, gi));
}


using namespace boost::python;

void export_parallel()
{
    def("label_parallel_edges", &do_label_parallel_edges);
    def("label_self_loops", &do_label_self_loops);
}
