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

#include "graph_parallel.hh"

#include <boost/python.hpp>

#include <iostream>

using namespace std;
using namespace boost;
using namespace graph_tool;

void do_label_parallel_edges(GraphInterface& gi, boost::any property,
                             bool mark_only, bool count_all)
{
    run_action<>()(gi, std::bind(label_parallel_edges(), placeholders::_1,
                                 placeholders::_2, mark_only, count_all),
                   writable_edge_scalar_properties())(property);
}

void do_label_self_loops(GraphInterface& gi, boost::any property,
                         bool mark_only)
{
    run_action<>()(gi, std::bind(label_self_loops(),
                                 placeholders::_1, placeholders::_2, mark_only),
                   writable_edge_scalar_properties())(property);
}

void do_remove_labeled_edges(GraphInterface& gi, boost::any property)
{
    run_action<graph_tool::detail::always_directed_never_reversed, mpl::true_>()
        (gi, std::bind(remove_labeled_edges(), placeholders::_1,
                       placeholders::_2),
         edge_scalar_properties())(property);
}

using namespace boost::python;

void export_parallel()
{
    def("label_parallel_edges", &do_label_parallel_edges);
    def("label_self_loops", &do_label_self_loops);
    def("remove_labeled_edges", &do_remove_labeled_edges);
}
