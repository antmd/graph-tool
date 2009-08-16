// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007  Tiago de Paula Peixoto <tiago@forked.de>
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

#include <boost/python.hpp>
#include "graph.hh"

using namespace boost;
using namespace boost::python;
using namespace graph_tool;

bool check_isomorphism(GraphInterface& gi1, GraphInterface& gi2,
                       boost::any iso_map);

bool get_kruskal_spanning_tree(GraphInterface& gi, boost::any weight_map,
                               boost::any tree_map);

bool get_prim_spanning_tree(GraphInterface& gi, size_t root,
                            boost::any weight_map, boost::any tree_map);

void topological_sort(GraphInterface& gi, vector<int32_t>& sort);

bool denominator_tree(GraphInterface& gi, size_t entry, boost::any pred_map);

void transitive_closure(GraphInterface& gi, GraphInterface& tcgi);

void export_components();

BOOST_PYTHON_MODULE(libgraph_tool_topology)
{
    def("check_isomorphism", &check_isomorphism);
    def("get_kruskal_spanning_tree", &get_kruskal_spanning_tree);
    def("get_prim_spanning_tree", &get_prim_spanning_tree);
    def("topological_sort", &topological_sort);
    def("denominator_tree", &denominator_tree);
    def("transitive_closure", &transitive_closure);
    export_components();
}
