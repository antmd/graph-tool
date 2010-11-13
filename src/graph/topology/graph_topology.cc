// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007-2010 Tiago de Paula Peixoto <tiago@skewed.de>
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
void get_kruskal_spanning_tree(GraphInterface& gi, boost::any weight_map,
                               boost::any tree_map);
void get_prim_spanning_tree(GraphInterface& gi, size_t root,
                            boost::any weight_map, boost::any tree_map);
void topological_sort(GraphInterface& gi, vector<int32_t>& sort);
void dominator_tree(GraphInterface& gi, size_t entry, boost::any pred_map);
void transitive_closure(GraphInterface& gi, GraphInterface& tcgi);
bool is_planar(GraphInterface& gi, boost::any embed_map, boost::any kur_map);
void subgraph_isomorphism(GraphInterface& gi1, GraphInterface& gi2,
                          boost::any vertex_label1, boost::any vertex_label2,
                          boost::any edge_label1, boost::any edge_label2,
                          python::list vmapping, python::list emapping,
                          size_t n_max, size_t seed);

void export_components();
void export_dists();
void export_all_dists();

BOOST_PYTHON_MODULE(libgraph_tool_topology)
{
    def("check_isomorphism", &check_isomorphism);
    def("subgraph_isomorphism", &subgraph_isomorphism);
    def("get_kruskal_spanning_tree", &get_kruskal_spanning_tree);
    def("get_prim_spanning_tree", &get_prim_spanning_tree);
    def("topological_sort", &topological_sort);
    def("dominator_tree", &dominator_tree);
    def("transitive_closure", &transitive_closure);
    def("is_planar", &is_planar);
    export_components();
    export_dists();
    export_all_dists();
}
