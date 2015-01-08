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

#include <boost/python.hpp>
#include "graph.hh"
#include "random.hh"

using namespace boost;
using namespace boost::python;
using namespace graph_tool;

bool check_isomorphism(GraphInterface& gi1, GraphInterface& gi2,
                       boost::any ainv_map1, boost::any ainv_map2,
                       int64_t max_inv, boost::any aiso_map);
void get_kruskal_spanning_tree(GraphInterface& gi, boost::any weight_map,
                               boost::any tree_map);
void get_prim_spanning_tree(GraphInterface& gi, size_t root,
                            boost::any weight_map, boost::any tree_map);
bool topological_sort(GraphInterface& gi, vector<int32_t>& sort);
void dominator_tree(GraphInterface& gi, size_t entry, boost::any pred_map);
void transitive_closure(GraphInterface& gi, GraphInterface& tcgi);
bool is_planar(GraphInterface& gi, boost::any embed_map, boost::any kur_map);
void maximal_planar(GraphInterface& gi);
void subgraph_isomorphism(GraphInterface& gi1, GraphInterface& gi2,
                          boost::any vertex_label1, boost::any vertex_label2,
                          boost::any edge_label1, boost::any edge_label2,
                          python::list vmapping, size_t max_n, bool induced,
                          bool iso);
double reciprocity(GraphInterface& gi);
size_t sequential_coloring(GraphInterface& gi, boost::any order,
                           boost::any color);
bool is_bipartite(GraphInterface& gi, boost::any part_map);
void get_random_spanning_tree(GraphInterface& gi, size_t root,
                              boost::any weight_map, boost::any tree_map,
                              rng_t& rng);
vector<int32_t> get_tsp(GraphInterface& gi, size_t src, boost::any weight_map);

void export_components();
void export_kcore();
void export_similarity();
void export_dists();
void export_all_dists();
void export_diam();
void export_random_matching();
void export_maximal_vertex_set();


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
    def("maximal_planar", &maximal_planar);
    def("reciprocity", &reciprocity);
    def("sequential_coloring", &sequential_coloring);
    def("is_bipartite", &is_bipartite);
    def("random_spanning_tree", &get_random_spanning_tree);
    def("get_tsp", &get_tsp);
    export_components();
    export_kcore();
    export_similarity();
    export_dists();
    export_all_dists();
    export_diam();
    export_random_matching();
    export_maximal_vertex_set();
}
