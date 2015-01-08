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

using namespace graph_tool;
using namespace boost;

void edmonds_karp_max_flow(GraphInterface& gi, size_t src, size_t sink,
                           boost::any capacity, boost::any res);
void push_relabel_max_flow(GraphInterface& gi, size_t src, size_t sink,
                           boost::any capacity, boost::any res);
void kolmogorov_max_flow(GraphInterface& gi, size_t src, size_t sink,
                         boost::any capacity, boost::any res);
bool max_cardinality_matching(GraphInterface& gi, boost::any match);
double min_cut(GraphInterface& gi, boost::any weight, boost::any part_map);
void get_residual_graph(GraphInterface& gi, boost::any capacity, boost::any res,
                        boost::any oaugment);

#include <boost/python.hpp>
using namespace boost::python;

BOOST_PYTHON_MODULE(libgraph_tool_flow)
{
    def("edmonds_karp_max_flow", &edmonds_karp_max_flow);
    def("push_relabel_max_flow", &push_relabel_max_flow);
    def("kolmogorov_max_flow", &kolmogorov_max_flow);
    def("max_cardinality_matching", &max_cardinality_matching);
    def("min_cut", &min_cut);
    def("residual_graph", &get_residual_graph);
}
