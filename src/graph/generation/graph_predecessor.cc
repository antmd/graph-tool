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

#include "graph_predecessor.hh"

using namespace graph_tool;
using namespace boost;

void predecessor_graph(GraphInterface& gi, GraphInterface& gpi,
                       boost::any pred_map)
{
    run_action<>()(gi, std::bind(get_predecessor_graph(), placeholders::_1,
                                 gi.GetVertexIndex(), std::ref(gpi.GetGraph()),
                                 placeholders::_2),
                   vertex_scalar_properties())(pred_map);
}
