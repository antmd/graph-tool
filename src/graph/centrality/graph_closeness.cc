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

#include "graph_filtering.hh"
#include "graph.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"

#include "graph_closeness.hh"

#include <functional>

using namespace std;
using namespace graph_tool;

void do_get_closeness(GraphInterface& gi, boost::any weight,
                      boost::any closeness, bool harmonic, bool norm)
{
    if (weight.empty())
    {
        run_action<>()(gi,
                       std::bind(get_closeness(), placeholders::_1,
                                 gi.GetVertexIndex(), no_weightS(),
                                 placeholders::_2, harmonic, norm),
                       writable_vertex_scalar_properties())(closeness);
    }
    else
    {
        run_action<>()(gi,
                       std::bind(get_closeness(), placeholders::_1,
                                 gi.GetVertexIndex(), placeholders::_2,
                                 placeholders::_3, harmonic, norm),
                       edge_scalar_properties(),
                       writable_vertex_scalar_properties())(weight, closeness);
    }
}

void export_closeness()
{
    boost::python::def("closeness", &do_get_closeness);
}
