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
#include "graph_selectors.hh"

#include "graph_kcore.hh"

#include <boost/python.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

void do_kcore_decomposition(GraphInterface& gi, boost::any prop,
                            GraphInterface::deg_t deg)
{
    run_action<>()(gi, std::bind(kcore_decomposition(), placeholders::_1,
                                 gi.GetVertexIndex(), placeholders::_2,
                                 placeholders::_3),
                   writable_vertex_scalar_properties(),
                   degree_selectors())(prop, degree_selector(deg));
}

void export_kcore()
{
    python::def("kcore_decomposition", &do_kcore_decomposition);
};
