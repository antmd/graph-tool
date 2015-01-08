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

#include <boost/python.hpp>

#include "graph.hh"
#include "graph_selectors.hh"

#include "graph_trust_transitivity.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

void trust_transitivity(GraphInterface& g, int64_t source, int64_t target,
                    boost::any c, boost::any t)
{
    if (!belongs<edge_floating_properties>()(c))
        throw ValueException("edge property must be of floating point value type");
    if (!belongs<vertex_floating_vector_properties>()(t))
        throw ValueException("vertex property must be of floating point valued vector type");

    run_action<>()(g,
                   bind<void>(get_trust_transitivity(), _1, g.GetVertexIndex(),
                              source, target, _2, _3),
                   edge_floating_properties(),
                   vertex_floating_vector_properties())(c,t);
}

void export_trust_transitivity()
{
    using namespace boost::python;
    def("get_trust_transitivity", &trust_transitivity);
}
