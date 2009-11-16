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

#include "graph_filtering.hh"

#include <boost/python.hpp>

#include "graph.hh"
#include "graph_selectors.hh"

#include "graph_absolute_trust.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

void absolute_trust(GraphInterface& g, int64_t source, boost::any c,
                    boost::any t, size_t n_paths, size_t n_paths_vertex,
                    double epsilon, bool reversed)
{
    if (!belongs<edge_floating_properties>()(c))
        throw ValueException("edge property must be of floating point value type");
    if (!belongs<vertex_floating_vector_properties>()(t))
        throw ValueException("vertex property must be of floating point vector value type");

    run_action<>()(g,
                   bind<void>(get_absolute_trust(), _1, g.GetVertexIndex(),
                              g.GetEdgeIndex(), g.GetMaxEdgeIndex(), source, _2,
                              _3, make_tuple(n_paths, n_paths_vertex, epsilon),
                              reversed), edge_floating_properties(),
                   vertex_floating_vector_properties())(c, t);
}

void export_absolute_trust()
{
    using namespace boost::python;
    def("get_absolute_trust", &absolute_trust);
}
