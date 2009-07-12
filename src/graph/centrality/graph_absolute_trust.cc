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

#include <tr1/random>
typedef std::tr1::mt19937 rng_t;

#include "graph_absolute_trust.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

void absolute_trust(GraphInterface& g, int64_t source, boost::any c,
                    boost::any t, double epslon, size_t min_iter,
                    size_t max_iter, bool reversed, size_t seed)
{
    rng_t rng(static_cast<rng_t::result_type>(seed));

    if (!belongs<edge_floating_properties>()(c))
        throw GraphException("edge property must be of floating point value type");
    if (!belongs<vertex_floating_vector_properties>()(t))
        throw GraphException("vertex property must be of floating point vector value type");

    run_action<>()(g,
                   bind<void>(get_absolute_trust(), _1, g.GetVertexIndex(),
                              source, _2, _3, epslon,
                              make_pair(min_iter, max_iter), reversed,
                              ref(rng)),
                   edge_floating_properties(),
                   vertex_floating_vector_properties())(c, t);
}

void export_absolute_trust()
{
    using namespace boost::python;
    def("get_absolute_trust", &absolute_trust);
}
