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

#include <boost/lambda/bind.hpp>

#include "graph.hh"
#include "graph_filtering.hh"

#include <tr1/random>
typedef std::tr1::mt19937 rng_t;

#include "graph_rewiring.hh"

using namespace graph_tool;
using namespace boost;
using namespace boost::lambda;

void random_rewire(GraphInterface& gi, string strat, bool self_loops,
                   bool parallel_edges, size_t seed)
{
    rng_t rng(static_cast<rng_t::result_type>(seed));

    GraphInterface::edge_index_map_t edge_index =
        any_cast<GraphInterface::edge_index_map_t>(gi.GetEdgeIndex());
    if (strat == "uncorrelated")
        run_action<graph_tool::detail::never_reversed>()
            (gi, bind<void>(graph_rewire<RandomRewireStrategy>(),
                            _1, edge_index, ref(rng), self_loops,
                            parallel_edges))();
    else if (strat == "correlated")
        run_action<graph_tool::detail::never_reversed>()
            (gi, bind<void>(graph_rewire<CorrelatedRewireStrategy>(),
                            _1, edge_index, ref(rng), self_loops,
                            parallel_edges))();
    else
        throw ValueException("invalid random rewire strategy: " + strat);
}
