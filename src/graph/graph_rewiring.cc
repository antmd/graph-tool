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
#include "graph_rewiring.hh"

using namespace graph_tool;
using namespace boost;
using namespace boost::lambda;

void GraphInterface::RandomRewire(string strat, bool self_loops,
                                  bool parallel_edges, size_t seed)
{
    bool reversed = GetReversed();
    SetReversed(false);

    if (strat == "uncorrelated")
        run_action<detail::never_reversed>()
            (*this, bind<void>(graph_rewire<RandomRewireStrategy>(),
                               _1, _edge_index, seed, self_loops,
                               parallel_edges))();
    else if (strat == "correlated")
        run_action<detail::never_reversed>()
            (*this, bind<void>(graph_rewire<CorrelatedRewireStrategy>(),
                               _1, _edge_index, seed, self_loops,
                               parallel_edges))();
    else
        throw GraphException("invalid random rewire stategy: " + strat);
    SetReversed(reversed);
}
