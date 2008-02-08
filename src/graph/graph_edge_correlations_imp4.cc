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

#include "graph.hh"
#include "histogram.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"

#include "graph_edge_correlations.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

typedef mpl::vector<double> edge_scalar_types;

struct edge_map_types: property_map_types::apply<
    edge_scalar_types,
    GraphInterface::edge_index_map_t,
    mpl::bool_<false> >::type {};

void graph_edge_correlations_imp4(const GraphInterface& g, 
                                  hist3d_t& hist,
                                  boost::any deg1,
                                  boost::any edge_prop,
                                  boost::any deg2)
{
    run_action<>()(g, get_edge_correlation_histogram<hist3d_t>(hist),
                   all_selectors(), edge_map_types(), all_selectors())
        (deg1, edge_prop, deg2);    
}
