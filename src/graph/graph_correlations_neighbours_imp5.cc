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

#include "graph_correlations_neighbours.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

typedef mpl::vector<long double> weight_scalar_types;

struct weight_map_types: property_map_types::apply<
    weight_scalar_types,
    GraphInterface::edge_index_map_t,
    mpl::bool_<false> >::type {};

void graph_correlations_neighbours_imp5(const GraphInterface& g, 
                                        avg_corr_t& avg_corr,
                                        boost::any origin_deg,
                                        boost::any neighbours_deg,
                                        boost::any weight)
{
    run_action<>()(g, get_average_nearest_neighbours_correlation<avg_corr_t>
                   (avg_corr), all_selectors(), all_selectors(), 
                   weight_map_types())(origin_deg, neighbours_deg, weight);    
}
