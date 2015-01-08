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
#include "histogram.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"

#include "graph_corr_hist.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;


void graph_correlations_imp1(GraphInterface& g, boost::python::object& hist,
                             boost::python::object& ret_bins,
                             boost::any deg1, boost::any deg2,
                             boost::any weight,
                             const std::array<vector<long double>,2>& bins)
{
    typedef DynamicPropertyMapWrap<long double, GraphInterface::edge_t>
        wrapped_weight_t;
    run_action<>()(g, get_correlation_histogram<GetNeighboursPairs>
                   (hist, bins, ret_bins),
                   scalar_selectors(), scalar_selectors(),
                   boost::mpl::vector<wrapped_weight_t>())
        (deg1, deg2, weight);
}
