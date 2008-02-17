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
#include "shared_map.hh"

#include <boost/lambda/bind.hpp>

#include "graph_correlations_combined.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

void graph_correlations_combined_imp1(const GraphInterface& g,
                                      hist2d_t& hist,
                                      boost::any deg1, boost::any deg2)
{
    run_action<>()(g, bind<void>(get_combined_degree_histogram(),
                                 _1, _2, _3, var(hist)),
                   all_selectors(),
                   graph_tool::detail::split::apply<all_selectors>::type
                   ::second())
        (deg1, deg2);
}
