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

#include <boost/lambda/bind.hpp>

#include "graph_edge_correlations.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

hist3d_t
GraphInterface::GetEdgeVertexCorrelationHistogram(GraphInterface::deg_t deg1,
                                                  string edge_scalar,
                                                  GraphInterface::deg_t deg2)
    const
{
    hist3d_t hist;

    try
    {
        run_action<>()(*this, get_edge_correlation_histogram<hist3d_t>(hist),
                       all_selectors(), mpl::vector<edge_index_map_t>(),
                       all_selectors())
            (degree_selector(deg1, _properties),
             prop(edge_scalar, _edge_index, _properties),
             degree_selector(deg2, _properties));
    }
    catch (dynamic_get_failure& e)
    {
        throw GraphException("error getting scalar property: " +
                             string(e.what()));
    }
    return hist;
}
