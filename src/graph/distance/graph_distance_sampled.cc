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

#include "graph_distance_sampled.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

hist_t
GraphInterface::GetSampledDistanceHistogram(string weight, size_t samples,
                                            size_t seed) const
{
    hist_t hist;

    if (weight == "")
    {
        run_action<>()(*this, bind<void>(get_sampled_distances(), _1,
                                         _vertex_index, no_weightS(), var(hist),
                                         samples, seed))();
    }
    else
    {
        try
        {
            run_action<>()(*this,
                           bind<void>(get_sampled_distances(), _1,
                                      _vertex_index, _2, var(hist), samples,
                                      seed), edge_scalar_properties())
                (prop(weight, _edge_index, _properties));
        }
        catch (property_not_found& e)
        {
            throw GraphException("error getting scalar property: " +
                                 string(e.what()));
        }
    }
    return hist;
}
