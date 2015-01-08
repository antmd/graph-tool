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
#include "graph_properties.hh"

#include <boost/graph/metric_tsp_approx.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

struct get_tsp_approx
{
    template <class Graph, class WeightMap, class IntType>
    void operator()(Graph& g, size_t src, WeightMap weights,
                    vector<IntType>& tour) const
    {
        back_insert_iterator<vector<IntType> > back_it(tour);
        metric_tsp_approx_tour_from_vertex(g, vertex(src, g), weights,
                                           back_it);
    }
};

vector<int32_t> get_tsp(GraphInterface& gi, size_t src, boost::any weight_map)
{
    vector<int32_t> tour;

    typedef ConstantPropertyMap<size_t,GraphInterface::edge_t> cweight_t;

    if (weight_map.empty())
        weight_map = cweight_t(1);

    typedef mpl::push_back<edge_scalar_properties, cweight_t>::type
        weight_maps;

    run_action<graph_tool::detail::never_directed>()
        (gi, std::bind(get_tsp_approx(), placeholders::_1, src,
                       placeholders::_2, std::ref(tour)),
         weight_maps())(weight_map);

    return tour;
}
