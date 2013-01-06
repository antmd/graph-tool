// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2013 Tiago de Paula Peixoto <tiago@skewed.de>
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

#include <boost/graph/stoer_wagner_min_cut.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

struct get_min_cut
{
    template <class Graph, class EdgeWeight, class PartMap>
    void operator()(Graph& g, EdgeWeight eweight, PartMap part_map, double& mc) const
    {
        try
        {
            mc = stoer_wagner_min_cut(g, eweight, parity_map(part_map));
        }
        catch (bad_graph&)
        {
            throw ValueException("Graph has less than 2 vertices.");
        }
    }

};

double min_cut(GraphInterface& gi, boost::any weight, boost::any part_map)
{
    double mc = 0;

    typedef ConstantPropertyMap<size_t,GraphInterface::edge_t> cweight_t;

    if (weight.empty())
        weight = cweight_t(1);

    typedef mpl::push_back<writable_edge_scalar_properties, cweight_t>::type
        weight_maps;

    run_action<graph_tool::detail::never_directed>()
        (gi, bind<void>(get_min_cut(), _1, _2, _3, ref(mc)),
         weight_maps(), writable_vertex_scalar_properties())(weight, part_map);
    return mc;
}
