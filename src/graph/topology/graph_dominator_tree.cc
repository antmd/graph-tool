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

#include <boost/graph/dominator_tree.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

struct get_dominator_tree
{
    template <class Graph, class PredMap>
    void operator()(const Graph& g, size_t entry, PredMap pred_map) const
    {
        lengauer_tarjan_dominator_tree(g, vertex(entry,g), pred_map);
    }
};

typedef property_map_types::apply<mpl::vector<int32_t>,
                                  GraphInterface::vertex_index_map_t,
                                  mpl::bool_<false> >::type
    pred_properties;

void dominator_tree(GraphInterface& gi, size_t entry, boost::any pred_map)
{
    run_action<graph_tool::detail::always_directed>()
        (gi, std::bind(get_dominator_tree(), placeholders::_1, entry,
                       placeholders::_2),
         pred_properties())(pred_map);
}
