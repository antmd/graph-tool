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

#ifndef GRAPH_INCIDENCE_HH
#define GRAPH_INCIDENCE_HH

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_util.hh"

namespace graph_tool
{
using namespace boost;

struct get_incidence
{
    template <class Graph, class VIndex, class EIndex>
    void operator()(Graph& g, VIndex vindex, EIndex eindex,
                    multi_array_ref<double,1>& data,
                    multi_array_ref<int32_t,1>& i,
                    multi_array_ref<int32_t,1>& j) const
    {
        int pos = 0;
        for (auto v : vertices_range(g))
        {
            for (const auto& e : out_edges_range(v, g))
            {
                if (is_directed::apply<Graph>::type::value)
                    data[pos] = -1;
                else
                    data[pos] = 1;
                i[pos] = get(vindex, v);
                j[pos] = get(eindex, e);
                ++pos;
            }

            for (const auto& e : in_edges_range(v, g))
            {
                data[pos] = 1;
                i[pos] = get(vindex, v);
                j[pos] = get(eindex, e);
                ++pos;
            }
        }
    }
};

} // namespace graph_tool

#endif // GRAPH_INCIDENCE_HH
