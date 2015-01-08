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

#ifndef GRAPH_ADJACENCY_MATRIX_HH
#define GRAPH_ADJACENCY_MATRIX_HH

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_util.hh"

namespace graph_tool
{
using namespace boost;

struct get_adjacency
{
    template <class Graph, class Index, class Weight>
    void operator()(Graph& g, Index index, Weight weight,
                    multi_array_ref<double,1>& data,
                    multi_array_ref<int,1>& i,
                    multi_array_ref<int,1>& j) const
    {
        int pos = 0;
        for (const auto& e : edges_range(g))
        {
            data[pos] = get(weight, e);
            i[pos] = get(index, target(e, g));
            j[pos] = get(index, source(e, g));

            ++pos;
            if (!is_directed::apply<Graph>::type::value)
            {
                data[pos] = get(weight, e);
                i[pos] = get(index, source(e, g));
                j[pos] = get(index, target(e, g));
                ++pos;
            }
        }
    }
};

} // namespace graph_tool

#endif // GRAPH_ADJACENCY_MATRIX_HH
