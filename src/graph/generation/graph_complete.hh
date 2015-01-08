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

#ifndef GRAPH_COMPLETE_HH
#define GRAPH_COMPLETE_HH

#include <iostream>
#include "graph_util.hh"

namespace graph_tool
{
using namespace std;
using namespace boost;

struct get_complete
{
    template <class Graph>
    void operator()(Graph& g, size_t N, bool directed, bool self_loops) const
    {
        for (size_t i = 0; i < N; ++i)
            add_vertex(g);

        for (size_t i = 0; i < N; ++i)
        {
            for (size_t j = directed ? 0 : i; j < N; ++j)
            {
                if (!self_loops && j == i)
                    continue;
                add_edge(vertex(i, g),
                         vertex(j, g), g);
            }
        }
    }
};

struct get_circular
{
    template <class Graph>
    void operator()(Graph& g, size_t N, size_t k, bool directed,
                    bool self_loops) const
    {
        for (size_t i = 0; i < N; ++i)
            add_vertex(g);

        for (size_t i = 0; i < N; ++i)
        {
            for (size_t j = i; j < i + k + 1; ++j)
            {
                if (!self_loops && j == i)
                    continue;
                add_edge(vertex(i, g),
                         vertex(j % N, g), g);
                if (directed && j != i)
                    add_edge(vertex(j % N, g),
                             vertex(i, g), g);
            }
        }
    }
};

} // namespace graph_tool

#endif // GRAPH_COMPLETE_HH
