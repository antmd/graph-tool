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

#ifndef GRAPH_PREDECESSOR_HH
#define GRAPH_PREDECESSOR_HH

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_util.hh"

namespace graph_tool
{
using namespace std;
using namespace boost;

struct get_predecessor_graph
{
    template <class Graph, class IndexMap, class PredGraph, class PredMap>
    void operator()(Graph& g, IndexMap vertex_index, PredGraph& pg,
                    PredMap pred_map) const
    {
        unchecked_vector_property_map<size_t,IndexMap>
            index_map(vertex_index, num_vertices(g));

        size_t count = 0;
        typename graph_traits<Graph>::vertex_iterator v,v_end;
        for (tie(v,v_end) = vertices(g); v != v_end; ++v)
        {
            index_map[*v] = count++;
            add_vertex(pg);
        }

        for (tie(v,v_end) = vertices(g); v != v_end; ++v)
        {
            size_t pred_i = get(pred_map, *v);
            if (pred_i >= num_vertices(g))
                continue;

            typename graph_traits<Graph>::vertex_descriptor pred =
                vertex(pred_i, g);
            if (pred == graph_traits<Graph>::null_vertex())
                continue;

            if (pred != *v)
            {
                typename graph_traits<PredGraph>::vertex_descriptor s, t;
                s = vertex(index_map[pred], pg);
                t = vertex(index_map[*v], pg);
                add_edge(s, t, pg);
            }
        }
    }
};

} // graph_tool namespace

#endif // GRAPH_PREDECESSOR_HH
