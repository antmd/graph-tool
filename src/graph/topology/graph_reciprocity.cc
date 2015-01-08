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


#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

struct get_reciprocity
{
    template <class Graph>
    void operator()(const Graph& g, double& reciprocity) const
    {
        size_t L = 0, Lbd = 0;

        int i, NV = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) reduction(+:L,Lbd) \
            schedule(runtime) if (NV > 100)
        for (i = 0; i < NV; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            typename graph_traits<Graph>::out_edge_iterator e, e_begin, e_end;
            tie(e_begin,e_end) = out_edges(v,g);
            for(e = e_begin; e != e_end; ++e)
            {
                typename graph_traits<Graph>::vertex_descriptor t;
                t = target(*e, g);

                typename graph_traits<Graph>::adjacency_iterator a, a_end;
                for (tie(a, a_end) = adjacent_vertices(t, g); a != a_end; ++a)
                    if (*a == v)
                    {
                        Lbd += 1;
                        break;
                    }
                L++;
            }
        }

        reciprocity = Lbd / double(L);
    }
};

double reciprocity(GraphInterface& gi)
{
    double reciprocity;
    run_action<>()(gi, std::bind(get_reciprocity(), placeholders::_1,
                                 std::ref(reciprocity)))();
    return reciprocity;
}
