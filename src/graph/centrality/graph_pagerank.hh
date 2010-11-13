// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007-2010 Tiago de Paula Peixoto <tiago@skewed.de>
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

#ifndef GRAPH_PAGERANK_HH
#define GRAPH_PAGERANK_HH

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_util.hh"

namespace graph_tool
{
using namespace std;
using namespace boost;

struct get_pagerank
{
    template <class Graph, class VertexIndex, class RankMap>
    void operator()(Graph& g, VertexIndex vertex_index, RankMap rank,
                    double damping, double epslon, size_t max_iter,
                    size_t& iter) const
    {
        typedef typename property_traits<RankMap>::value_type rank_type;

        RankMap r_temp(vertex_index,num_vertices(g));

        // init ranks
        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i)     \
                schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            rank[v] = 1.0/N;
        }

        rank_type delta = 2*epslon;
        rank_type d = damping;
        iter = 0;
        while (delta >= epslon)
        {
            delta = 0;
            int i, N = num_vertices(g);
            #pragma omp parallel for default(shared) private(i)     \
                schedule(dynamic) reduction(+:delta)
            for (i = 0; i < N; ++i)
            {
                typename graph_traits<Graph>::vertex_descriptor v =
                    vertex(i, g);
                if (v == graph_traits<Graph>::null_vertex())
                    continue;

                rank_type r = 0;
                typename in_edge_iteratorS<Graph>::type e, e_end;
                for (tie(e, e_end) = in_edge_iteratorS<Graph>::get_edges(v, g);
                     e != e_end; ++e)
                {
                    typename graph_traits<Graph>::vertex_descriptor s =
                        source(*e, g);
                    r += get(rank, s) / out_degree(s, g);
                }
                put(r_temp, v, (1.0 - d) + d * r);

                delta += abs(get(r_temp, v) - get(rank,v));
            }
            swap(r_temp, rank);
            ++iter;
            if (max_iter > 0 && iter == max_iter)
                break;
        }

        if (iter % 2 != 0)
        {
            #pragma omp parallel for default(shared) private(i)     \
                schedule(dynamic)
            for (int i = 0; i < N; ++i)
            {
                typename graph_traits<Graph>::vertex_descriptor v =
                    vertex(i, g);
                if (v == graph_traits<Graph>::null_vertex())
                    continue;
                put(rank, v, get(r_temp, v));
            }
        }
    }
};

}
#endif // GRAPH_PAGERANK_HH
