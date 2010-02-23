// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007-2010 Tiago de Paula Peixoto <tiago@forked.de>
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

#ifndef GRAPH_TRUST_HH
#define GRAPH_TRUST_HH

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_util.hh"

namespace graph_tool
{
using namespace std;
using namespace boost;

struct get_eigentrust
{
    template <class Graph, class VertexIndex, class EdgeIndex, class TrustMap,
              class InferredTrustMap>
    void operator()(Graph& g, VertexIndex vertex_index,
                    EdgeIndex edge_index, TrustMap c, InferredTrustMap t,
                    double epslon, size_t max_iter, size_t& iter) const
    {
        typedef typename property_traits<TrustMap>::value_type c_type;
        typedef typename property_traits<InferredTrustMap>::value_type t_type;

        InferredTrustMap t_temp(vertex_index, num_vertices(g));

        // Norm c values
        InferredTrustMap c_sum(vertex_index);
        if (typename is_directed::apply<Graph>::type())
        {
            TrustMap c_temp(edge_index, c.get_storage().size());

            int i, N = num_vertices(g);
            #pragma omp parallel for default(shared) private(i)  \
                schedule(dynamic)
            for (i = 0; i < N; ++i)
            {
                typename graph_traits<Graph>::vertex_descriptor v =
                    vertex(i, g);
                if (v == graph_traits<Graph>::null_vertex())
                    continue;

                c_type sum = 0;
                typename graph_traits<Graph>::out_edge_iterator e, e_end;
                for (tie(e, e_end) = out_edges(v, g); e != e_end; ++e)
                    sum += get(c, *e);

                if (sum > 0)
                    for (tie(e, e_end) = out_edges(v, g); e != e_end; ++e)
                        put(c_temp, *e, get(c, *e)/sum);
            }
            c = c_temp;
        }
        else
        {
            c_sum.reserve(num_vertices(g));
            int i, N = num_vertices(g);
            #pragma omp parallel for default(shared) private(i)     \
                schedule(dynamic)
            for (i = 0; i < N; ++i)
            {
                typename graph_traits<Graph>::vertex_descriptor v =
                    vertex(i, g);
                if (v == graph_traits<Graph>::null_vertex())
                    continue;

                c_sum[v] = 0;
                typename graph_traits<Graph>::out_edge_iterator e, e_end;
                for (tie(e, e_end) = out_edges(v, g); e != e_end; ++e)
                    c_sum[v] += c[*e];
            }
        }

        // init inferred trust t
        int i, N = num_vertices(g), V = HardNumVertices()(g);
        #pragma omp parallel for default(shared) private(i)     \
                schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            t[v] = 1.0/V;
        }

        t_type delta = epslon + 1;
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

                t_temp[v] = 0;
                typename in_edge_iteratorS<Graph>::type e, e_end;
                for (tie(e, e_end) = in_edge_iteratorS<Graph>::get_edges(v, g);
                     e != e_end; ++e)
                {
                    typename graph_traits<Graph>::vertex_descriptor s =
                        source(*e,g);
                    if (!typename is_directed::apply<Graph>::type())
                        t_temp[v] += get(c, *e)*t[s]/abs(c_sum[s]);
                    else
                        t_temp[v] += get(c, *e)*t[s];
                }
                delta += abs(t_temp[v] - t[v]);
            }
            swap(t_temp, t);

            ++iter;
            if (max_iter > 0 && iter== max_iter)
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
                t[v] = t_temp[v];
            }
        }
    }
};

}

#endif
