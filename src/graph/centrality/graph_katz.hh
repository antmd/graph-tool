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

#ifndef GRAPH_EIGENVECTOR_HH
#define GRAPH_EIGENVECTOR_HH

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_util.hh"

#ifndef __clang__
#include <ext/numeric>
using __gnu_cxx::power;
#else
template <class Value>
Value power(Value value, int n)
{
    return pow(value, n);
}
#endif

namespace graph_tool
{
using namespace std;
using namespace boost;

struct get_katz
{
    template <class Graph, class VertexIndex, class WeightMap,
              class CentralityMap, class PersonalizationMap>
    void operator()(Graph& g, VertexIndex vertex_index, WeightMap w,
                    CentralityMap c, PersonalizationMap beta, long double alpha,
                    long double epsilon, size_t max_iter) const
    {
        typedef typename property_traits<WeightMap>::value_type c_type;
        typedef typename property_traits<CentralityMap>::value_type t_type;

        CentralityMap c_temp(vertex_index, num_vertices(g));

        t_type delta = epsilon + 1;
        size_t iter = 0;
        int i, N = num_vertices(g);
        while (delta >= epsilon)
        {
            #pragma omp parallel for default(shared) private(i) \
                schedule(runtime) if (N > 100)
            for (i = 0; i < N; ++i)
            {
                typename graph_traits<Graph>::vertex_descriptor v =
                    vertex(i, g);
                if (v == graph_traits<Graph>::null_vertex())
                    continue;

                c_temp[v] = get(beta, v);
                typename in_or_out_edge_iteratorS<Graph>::type e, e_end;
                for (tie(e, e_end) = in_or_out_edge_iteratorS<Graph>::get_edges(v, g);
                     e != e_end; ++e)
                {
                    typename graph_traits<Graph>::vertex_descriptor s;
                    if (is_directed::apply<Graph>::type::value)
                        s = source(*e, g);
                    else
                        s = target(*e, g);
                    c_temp[v] += alpha * get(w, *e) * c[s];
                }
            }

            delta = 0;
            #pragma omp parallel for default(shared) private(i) \
                schedule(runtime) if (N > 100) reduction(+:delta)
            for (i = 0; i < N; ++i)
            {
                typename graph_traits<Graph>::vertex_descriptor v =
                    vertex(i, g);
                if (v == graph_traits<Graph>::null_vertex())
                    continue;
                delta += abs(c_temp[v] - c[v]);
            }
            swap(c_temp, c);

            ++iter;
            if (max_iter > 0 && iter == max_iter)
                break;
        }

        if (iter % 2 != 0)
        {
            #pragma omp parallel for default(shared) private(i)     \
                schedule(runtime) if (N > 100)
            for (i = 0; i < N; ++i)
            {
                typename graph_traits<Graph>::vertex_descriptor v =
                    vertex(i, g);
                if (v == graph_traits<Graph>::null_vertex())
                    continue;
                c_temp[v] = c[v];
            }
        }
    }
};

}

#endif
