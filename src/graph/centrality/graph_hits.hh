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

struct get_hits
{
    template <class Graph, class VertexIndex, class WeightMap,
              class CentralityMap>
    void operator()(Graph& g, VertexIndex vertex_index, WeightMap w,
                    CentralityMap x, CentralityMap y, double epsilon,
                    size_t max_iter, long double& eig) const
    {
        typedef typename property_traits<WeightMap>::value_type c_type;
        typedef typename property_traits<CentralityMap>::value_type t_type;

        CentralityMap x_temp(vertex_index, num_vertices(g));
        CentralityMap y_temp(vertex_index, num_vertices(g));

        // init centrality
        int i, N = num_vertices(g), V = HardNumVertices()(g);
        #pragma omp parallel for default(shared) private(i)     \
                schedule(runtime) if (N > 100)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            x[v] = 1.0 / V;
            y[v] = 1.0 / V;
        }

        t_type x_norm = 0;

        t_type delta = epsilon + 1;
        size_t iter = 0;
        while (delta >= epsilon)
        {
            x_norm = 0;
            #pragma omp parallel for default(shared) private(i) \
                schedule(runtime) if (N > 100) reduction(+:x_norm)
            for (i = 0; i < N; ++i)
            {
                typename graph_traits<Graph>::vertex_descriptor v =
                    vertex(i, g);
                if (v == graph_traits<Graph>::null_vertex())
                    continue;

                x_temp[v] = 0;
                typename in_or_out_edge_iteratorS<Graph>::type ie, ie_end;
                for (tie(ie, ie_end) = in_or_out_edge_iteratorS<Graph>::get_edges(v, g);
                     ie != ie_end; ++ie)
                {
                    typename graph_traits<Graph>::vertex_descriptor s =
                        source(*ie, g);
                    if (is_directed::apply<Graph>::type::value)
                        s = source(*ie, g);
                    else
                        s = target(*ie,g);
                    x_temp[v] += get(w, *ie) * y[s];
                }
                x_norm += power(x_temp[v], 2);

                y_temp[v] = 0;
                typename graph_traits<Graph>::out_edge_iterator e, e_end;
                for (tie(e, e_end) = out_edges(v, g); e != e_end; ++e)
                {
                    typename graph_traits<Graph>::vertex_descriptor s = target(*e, g);
                    y_temp[v] += get(w, *e) * x[s];
                }
            }
            x_norm = sqrt(x_norm);

            delta = 0;
            #pragma omp parallel for default(shared) private(i) \
                schedule(runtime) if (N > 100) reduction(+:delta)
            for (i = 0; i < N; ++i)
            {
                typename graph_traits<Graph>::vertex_descriptor v =
                    vertex(i, g);
                if (v == graph_traits<Graph>::null_vertex())
                    continue;
                x_temp[v] /= x_norm;
                delta += abs(x_temp[v] - x[v]);
                delta += abs(y_temp[v] - y[v]);
            }
            swap(x_temp, x);
            swap(y_temp, y);

            ++iter;
            if (max_iter > 0 && iter== max_iter)
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
                x[v] = x_temp[v];
                y[v] = y_temp[v];
            }
        }

        eig = x_norm;
    }
};

}

#endif
