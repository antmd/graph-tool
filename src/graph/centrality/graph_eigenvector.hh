// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007-2011 Tiago de Paula Peixoto <tiago@skewed.de>
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

#include <ext/numeric>
using __gnu_cxx::power;

namespace graph_tool
{
using namespace std;
using namespace boost;

struct get_eigenvector
{
    template <class Graph, class VertexIndex, class EdgeIndex, class WeightMap,
              class CentralityMap>
    void operator()(Graph& g, VertexIndex vertex_index,
                    EdgeIndex edge_index, WeightMap w, CentralityMap c,
                    double epsilon, size_t max_iter, long double& eig) const
    {
        typedef typename property_traits<WeightMap>::value_type c_type;
        typedef typename property_traits<CentralityMap>::value_type t_type;

        CentralityMap c_temp(vertex_index, num_vertices(g));

        // init centrality
        int i, N = num_vertices(g), V = HardNumVertices()(g);
        #pragma omp parallel for default(shared) private(i)     \
                schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            c[v] = 1.0 / V;
        }

        t_type norm = 0;


        t_type delta = epsilon + 1;
        size_t iter = 0;
        while (delta >= epsilon)
        {
            norm = 0;
            int i, N = num_vertices(g);
            #pragma omp parallel for default(shared) private(i) \
                schedule(dynamic) reduction(+:norm)
            for (i = 0; i < N; ++i)
            {
                typename graph_traits<Graph>::vertex_descriptor v =
                    vertex(i, g);
                if (v == graph_traits<Graph>::null_vertex())
                    continue;

                c_temp[v] = 0;
                typename in_edge_iteratorS<Graph>::type e, e_end;
                for (tie(e, e_end) = in_edge_iteratorS<Graph>::get_edges(v, g);
                     e != e_end; ++e)
                {
                    typename graph_traits<Graph>::vertex_descriptor s =
                        source(*e,g);
                    c_temp[v] += get(w, *e) * c[s];
                }
                norm += power(c_temp[v], 2);
            }
            norm = sqrt(norm);

            delta = 0;
            #pragma omp parallel for default(shared) private(i) \
                schedule(dynamic) reduction(+:delta)
            for (i = 0; i < N; ++i)
            {
                typename graph_traits<Graph>::vertex_descriptor v =
                    vertex(i, g);
                if (v == graph_traits<Graph>::null_vertex())
                    continue;
                c_temp[v] /= norm;
                delta += abs(c_temp[v] - c[v]);
            }
            swap(c_temp, c);

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
                c[v] = c_temp[v];
            }
        }

        eig = 1. / norm;
    }
};

}

#endif
