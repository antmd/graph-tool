// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007  Tiago de Paula Peixoto <tiago@forked.de>
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

#ifndef GRAPH_ABSOLUTE_TRUST_HH
#define GRAPH_ABSOLUTE_TRUST_HH

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_util.hh"

#include <boost/random/uniform_int.hpp>

namespace graph_tool
{
using namespace std;
using namespace boost;

struct get_eigentrust
{
    template <class Graph, class VertexIndex, class EdgeIndex, class TrustMap,
              class InferredTrustMap>
    void operator()(Graph* gp, VertexIndex vertex_index, EdgeIndex edge_index,
                    TrustMap c, InferredTrustMap t, double epslon,
                    size_t max_iter, rng_t& rng)
        const
    {
        Graph& g = *gp;
        typedef typename property_traits<TrustMap>::value_type c_type;
        typedef typename property_traits<InferredTrustMap>::value_type
            ::value_type t_type;

        unchecked_vector_property_map<vector<size_t>, VertexIndex>
            t_count(vertex_index, num_vertices(g));
        unchecked_vector_property_map<vector<size_t>, VertexIndex>
            v_mark(vertex_index, num_vertices(g));

        // init inferred trust t
        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i)     \
                schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            t[v].resize(N);
            t_count[v].resize(N,0);
            v_mark[v].resize(N,0);
        }

        t_type delta = 2*epslon;
        size_t iter = 0;
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

                typename graph_traits<Graph>::vertex_descriptor pos = v;
                t_type pos_t = 1.0;
                v_mark[v][vertex_index[v]] = iter+1;

                // start a self-avoiding walk from vertex v
                vector<typename graph_traits<Graph>::edge_descriptor> out_es;
                vector<t_type> out_prob;
                do
                {
                    out_es.clear();
                    out_prob.clear();

                    // obtain list of candidate edges to follow
                    typename graph_traits<Graph>::out_edge_iterator e, e_end;
                    for (tie(e, e_end) = out_edges(pos, g); e != e_end; ++e)
                    {
                        typename graph_traits<Graph>::vertex_descriptor t =
                            target(*e,g);
                        if (v_mark[v][vertex_index[t]] <= iter)
                        {
                            out_es.push_back(*e);
                            if (out_prob.empty())
                                out_prob.push_back(c[*e]);
                            else
                                out_prob.push_back(out_prob.back()+c[*e]);
                        }
                    }
                    if (!out_es.empty())
                    {
                        // select edge according to its probability
                        typename graph_traits<Graph>::edge_descriptor e;
                        uniform_real<t_type> random(0,out_prob.back());

                        t_type u;
                        {
                            #pragma omp critical
                            u = random(rng);
                        }
                        e = out_es[lower_bound(out_prob.begin(),
                                               out_prob.end(), u) -
                                   out_prob.begin()];
                        pos = target(e,g);
                        size_t posi = vertex_index[pos];

                        //update current path trust, and update new vertex
                        pos_t *= c[e];
                        t_type old = 0;
                        if (t_count[v][posi] > 0)
                            old = t[v][posi]/t_count[v][posi];
                        t[v][posi] += pos_t;
                        t_count[v][posi]++;
                        delta += abs(old - t[v][posi]/t_count[v][posi]);

                        v_mark[v][posi] = iter+1; // mark vertex
                    }
                }
                while (!out_es.empty());
            }
            ++iter;
            if (max_iter > 0 && iter == max_iter)
                break;
        }

        #pragma omp parallel for default(shared) private(i)     \
                schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            for (size_t j = 0; j < N; ++j)
                if (t_count[v][j] > 0)
                    t[v][j] /= t_count[v][j];
        }
    }
};

}

#endif
