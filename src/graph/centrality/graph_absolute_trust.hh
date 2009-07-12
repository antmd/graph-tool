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

#include <tr1/unordered_set>
#include <tr1/random>
#include <boost/functional/hash.hpp>

#include <iostream>

namespace graph_tool
{
using namespace std;
using namespace boost;

struct get_absolute_trust
{
    template <class Graph, class VertexIndex, class TrustMap,
              class InferredTrustMap>
    void operator()(Graph& g, VertexIndex vertex_index, int64_t source,
                    TrustMap c, InferredTrustMap t, double epslon,
                    pair<size_t,size_t> iter_range, bool reversed,
                    rng_t& rng) const
    {
        typedef typename property_traits<TrustMap>::value_type c_type;
        typedef typename property_traits<InferredTrustMap>::value_type
            ::value_type t_type;

        size_t min_iter = iter_range.first;
        size_t max_iter = iter_range.second;

        unchecked_vector_property_map<vector<t_type>, VertexIndex>
            t_count(vertex_index, num_vertices(g));
        unchecked_vector_property_map<vector<size_t>, VertexIndex>
            v_mark(vertex_index, num_vertices(g));

        // init inferred trust t
        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) schedule(dynamic)
        for (i = (source == -1) ? 0 : source;
             i < ((source == -1) ? N : source + 1); ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v =
                vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            t[v].resize(N);
            t_count[v].resize(N,0);
            v_mark[v].resize(N,0);
        }

        #pragma omp parallel for default(shared) private(i) schedule(dynamic)
        for (i = (source == -1) ? 0 : source;
             i < ((source == -1) ? N : source + 1); ++i)
        {
            // walk hash set
            tr1::unordered_set<size_t> path_set;
            tr1::uniform_int<size_t>
                random_salt(0, numeric_limits<size_t>::max());
            size_t salt = random_salt(rng);

            t_type delta = epslon + 1;
            size_t iter = 0;
            while (delta > epslon || iter < min_iter)
            {
                delta = 0;
                typename graph_traits<Graph>::vertex_descriptor v =
                    vertex(i, g);
                if (v == graph_traits<Graph>::null_vertex())
                    continue;

                typename graph_traits<Graph>::vertex_descriptor pos = v;
                t_type pos_t = 1.0;
                t_type pweight = 1.0;
                v_mark[v][vertex_index[v]] = iter+1;

                size_t path_hash = salt;
                hash_combine(path_hash, i);

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
                            if (c[*e] < 0)
                                continue;
                            out_es.push_back(*e);
                            if (out_prob.empty())
                                out_prob.push_back(c[*e]+0.1);
                            else
                                out_prob.push_back(out_prob.back()+c[*e]+0.1);
                        }
                    }
                    if (!out_es.empty())
                    {
                        // select edge according to its probability
                        typename graph_traits<Graph>::edge_descriptor e;
                        typedef tr1::uniform_real<t_type> dist_t;
                        tr1::variate_generator<rng_t, dist_t>
                            random(rng, dist_t(0, out_prob.back()));

                        t_type u;
                        {
                            #pragma omp critical
                            u = random();
                        }
                        e = out_es[lower_bound(out_prob.begin(),
                                               out_prob.end(), u) -
                                   out_prob.begin()];

                        // new position in random walk
                        pos = target(e,g);
                        size_t posi = vertex_index[pos];

                        //update current path trust
                        pos_t *= c[e];

                        if (reversed && boost::source(e, g) != v)
                            pweight *= c[e];

                        // get path hash
                        hash_combine(path_hash, posi);
                        if (path_set.find(path_hash) == path_set.end())
                        {
                            // if new path, modify vertex trust score
                            path_set.insert(path_hash);

                            t_type old = 0;
                            if (t_count[v][posi] > 0)
                                old = t[v][posi]/t_count[v][posi];
                            t[v][posi] += pos_t*pweight;
                            t_count[v][posi] += pweight;
                            delta += abs(old - t[v][posi]/t_count[v][posi]);
                        }

                        if (!reversed)
                            pweight *= c[e];

                        // stop if weight drops to zero
                        if (pweight == 0)
                            break;

                        v_mark[v][posi] = iter+1; // mark vertex
                    }
                }
                while (!out_es.empty());
                ++iter;
                if (max_iter > 0 && iter == max_iter)
                    break;
             }
        }

        #pragma omp parallel for default(shared) private(i) schedule(dynamic)
        for (i = (source == -1) ? 0 : source;
             i < ((source == -1) ? N : source + 1); ++i)
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
