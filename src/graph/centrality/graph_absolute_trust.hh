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
#include <tr1/tuple>
#include <algorithm>

#include "minmax.hh"

#include <iostream>

namespace graph_tool
{
using namespace std;
using namespace boost;
using std::tr1::get;
using std::tr1::tuple;

template <class Path>
struct path_cmp
{
    path_cmp(vector<Path>& paths): _paths(paths) {}
    vector<Path>& _paths;

    typedef size_t first_argument_type;
    typedef size_t second_argument_type;
    typedef bool result_type;
    inline bool operator()(size_t a, size_t b)
    {
        if (get<0>(_paths[a]).second == get<0>(_paths[b]).second)
            return get<1>(_paths[a]).size() > get<1>(_paths[b]).size();
        return get<0>(_paths[a]).second < get<0>(_paths[b]).second;
    }
};

struct get_absolute_trust
{
    template <class Graph, class VertexIndex, class EdgeIndex, class TrustMap,
              class InferredTrustMap>
    void operator()(Graph& g, VertexIndex vertex_index, EdgeIndex edge_index,
                    size_t max_edge_index, int64_t source, TrustMap c,
                    InferredTrustMap t, size_t n_paths, bool reversed) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;
        typedef typename property_traits<TrustMap>::value_type c_type;
        typedef typename property_traits<InferredTrustMap>::value_type::
            value_type t_type;

        // the path type: the first value is the (trust,weight) pair, the second
        // the set of vertices in the path and the third is the list of edges,
        // in the path sequence.
        typedef tuple<pair<t_type, t_type>, tr1::unordered_set<vertex_t>,
                      vector<edge_t> > path_t;

        int i, N = (source == -1) ? num_vertices(g) : source + 1;
        #pragma omp parallel for default(shared) private(i) schedule(dynamic)
        for (i= (source == -1) ? 0 : source; i < N; ++i)
        {
            vertex_t v = vertex(i, g);
            t[v].resize(num_vertices(g));

            // path priority queue
            vector<path_t> paths(1);
            typedef double_priority_queue<size_t, path_cmp<path_t> > queue_t;
            queue_t queue = queue_t(path_cmp<path_t>(paths));
            get<0>(paths.back()).first = get<0>(paths.back()).second = 1;
            get<1>(paths.back()).insert(v);
            queue.push(0);

            // this is the actual queue of paths which will be used to compute
            // the trust values
            queue_t final_queue = queue_t(path_cmp<path_t>(paths));

            // store all paths which reach a given vertex
            unchecked_vector_property_map<tr1::unordered_set<size_t>,
                                          VertexIndex>
                path_map(vertex_index, num_vertices(g));

            while (!queue.empty())
            {
                size_t pi = queue.top();
                queue.pop_top();
                vertex_t w;

                // push queue top into final queue
                if (get<2>(paths[pi]).size() > 0)
                {
                    w = target(get<2>(paths[pi]).back(), g);
                    final_queue.push(pi);

                    // augment path map
                    path_map[target(get<2>(paths[pi]).back(),g)].insert(pi);
                }
                else
                {
                    w = v; // the first path
                }

                // if maximum size is reached, remove the bottom
                if ((n_paths > 0) && (final_queue.size() > n_paths))
                {
                    size_t bi = final_queue.bottom();
                    final_queue.pop_bottom();

                    // remove path from path map
                    path_map[boost::target(get<2>(paths[bi]).back(),g)].
                        erase(bi);

                    if (bi == pi)
                        continue;
                }

                // augment paths and put them in the queue
                typename graph_traits<Graph>::out_edge_iterator e, e_end;
                for (tie(e, e_end) = out_edges(w, g); e != e_end; ++e)
                {
                    vertex_t a = target(*e, g);
                    // no loops
                    if (get<1>(paths[pi]).find(a) == get<1>(paths[pi]).end())
                    {
                        // only follow non-zero paths
                        if (c[*e] > 0 || (reversed &&
                                          get<1>(paths[pi]).size() == 1))
                        {
                            size_t npi;
                            paths.push_back(paths[pi]); // clone last path
                            npi = paths.size()-1;

                            path_t& np = paths[npi]; // new path

                            if (!reversed)
                            {
                                // path weight
                                get<0>(np).second = get<0>(np).first;

                                // path value
                                get<0>(np).first *= c[*e];
                            }
                            else
                            {
                                if (get<1>(np).size() > 1)
                                    get<0>(np).second *= c[*e];
                                get<0>(np).first *= c[*e];
                            }
                            get<1>(np).insert(a);
                            get<2>(np).push_back(*e);

                            // keep following paths only if there is a chance
                            // they will make it into the final queue
                            if ((n_paths > 0 && final_queue.size() < n_paths) ||
                                (final_queue.size() == 0 ||
                                 (get<0>(np).second >=
                                  get<0>(paths[final_queue.bottom()]).second)))
                                queue.push(npi);
                            else
                                paths.pop_back();
                        }
                    }
                }
            }

            unchecked_vector_property_map<t_type, VertexIndex>
                weight_sum(vertex_index, num_vertices(g));

            // paths which were already calculated and can be skipped
            tr1::unordered_set<size_t> skip;

            // calculate trust from paths in the final queue
            while (!final_queue.empty())
            {
                size_t pi = final_queue.top();
                final_queue.pop_top();
                path_t& p = paths[pi];

                if (skip.find(pi) != skip.end())
                    continue;

                // mark edges which are to be considered
                unchecked_vector_property_map<uint8_t, EdgeIndex>
                    mark(edge_index, max_edge_index+1);
                tr1::unordered_set<size_t>& apaths =
                    path_map[target(get<2>(p).back(), g)]; // all paths with the
                                                           // same final target

                tr1::unordered_set<size_t> vlist; // all vertices involved
                for (typeof(apaths.begin()) iter = apaths.begin();
                     iter != apaths.end(); ++iter)
                {
                    for (size_t j = 0; j < get<2>(paths[*iter]).size(); ++j)
                    {
                        edge_t e = get<2>(paths[*iter])[j];
                        mark[e] = 1;
                        vlist.insert(target(e,g));
                    }
                    vlist.insert(boost::source(get<2>(paths[*iter]).front(),g));
                }

                // compute out trust
                unchecked_vector_property_map<t_type, VertexIndex>
                    out_trust(vertex_index, num_vertices(g));
                for (typeof(vlist.begin()) viter = vlist.begin();
                     viter != vlist.end(); ++viter)
                {
                    vertex_t u = *viter;
                    if (!reversed)
                    {
                        typename graph_traits<Graph>::out_edge_iterator e,e_end;
                        for (tie(e, e_end) = out_edges(u, g); e != e_end; ++e)
                            if (mark[*e] == 1)
                                out_trust[u] += c[*e];
                    }
                    else
                    {
                        // if reversed, use "in-trust" instead
                        typename in_edge_iteratorS<Graph>::type e, e_end;
                        for (tie(e, e_end) =
                                 in_edge_iteratorS<Graph>::get_edges(v, g);
                             e != e_end; ++e)
                        {
                            if (mark[*e] == 1)
                                out_trust[u] += c[*e];
                        }
                    }
                }

                for (typeof(apaths.begin()) iter = apaths.begin();
                     iter != apaths.end(); ++iter)
                {
                    size_t pi = *iter;
                    path_t& p = paths[pi];

                    // calculate the trust value and weight of the path
                    t_type w = 1, val = 1;
                    for (size_t i = 0; i < get<2>(p).size(); ++i)
                    {
                        edge_t e = get<2>(p)[i];
                        vertex_t u = (!reversed) ?
                            boost::source(e,g) : target(e,g);
                        if (out_trust[u] > 0)
                        {
                            if ((!reversed && i < get<2>(p).size()-1) ||
                                (reversed && i > 0))
                                w *= c[e]*c[e]/out_trust[u];
                            else
                                w *= c[e]/out_trust[u];
                        }
                        val *= c[e];
                    }
                    vertex_t u = target(get<2>(p).back(), g);
                    weight_sum[u] += w;
                    t[v][u] += w*val;

                    skip.insert(pi);
                }
            }

            int j, N = num_vertices(g);
            #pragma omp parallel for default(shared) private(j) \
                schedule(dynamic)
            for (j = 0; j < N; ++j)
            {
                vertex_t w = vertex(i, g);
                if (w == graph_traits<Graph>::null_vertex())
                    continue;
                if (weight_sum[w] > 0)
                    t[v][w] /= weight_sum[v];
            }
        }
    }

};

}

#endif
