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
                    InferredTrustMap t, boost::tuple<size_t, size_t,
                    double> path_limits, bool reversed) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;
        typedef typename property_traits<TrustMap>::value_type c_type;
        typedef typename property_traits<InferredTrustMap>::value_type::
            value_type t_type;

        size_t n_paths = get<0>(path_limits);
        size_t n_paths_vertex = get<1>(path_limits);
        t_type epsilon = get<2>(path_limits);

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
            vector<size_t> free_paths; // free path indexes

            // this is the actual queue of paths which will be used to compute
            // the trust values
            queue_t final_queue = queue_t(path_cmp<path_t>(paths));

            // store all paths which reach a given vertex
            unchecked_vector_property_map<tr1::unordered_set<size_t>,
                                          VertexIndex>
                path_map(vertex_index, num_vertices(g));

            unchecked_vector_property_map<bool, VertexIndex>
                saturated(vertex_index, num_vertices(g));

            while (!queue.empty())
            {
                size_t pi = queue.top();
                queue.pop_top();
                vertex_t w;

                vertex_t ptarget;
                if (get<2>(paths[pi]).size() > 0)
                    ptarget = target(get<2>(paths[pi]).back(),g);

                // push queue top into final queue
                if (get<2>(paths[pi]).size() > 0 &&
                    (n_paths_vertex == 0 || !saturated[ptarget]))
                {
                    w = target(get<2>(paths[pi]).back(), g);
                    final_queue.push(pi);

                    // augment path map
                    path_map[ptarget].insert(pi);

                    if (path_map[ptarget].size() == n_paths_vertex)
                        saturated[ptarget] = true;
                }
                else
                {
                    w = v; // the first path
                }

                // if maximum size is reached, remove the bottom
                if ((n_paths > 0) && (final_queue.size() > n_paths))
                {
                    size_t bi = final_queue.bottom();
                    ptarget = target(get<2>(paths[bi]).back(),g);
                    if (!saturated[ptarget])
                    {
                        final_queue.pop_bottom();
                        free_paths.push_back(bi);

                        // remove path from path map
                        path_map[ptarget].erase(bi);

                        if (bi == pi)
                            continue;
                    }
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
                            // clone last path
                            if (free_paths.empty())
                            {
                                paths.push_back(paths[pi]);
                                npi = paths.size()-1;
                            }
                            else
                            {
                                npi = free_paths.back();
                                free_paths.pop_back();
                                paths[npi] = paths[pi];
                            }

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
                            {
                                // drop paths with weight smaller than epsilon
                                if (get<0>(np).second > epsilon)
                                    queue.push(npi);
                            }
                            else
                            {
                                if (npi == paths.size() - 1)
                                    paths.pop_back();
                                else
                                    free_paths.push_back(npi);
                            }
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

                tr1::unordered_set<size_t>& apaths =
                    path_map[target(get<2>(p).back(), g)]; // all paths with the
                                                           // same final target

                tr1::unordered_set<size_t> vlist; // all vertices involved

                // compute cumulative edge weight
                tr1::unordered_map<size_t,
                                   tr1::unordered_map<vertex_t, t_type> >
                    cum_weight;

                for (typeof(apaths.begin()) iter = apaths.begin();
                     iter != apaths.end(); ++iter)
                {
                    path_t& path = paths[*iter];
                    size_t path_size = get<2>(path).size();
                    vertex_t ptarget = target(get<2>(path).back(),g);
                    t_type w = 1;
                    for (size_t j = 0; j < path_size; ++j)
                    {
                        edge_t e;
                        if (!reversed)
                        {
                            e = get<2>(path)[path_size - 1 - j];
                            if (j < path_size - 1)
                                w *= c[e];
                        }
                        else
                        {
                            e = get<2>(path)[j];
                            if (j > 0)
                                w *= c[e];
                        }
                        cum_weight[edge_index[e]][ptarget] += w;
                        vlist.insert(target(e,g));
                    }
                    vlist.insert(boost::source(get<2>(paths[*iter]).front(),g));
                }

                // compute out-weight
                tr1::unordered_map<vertex_t,
                                   tr1::unordered_map<vertex_t, t_type> >
                    out_weight;
                for (typeof(vlist.begin()) viter = vlist.begin();
                     viter != vlist.end(); ++viter)
                {
                    vertex_t u = *viter;
                    if (!reversed)
                    {
                        typename graph_traits<Graph>::out_edge_iterator e,e_end;
                        for (tie(e, e_end) = out_edges(u, g); e != e_end; ++e)
                        {
                            size_t ei = edge_index[*e];
                            for (typeof(cum_weight[ei].begin()) witer =
                                     cum_weight[ei].begin();
                                 witer != cum_weight[ei].end(); ++witer)
                                out_weight[u][witer->first] += witer->second;
                        }
                    }
                    else
                    {
                        // if reversed, use "in-trust" instead
                        typename in_edge_iteratorS<Graph>::type e, e_end;
                        for (tie(e, e_end) =
                                 in_edge_iteratorS<Graph>::get_edges(v, g);
                             e != e_end; ++e)
                        {
                            size_t ei = edge_index[*e];
                            for (typeof(cum_weight[ei].begin()) witer =
                                     cum_weight[ei].begin();
                                 witer != cum_weight[ei].end(); ++witer)
                                out_weight[u][witer->first] += witer->second;
                        }
                    }
                }

                for (typeof(apaths.begin()) iter = apaths.begin();
                     iter != apaths.end(); ++iter)
                {
                    size_t pi = *iter;
                    path_t& p = paths[pi];
                    vertex_t ptarget = target(get<2>(paths[*iter]).back(),g);

                    // calculate the trust value and weight of the path
                    t_type w = 1, val = 1;
                    for (size_t i = 0; i < get<2>(p).size(); ++i)
                    {
                        edge_t e = get<2>(p)[i];
                        vertex_t u = (!reversed) ?
                            boost::source(e,g) : target(e,g);
                        if (out_weight[u][ptarget] > 0)
                        {
                            if ((!reversed && i < get<2>(p).size()-1) ||
                                (reversed && i > 0))
                                w *= c[e]*cum_weight[edge_index[e]][ptarget]/
                                    out_weight[u][ptarget];
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
                vertex_t w = vertex(j, g);
                if (w == graph_traits<Graph>::null_vertex())
                    continue;
                if (weight_sum[w] > 0)
                    t[v][w] /= weight_sum[w];
            }
        }
    }

};

}

#endif
