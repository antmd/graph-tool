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
        if (get<0>(_paths[a]).first == get<0>(_paths[b]).first)
            return get<1>(_paths[a]).size() > get<1>(_paths[b]).size();
        return get<0>(_paths[a]).first < get<0>(_paths[b]).first;
    }
};

struct get_absolute_trust
{
    template <class Graph, class VertexIndex, class TrustMap,
              class InferredTrustMap>
    void operator()(Graph& g, VertexIndex vertex_index, int64_t source,
                    TrustMap c, InferredTrustMap t, size_t n_paths,
                    double epsilon, bool reversed) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        typedef typename property_traits<TrustMap>::value_type c_type;
        typedef typename property_traits<InferredTrustMap>::value_type
            ::value_type t_type;

        typedef tuple<pair<t_type, t_type>, tr1::unordered_set<vertex_t>,
                      size_t > path_t;

        double delta = epsilon+1;
        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) schedule(dynamic)
        for (i = (source == -1) ? 0 : source;
             i < ((source == -1) ? N : source + 1); ++i)
        {
            vertex_t v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            // path priority queue
            vector<path_t> paths(1);
            vector<size_t> free_indexes;
            typedef double_priority_queue<size_t, path_cmp<path_t> > queue_t;
            queue_t queue = queue_t(path_cmp<path_t>(paths));

            get<0>(paths.back()).first = get<0>(paths.back()).second = 1;
            get<1>(paths.back()).insert(v);
            get<2>(paths.back()) = vertex_index[v];
            queue.push(0);

            t[v].resize(num_vertices(g));
            unchecked_vector_property_map<t_type, VertexIndex>
                sum_weight(vertex_index, num_vertices(g));

            size_t count = 1;
            while (!queue.empty())
            {
                size_t pi = queue.top();
                queue.pop_top();
                vertex_t w = vertex(get<2>(paths[pi]), g);

                typename graph_traits<Graph>::out_edge_iterator e, e_end;
                for (tie(e, e_end) = out_edges(w, g); e != e_end; ++e)
                {
                    vertex_t a = target(*e, g);
                    // no loops
                    if (get<1>(paths[pi]).find(a) == get<1>(paths[pi]).end())
                    {
                        // new path;  only follow non-zero paths
                        t_type old = (sum_weight[a] > 0) ?
                            t[v][a]/sum_weight[a] : 0;
                        if (c[*e] > 0 || (reversed &&
                                          get<1>(paths[pi]).size() == 1))
                        {
                            size_t npi;
                            if (free_indexes.empty())
                            {
                                paths.push_back(paths[pi]); // clone last path
                                npi = paths.size()-1;
                            }
                            else
                            {
                                npi = free_indexes.back();
                                free_indexes.pop_back();
                                paths[npi] = paths[pi];
                            }

                            path_t& np = paths[npi];
                            if (!reversed)
                            {
                                get<0>(np).second = get<0>(np).first;
                                get<0>(np).first *= c[*e];
                            }
                            else
                            {
                                if (get<1>(np).size() > 1)
                                    get<0>(np).second *= c[*e];
                                get<0>(np).first *= c[*e];
                            }
                            get<1>(np).insert(a);
                            get<2>(np) = vertex_index[a];

                            t[v][a] += get<0>(np).first*get<0>(np).second;
                            sum_weight[a] += get<0>(np).second;
                            queue.push(npi);
                            if (n_paths > 0 && queue.size() > n_paths)
                            {
                                size_t bi = queue.bottom();
                                queue.pop_bottom();
                                free_indexes.push_back(bi);
                            }
                        }
                        else
                        {
                            sum_weight[a] +=  get<0>(paths[pi]).first;
                        }
                        if (sum_weight[a] > 0)
                            delta += abs(old-t[v][a]/sum_weight[a]);
                    }
                }
                free_indexes.push_back(pi);

                if ((count % N) == 0)
                {
                    if (delta < epsilon)
                        break;
                    else
                        delta = 0;
                }
                count++;
            }

            typename graph_traits<Graph>::vertex_iterator w, w_end;
            for (tie(w, w_end) = vertices(g); w != w_end; ++w)
            {
                if (sum_weight[*w] > 0)
                    t[v][*w] /= sum_weight[*w];
            }
        }
    }
};

}

#endif
