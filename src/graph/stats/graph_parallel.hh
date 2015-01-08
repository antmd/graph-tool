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

#ifndef GRAPH_PARALLEL_HH
#define GRAPH_PARALLEL_HH

#include <unordered_map>

#ifdef HAVE_SPARSEHASH
#include SPARSEHASH_INCLUDE(dense_hash_map)
#endif

#include "graph_util.hh"

namespace graph_tool
{
using namespace std;
using namespace boost;

// label parallel edges in the order they are found, starting from 1
struct label_parallel_edges
{
    template <class Graph, class ParallelMap>
    void operator()(const Graph& g, ParallelMap parallel, bool mark_only,
                    bool count_all) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;
        typename property_map<Graph, edge_index_t>::type eidx = get(edge_index, g);

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) schedule(runtime) if (N > 100)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

#ifdef HAVE_SPARSEHASH
            google::dense_hash_map<vertex_t, edge_t> vset;
            vset.set_empty_key(graph_traits<Graph>::null_vertex());
            google::dense_hash_map<size_t, bool> self_loops;
            self_loops.set_empty_key(numeric_limits<size_t>::max());
#else
            unordered_map<vertex_t, edge_t> vset;
            unordered_map<size_t, bool> self_loops;
#endif

            typename graph_traits<Graph>::out_edge_iterator e, e_end;
            for (tie(e, e_end) = out_edges(v, g); e != e_end; ++e)
            {
                vertex_t u = target(*e, g);

                // do not visit edges twice in undirected graphs
                if (!is_directed::apply<Graph>::type::value && u < v)
                    continue;

                if (u == v)
                {
                    if (self_loops[eidx[*e]])
                        continue;
                    self_loops[eidx[*e]] = true;
                }

                typeof(vset.begin()) iter = vset.find(u);
                if (iter == vset.end())
                {
                    vset[u] = *e;
                }
                else
                {
                    if (mark_only)
                    {
                        parallel[*e] = true;
                    }
                    else
                    {
                        parallel[*e] = parallel[iter->second] + 1;
                        vset[u] = *e;
                    }
                }
            }
        }
    }
};

// label self loops edges in the order they are found, starting from 1
struct label_self_loops
{
    template <class Graph, class SelfMap>
    void operator()(const Graph& g, SelfMap self, bool mark_only) const
    {
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) schedule(runtime) if (N > 100)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            size_t n = 1;
            typename graph_traits<Graph>::out_edge_iterator e, e_end;
            for (tie(e, e_end) = out_edges(v, g); e != e_end; ++e)
            {
                if (target(*e, g) == v)
                    put(self, *e, mark_only ? 1 : n++);
                else
                    put(self, *e, 0);
            }
        }
    }
};

// remove edges with label larger than 0
struct remove_labeled_edges
{
    template <class Graph, class LabelMap>
    void operator()(Graph& g, LabelMap label) const
    {
        int i, N = num_vertices(g);

        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            typedef typename graph_traits<Graph>::edge_descriptor edge_t;
            vector<edge_t> r_edges;

            typename graph_traits<Graph>::out_edge_iterator e, e_end;
            for (tie(e, e_end) = out_edges(v, g); e != e_end; ++e)
            {
                if (label[*e] > 0)
                    r_edges.push_back(*e);
            }
            for (size_t j = 0; j < r_edges.size(); ++j)
                remove_edge(r_edges[j], g);
        }
    }
};

} // graph_tool namespace

#endif //GRAPH_PARALLEL_HH
