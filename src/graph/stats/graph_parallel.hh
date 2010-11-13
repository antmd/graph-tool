// Copyright (C) 2008  Tiago de Paula Peixoto <tiago@skewed.de>
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

#if (GCC_VERSION >= 40400)
#   include <tr1/unordered_set>
#else
#   include <boost/tr1/unordered_set.hpp>
#endif
#include "graph_util.hh"

namespace graph_tool
{
using namespace std;
using namespace boost;

// label parallel edges in the order they are found, starting from 1
struct label_parallel_edges
{
    template <class Graph, class EdgeIndexMap, class ParallelMap>
    void operator()(const Graph& g, EdgeIndexMap edge_index,
                    ParallelMap parallel) const
    {
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            tr1::unordered_set<edge_t,DescriptorHash<EdgeIndexMap> >
                p_edges(0, DescriptorHash<EdgeIndexMap>(edge_index));

            typename graph_traits<Graph>::out_edge_iterator e1, e2,
                e_end1, e_end2;
            for (tie(e1, e_end1) = out_edges(v, g); e1 != e_end1; ++e1)
            {
                if (p_edges.find(*e1) != p_edges.end())
                    continue;

                // do not visit edges twice in undirected graphs
                if (!is_directed::apply<Graph>::type::value &&
                    target(*e1, g) < v)
                    continue;

                size_t n = 0;
                put(parallel, *e1, n);
                for (tie(e2, e_end2) = out_edges(v, g); e2 != e_end2; ++e2)
                    if (*e2 != *e1 && target(*e1, g) == target(*e2, g))
                    {
                        put(parallel, *e2, ++n);
                        p_edges.insert(*e2);
                    }
            }
        }
    }
};

// label self loops edges in the order they are found, starting from 1
struct label_self_loops
{
    template <class Graph, class EdgeIndexMap, class SelfMap>
    void operator()(const Graph& g, EdgeIndexMap edge_index,
                    SelfMap self) const
    {
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) schedule(dynamic)
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
                    put(self, *e, n++);
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
