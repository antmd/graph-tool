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

#ifndef GRAPH_KCORE_HH
#define GRAPH_KCORE_HH

namespace graph_tool
{
using namespace std;
using namespace boost;

struct kcore_decomposition
{
    template <class Graph, class VertexIndex, class CoreMap, class DegSelector>
    void operator()(Graph& g, VertexIndex vertex_index, CoreMap core_map,
                    DegSelector degS) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

        unchecked_vector_property_map<size_t, VertexIndex> deg(vertex_index,
                                                               num_vertices(g));
        unchecked_vector_property_map<size_t, VertexIndex> pos(vertex_index,
                                                               num_vertices(g));
        vector<vector<vertex_t> > bins;

        typename graph_traits<Graph>::vertex_iterator vi, vi_end;
        for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
        {
            size_t k = degS(*vi, g);
            deg[*vi] = k;
            if (k >= bins.size())
                bins.resize(k + 1);
            bins[k].push_back(*vi);
            pos[*vi] = bins[k].size() - 1;
        }

        for (size_t k = 0; k < bins.size(); ++k)
        {
            while (!bins[k].empty())
            {
                vertex_t v = bins[k].back();
                bins[k].pop_back();
                core_map[v] = k;
                typename graph_traits<Graph>::out_edge_iterator e, e_end;
                for (tie(e, e_end) = out_edges(v, g); e != e_end; ++e)
                {
                    vertex_t u = target(*e, g);
                    if (deg[u] > deg[v])
                    {
                        size_t ku = deg[u];
                        vertex_t w = bins[ku].back();
                        pos[w] = pos[u];
                        bins[ku][pos[w]] = w;
                        bins[ku].pop_back();
                        bins[ku - 1].push_back(u);
                        pos[u] = bins[ku - 1].size() - 1;
                        --deg[u];
                    }
                }
            }
        }
    }
};

} // graph_tool namespace

#endif // GRAPH_KCORE_HH
