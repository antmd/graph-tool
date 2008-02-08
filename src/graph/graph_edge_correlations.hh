
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

#ifndef GRAPH_EDGE_CORRELATIONS_HH
#define GRAPH_EDGE_CORRELATIONS_HH

#include <tr1/unordered_set>
#include <boost/tuple/tuple.hpp>
#include "shared_map.hh"

namespace graph_tool
{

using namespace std;
using namespace boost;

// retrieves the generalized vertex-edge-vertex correlation histogram

template <class Hist>
struct get_edge_correlation_histogram
{
    get_edge_correlation_histogram(Hist& hist): _hist(hist) {}

    template <class Graph, class DegreeSelector1, class EdgeProperty,
              class DegreeSelector2>
    void operator()(const Graph* gp, DegreeSelector1 deg1, 
                    EdgeProperty edge_prop, DegreeSelector2 deg2) const
    {
        const Graph& g = *gp;
        SharedMap<Hist> s_hist(_hist);

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) \
            firstprivate(s_hist) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            typename Hist::key_type key;
            tuples::get<0>(key) = deg1(v,g);

            typename graph_traits<Graph>::out_edge_iterator e, e_begin, e_end;
            tie(e_begin,e_end) = out_edges(v,g);
            for(e = e_begin; e != e_end; ++e)
            {
                tuples::get<1>(key) = get(edge_prop, *e);
                tuples::get<2>(key) = deg2(target(*e,g),g);
                s_hist[key]++;
            }
        }
        s_hist.Gather();
    }

    Hist& _hist;
};

} // graph_tool namespace

#endif // GRAPH_EDGE_CORRELATIONS_HH
