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

#ifndef GRAPH_CORRELATIONS_HH
#define GRAPH_CORRELATIONS_HH

#include <algorithm>
#include <tr1/unordered_set>
#include "shared_map.hh"

namespace graph_tool
{
using namespace std;
using namespace boost;
using namespace boost::lambda;

// retrieves the generalized vertex-vertex correlation histogram

template <class Hist>
struct get_correlation_histogram
{
    get_correlation_histogram(Hist& hist): _hist(hist) {}

    template <class Graph, class DegreeSelector1, class DegreeSelector2,
              class WeightMap>
    void operator()(Graph* gp, DegreeSelector1 deg1, DegreeSelector2 deg2,
                    WeightMap weight) const
    {
        Graph& g = *gp;
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
            key.first = deg1(v, g);
            typename graph_traits<Graph>::out_edge_iterator e, e_end;
            for (tie(e,e_end) = out_edges(v, g); e != e_end; ++e)
            {
                key.second = deg2(target(*e,g),g);
                s_hist[key] +=
                    typename Hist::value_type::second_type(get(weight, *e));
            }
        }
        s_hist.Gather();
    }
    Hist& _hist;
};


} // graph_tool namespace

#endif
