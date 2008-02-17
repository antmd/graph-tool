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

#ifndef GRAPH_CORRELATIONS_NEIGHBOURS_HH
#define GRAPH_CORRELATIONS_NEIGHBOURS_HH

#include "graph_filtering.hh"
#include "shared_map.hh"

#include <tr1/unordered_set>

namespace graph_tool
{
using namespace std;
using namespace boost;

// return generalized average nearest neighbours correlation

template <class AvgDeg>
struct get_average_nearest_neighbours_correlation
{
    get_average_nearest_neighbours_correlation(AvgDeg& avg_deg)
        : _avg_deg(avg_deg) {}

    template <class Graph, class DegreeSelectorOrigin,
              class DegreeSelectorNeighbours, class WeightMap>
    void operator()(const Graph* gp, DegreeSelectorOrigin& origin_deg,
                    DegreeSelectorNeighbours& neighbours_deg, WeightMap weight)
        const
    {
        const Graph& g = *gp;
        tr1::unordered_map<double,double> count;
        SharedMap<tr1::unordered_map<double,double> > s_count(count);
        SharedMap<AvgDeg> s_avg_deg(_avg_deg);

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) \
            firstprivate(s_count, s_avg_deg) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            typename AvgDeg::key_type orig_deg = origin_deg(v,g);

            typename graph_traits<Graph>::out_edge_iterator e, e_begin, e_end;
            tie(e_begin,e_end) = out_edges(v,g);
            for(e = e_begin; e != e_end; ++e)
            {
                typename AvgDeg::value_type::second_type::first_type deg =
                    neighbours_deg(target(*e,g),g);
                s_avg_deg[orig_deg].first += deg*get(weight, *e);
                s_avg_deg[orig_deg].second += deg*deg;
                s_count[orig_deg] += get(weight,*e);
            }
        }

        s_count.Gather();
        s_avg_deg.Gather();

        for (typeof(_avg_deg.begin()) iter = _avg_deg.begin();
             iter != _avg_deg.end(); ++iter)
        {
            double N = count[iter->first];
            iter->second.first /= N;
            if (N > 1)
                iter->second.second =
                    sqrt((iter->second.second - N*iter->second.first *
                          iter->second.first)/(N*(N-1)));
            else
                iter->second.second = 0.0;
        }
    }
    AvgDeg& _avg_deg;
};

inline void operator+=(pair<double, double>&a, const pair<double, double>&b)
{
    a.first += b.first;
    a.second += b.second;
}

} // graph_tool namespace

#endif // GRAPH_CORRELATIONS_NEIGHBOURS_HH
