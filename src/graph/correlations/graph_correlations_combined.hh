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

#ifndef GRAPH_CORRELATIONS_COMBINED_HH
#define GRAPH_CORRELATIONS_COMBINED_HH

#include <tr1/unordered_set>
#include "shared_map.hh"

namespace graph_tool
{
using namespace std;
using namespace boost;

// retrieves the distribution of combined (deg1,deg2) degrees
struct get_combined_degree_histogram
{
    template <class Graph, class DegreeSelector1, class DegreeSelector2,
              class Hist>
    void operator()(const Graph* gp, DegreeSelector1 deg1, DegreeSelector2 deg2,
                    Hist &hist) const
    {
        const Graph& g = *gp;
        SharedMap<Hist> s_hist(hist);

        typedef typename Hist::key_type::first_type first_type;
        typedef typename Hist::key_type::second_type second_type;

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) \
            firstprivate(s_hist) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            s_hist[make_pair(first_type(deg1(v,g)),
                             second_type(deg2(v,g)))]++;
        }

        s_hist.Gather();
    }
};


// retrieves the average of deg2 in function of deg1

struct get_average_combined_degree_correlation
{
    template <class Graph, class DegreeSelector1, class DegreeSelector2,
              class AvgDeg>
    void operator()(const Graph* gp, DegreeSelector1 deg1, DegreeSelector2 deg2,
                    AvgDeg& avg_deg) const
    {
        const Graph& g = *gp;
        tr1::unordered_map<double,size_t> count;

        typename graph_traits<Graph>::vertex_iterator v, v_begin, v_end;
        tie(v_begin, v_end) = vertices(g);
        for(v = v_begin; v != v_end; ++v)
        {
            typename AvgDeg::key_type d1 = deg1(*v,g);
            typename AvgDeg::value_type::second_type::first_type d2 =
                deg2(*v, g);
            avg_deg[d1].first += d2;
            avg_deg[d1].second += d2*d2;
            count[d1]++;
        }

        for (typeof(avg_deg.begin()) iter = avg_deg.begin();
             iter != avg_deg.end(); ++iter)
        {
            size_t N = count[iter->first];
            iter->second.first /= N;
            if (N > 1)
            {
                double err = (iter->second.second - N*iter->second.first *
                              iter->second.first)/(N*(N-1));
                iter->second.second = (err<0.0)?0.0:sqrt(err);
            }
            else
            {
                iter->second.second = 0.0;
            }
        }
    }
};

} // graph_tool namespace

#endif // GRAPH_CORRELATIONS_COMBINED_HH

