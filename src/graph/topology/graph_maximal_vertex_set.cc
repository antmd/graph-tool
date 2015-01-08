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

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_properties.hh"
#include "graph_selectors.hh"
#include "graph_util.hh"

#include "random.hh"

#include <boost/python.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

struct do_maximal_vertex_set
{
    template <class Graph, class VertexIndex, class VertexSet,
              class RNG>
    void operator()(const Graph& g, VertexIndex vertex_index, VertexSet mvs,
                    bool high_deg, RNG& rng) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;
        typedef typename property_traits<VertexSet>::value_type wval_t;

        uniform_real_distribution<> sample(0, 1);

        VertexSet marked(vertex_index, num_vertices(g));
        vector<vertex_t> vlist;
        double max_deg = 0, tmp_max_deg = 0;
        typename graph_traits<Graph>::vertex_iterator v, v_end;
        for (tie(v, v_end) = vertices(g); v != v_end; ++v)
        {
            vlist.push_back(*v);
            mvs[*v] = marked[*v] = false;
            max_deg = max(out_degree(*v, g), max_deg);
        }

        vector<vertex_t> selected, tmp;
        tmp.reserve(vlist.size());
        selected.reserve(vlist.size());
        while (!vlist.empty())
        {
            selected.clear();
            tmp.clear();
            tmp_max_deg = 0;

            int i, N = vlist.size();
            #pragma omp parallel for default(shared) private(i)
            for (i = 0; i < N; ++i)
            {
                typename graph_traits<Graph>::vertex_descriptor v =
                    vlist[i];

                marked[v] = false;
                bool include = true;
                typename graph_traits<Graph>::adjacency_iterator a, a_end;
                for(tie(a, a_end) = adjacent_vertices(v, g); a != a_end; ++a)
                {
                    if (mvs[*a])
                    {
                        include = false;
                        break;
                    }
                }
                if (!include)
                    continue;

                include = false;
                if (out_degree(v, g) > 0)
                {
                    double p, r;
                    if (high_deg)
                        p = out_degree(v, g) / max_deg;
                    else
                        p = 1. / (2 * out_degree(v, g));

                    #pragma omp critical
                    {
                        r = sample(rng);
                    }
                    if (r < p)
                        include = true;
                }
                else
                {
                    include = true;
                }

                if (include)
                {
                    marked[v] = true;
                    #pragma omp critical
                    {
                        selected.push_back(v);
                    }
                }
                else
                {
                    #pragma omp critical
                    {
                        tmp.push_back(v);
                        tmp_max_deg = max(tmp_max_deg, out_degree(v, g));
                    }
                }
            }

            N = selected.size();
            #pragma omp parallel for default(shared) private(i)
            for (i = 0; i < N; ++i)
            {
                typename graph_traits<Graph>::vertex_descriptor v =
                    selected[i];
                bool include = true;
                typename graph_traits<Graph>::adjacency_iterator a, a_end;
                for(tie(a, a_end) = adjacent_vertices(v, g); a != a_end; ++a)
                {
                    if (*a == v)  //skip self-loops
                        continue;
                    if (mvs[*a])
                    {
                        include = false;
                        break;
                    }

                    if (marked[*a])
                    {
                        bool inc = ((high_deg && (out_degree(v, g) >
                                                  out_degree(*a, g))) ||
                                    (!high_deg && (out_degree(v, g) <
                                                   out_degree(*a, g))));
                        if (out_degree(v, g) == out_degree(*a, g))
                            inc = v < *a;
                        include = include && inc;
                    }
                }

                if (include)
                {
                    mvs[v] = true;
                }
                else
                {
                    #pragma omp critical
                    tmp.push_back(v);
                    tmp_max_deg = max(tmp_max_deg, out_degree(v, g));
                }
                marked[v] = false;
            }

            vlist = tmp;
            max_deg = tmp_max_deg;
        }
    }
};

void maximal_vertex_set(GraphInterface& gi, boost::any mvs, bool high_deg,
                        rng_t& rng)
{
    run_action<>()
        (gi, std::bind(do_maximal_vertex_set(), placeholders::_1, gi.GetVertexIndex(),
                       placeholders::_2, high_deg, std::ref(rng)),
         writable_vertex_scalar_properties())(mvs);
}

void export_maximal_vertex_set()
{
    python::def("maximal_vertex_set", &maximal_vertex_set);
}
