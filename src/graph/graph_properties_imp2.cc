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

#include "graph_python_interface.hh"
#include "graph.hh"
#include "graph_properties.hh"
#include "graph_filtering.hh"
#include "graph_selectors.hh"
#include "graph_util.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

template <class Val1, class Val2>
void operator+=(std::vector<Val1>& v1, const std::vector<Val2>& v2)
{
    if (v2.size() > v1.size())
        v1.resize(v2.size());
    for (size_t i = 0; i < v2.size(); ++i)
        v1[i] += v2[i];
}

struct do_out_edges_sum
{
    template <class Graph, class EProp, class VProp>
    void operator()(Graph& g, EProp eprop, VProp vprop) const
    {

        typedef typename property_traits<EProp>::value_type eval_t;
        typedef typename property_traits<VProp>::value_type vval_t;

        convert<vval_t, eval_t> conv;

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i)     \
            schedule(runtime) if (N > 100)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            vprop[v] = vval_t();
            for (auto e : out_edges_range(v, g))
                vprop[v] += conv(eprop[e]);
        }
    }
};

void out_edges_sum(GraphInterface& gi, boost::any eprop, boost::any vprop)
{
    run_action<>()(gi, std::bind(do_out_edges_sum(), placeholders::_1,
                                 placeholders::_2, placeholders::_3),
                   edge_properties(), writable_vertex_properties())
        (eprop, vprop);
}
