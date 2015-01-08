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

#ifndef GRAPH_TRANSITION_HH
#define GRAPH_TRANSITION_HH

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_util.hh"

namespace graph_tool
{

using namespace boost;

template <class Graph, class Weight, class EdgeSelector>
typename property_traits<Weight>::value_type
sum_degree(Graph& g, typename graph_traits<Graph>::vertex_descriptor v,
           Weight w, EdgeSelector)
{
    typename property_traits<Weight>::value_type sum = 0;
    typename EdgeSelector::type e, e_end;
    for(tie(e, e_end) = EdgeSelector::get_edges(v, g); e != e_end; ++e)
        sum += get(w, *e);
    return sum;
}

template <class Graph, class EdgeSelector>
double
sum_degree(Graph& g, typename graph_traits<Graph>::vertex_descriptor v,
           ConstantPropertyMap<double, GraphInterface::edge_t> w, out_edge_iteratorS<Graph>)
{
    return out_degreeS()(v, g);
}

struct get_transition
{
    template <class Graph, class Index, class Weight>
    void operator()(const Graph& g, Index index, Weight weight,
                    multi_array_ref<double,1>& data,
                    multi_array_ref<int32_t,1>& i,
                    multi_array_ref<int32_t,1>& j) const
    {
        int pos = 0;
        for (auto v: vertices_range(g))
        {
            double k = sum_degree(g, v, weight, out_edge_iteratorS<Graph>());
            for (const auto& e: out_edges_range(v, g))
            {
                data[pos] = weight[e] / k;
                i[pos] = get(index, source(e, g));
                j[pos] = get(index, target(e, g));
                ++pos;
            }
        }
    }
};

} // namespace graph_tool

#endif // GRAPH_TRANSITION_HH
