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

#ifndef GRAPH_AUGMENT_HH
#define GRAPH_AUGMENT_HH

namespace graph_tool
{
using namespace std;
using namespace boost;

template <class Graph, class AugmentedMap, class CapacityMap,
          class ReversedMap,  class ResidualMap>
void augment_graph(Graph& g, AugmentedMap augmented, CapacityMap capacity,
                   ReversedMap rmap, ResidualMap res,
                   bool detect_reversed = false)
{
    for (auto e : edges_range(g))
        augmented[e] = 0;

    vector<typename graph_traits<Graph>::edge_descriptor> e_list;
    for (auto v : vertices_range(g))
    {
        e_list.clear();

        for (auto e : out_edges_range(v, g))
        {
            if (detect_reversed && augmented[e] == 0)
            {
                for (auto e2 : out_edges_range(target(e, g), g))
                {
                    if (augmented[e2] != 0)
                        continue;

                    if (target(e2, g) == v)
                    {
                        augmented[e] = 2;
                        augmented[e2] = 2;
                        rmap[e] = e2;
                        rmap[e2] = e;
                        break;
                    }
                }
            }

            if (augmented[e] == 0)
                e_list.push_back(e);
        }

        for (auto& e : e_list)
        {
            auto ae = add_edge(target(e, g), source(e, g), g).first;
            augmented[ae] = 1;
            capacity[ae] = 0;
            rmap[e] = ae;
            rmap[ae] = e;
            res[ae] = 0;
        }
    }
}

template <class Graph, class AugmentedMap>
void deaugment_graph(Graph& g, AugmentedMap augmented)
{
    vector<typename graph_traits<Graph>::edge_descriptor> e_list;
    for (auto v : vertices_range(g))
    {
        e_list.clear();
        for (auto e : out_edges_range(v, g))
        {
            if (augmented[e] == 1)
                e_list.push_back(e);
        }

        for (auto& e : e_list)
            remove_edge(e, g);
    }
}


template <class Graph, class CapacityMap, class ResidualMap,
          class AugmentedMap>
void residual_graph(Graph& g, CapacityMap capacity, ResidualMap res,
                    AugmentedMap augmented)
{
    vector<typename graph_traits<Graph>::edge_descriptor> e_list;
    for (auto e : edges_range(g))
    {
        if (capacity[e] - res[e] > 0)
            e_list.push_back(e);
    }

    for (auto& e : e_list)
    {
        auto ne = add_edge(target(e, g), source(e, g), g);
        augmented[ne.first] = true;
    }
}

} // graph_tool namespace

#endif // GRAPH_AUGMENT_HH
