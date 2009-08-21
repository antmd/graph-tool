// Copyright (C) 2008  Tiago de Paula Peixoto <tiago@forked.de>
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
                   ReversedMap rmap, ResidualMap res)
{
    typename graph_traits<Graph>::edge_iterator e, e_end;
    for (tie(e,e_end) = edges(g); e != e_end; ++e)
        augmented[*e] = false;

    typename graph_traits<Graph>::vertex_iterator v, v_end;
    vector<typename graph_traits<Graph>::edge_descriptor> e_list;
    for (tie(v,v_end) = vertices(g); v != v_end; ++v)
    {
        e_list.clear();

        typename graph_traits<Graph>::out_edge_iterator e, e_end;
        for (tie(e, e_end) = out_edges(*v, g); e != e_end; ++e)
        {
            if (!augmented[*e])
                e_list.push_back(*e);
        }

        for (size_t i = 0; i < e_list.size(); ++i)
        {
            typename graph_traits<Graph>::edge_descriptor ae;
            ae = add_edge(target(e_list[i],g), source(e_list[i],g), g).first;
            augmented[ae] = true;
            capacity[ae] = 0;
            rmap[e_list[i]] = ae;
            rmap[ae] = e_list[i];
            res[ae] = 0;
        }
    }
}

template <class Graph, class AugmentedMap>
void deaugment_graph(Graph& g, AugmentedMap augmented)
{
    vector<typename graph_traits<Graph>::edge_descriptor> e_list;
    typename graph_traits<Graph>::vertex_iterator v, v_end;
    for (tie(v,v_end) = vertices(g); v != v_end; ++v)
    {
        e_list.clear();
        typename graph_traits<Graph>::out_edge_iterator e, e_end;
        for (tie(e, e_end) = out_edges(*v, g); e != e_end; ++e)
        {
            if (augmented[*e])
                e_list.push_back(*e);
        }

        for (size_t i = 0; i < e_list.size(); ++i)
            remove_edge(e_list[i],g);
    }
}


template <class Type, class Index>
checked_vector_property_map<Type, Index>
get_checked(unchecked_vector_property_map<Type, Index> prop)
{
    return prop.get_checked();
}

template <class Prop>
Prop
get_checked(Prop prop)
{
    return prop;
}

} // graph_tool namespace

#endif // GRAPH_AUGMENT_HH
