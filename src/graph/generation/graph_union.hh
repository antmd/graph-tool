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

#ifndef GRAPH_REWIRING_HH
#define GRAPH_REWIRING_HH

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_util.hh"

namespace graph_tool
{
using namespace std;
using namespace boost;

struct graph_union
{
    template <class UnionGraph, class Graph, class VertexMap, class EdgeMap>
    void operator()(UnionGraph& ug, Graph* gp, VertexMap vmap, EdgeMap emap)
        const
    {
        Graph& g = *gp;
        typename graph_traits<Graph>::vertex_iterator v, v_end;
        for (tie(v,v_end) = vertices(g); v != v_end; ++v)
        {
            if (vmap[*v] == 0)
            {
                vmap[*v] = add_vertex(ug);
            }
            else
            {
                typename graph_traits<UnionGraph>::vertex_descriptor w =
                    vertex(vmap[*v] - 1, ug);
                if (w == graph_traits<UnionGraph>::null_vertex() ||
                    w >= num_vertices(g))
                    vmap[*v] = add_vertex(ug);
                else
                    vmap[*v] = w;
            }
        }

        typename graph_traits<Graph>::edge_iterator e, e_end;
        for (tie(e,e_end) = edges(g); e != e_end; ++e)
            emap[*e] = add_edge(vertex(vmap[source(*e,g)], g),
                                vertex(vmap[target(*e,g)], g), ug).first;
    }
};


struct property_union
{
    template <class UnionGraph, class Graph, class VertexMap, class EdgeMap,
              class UnionProp>
    void operator()(UnionGraph& ug, Graph* gp, VertexMap vmap, EdgeMap emap,
                    UnionProp uprop, boost::any aprop) const
    {
        Graph& g = *gp;
        typename UnionProp::checked_t prop =
            any_cast<typename UnionProp::checked_t>(aprop);
        dispatch(ug, g, vmap, emap, uprop, prop,
                 std::is_same<typename property_traits<UnionProp>::key_type,
                              typename graph_traits<Graph>::vertex_descriptor>());
    }

    template <class UnionGraph, class Graph, class VertexMap, class EdgeMap,
              class UnionProp, class Prop>
    void dispatch(UnionGraph&, Graph& g, VertexMap vmap, EdgeMap,
                  UnionProp uprop, Prop prop, std::true_type) const
    {
        typename graph_traits<Graph>::vertex_iterator v, v_end;
        for (tie(v,v_end) = vertices(g); v != v_end; ++v)
            uprop[vertex(vmap[*v], g)] = prop[*v];
    }

    template <class UnionGraph, class Graph, class VertexMap, class EdgeMap,
              class UnionProp, class Prop>
    void dispatch(UnionGraph&, Graph& g, VertexMap, EdgeMap emap,
                  UnionProp uprop, Prop prop, std::false_type) const
    {
        typename graph_traits<Graph>::edge_iterator e, e_end;
        for (tie(e,e_end) = edges(g); e != e_end; ++e)
            uprop[emap[*e]] = prop[*e];
    }

};

} // graph_tool namespace

#endif // GRAPH_REWIRING_HH
