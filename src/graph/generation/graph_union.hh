// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007-2010 Tiago de Paula Peixoto <tiago@skewed.de>
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
            vmap[*v] = add_vertex(ug);

        typename graph_traits<Graph>::edge_iterator e, e_end;
        for (tie(e,e_end) = edges(g); e != e_end; ++e)
            emap[*e] = add_edge(vmap[source(*e,g)],
                                vmap[target(*e,g)], ug).first;
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
                 is_same<typename property_traits<UnionProp>::key_type,
                         typename graph_traits<Graph>::vertex_descriptor>());
    }

    template <class UnionGraph, class Graph, class VertexMap, class EdgeMap,
              class UnionProp, class Prop>
    void dispatch(UnionGraph& ug, Graph& g, VertexMap vmap, EdgeMap emap,
                  UnionProp uprop, Prop prop, mpl::true_) const
    {
        typename graph_traits<Graph>::vertex_iterator v, v_end;
        for (tie(v,v_end) = vertices(g); v != v_end; ++v)
            uprop[vmap[*v]] = prop[*v];
    }

    template <class UnionGraph, class Graph, class VertexMap, class EdgeMap,
              class UnionProp, class Prop>
    void dispatch(UnionGraph& ug, Graph& g, VertexMap vmap, EdgeMap emap,
                  UnionProp uprop, Prop prop, mpl::false_) const
    {
        typename graph_traits<Graph>::edge_iterator e, e_end;
        for (tie(e,e_end) = edges(g); e != e_end; ++e)
            uprop[emap[*e]] = prop[*e];
    }

};

} // graph_tool namespace

#endif // GRAPH_REWIRING_HH
