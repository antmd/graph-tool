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

#include <boost/lambda/bind.hpp>

#include "graph.hh"
#include "graph_filtering.hh"

#include <boost/graph/isomorphism.hpp>

using namespace graph_tool;
using namespace boost;
using namespace boost::lambda;

struct check_iso
{
    template <class Graph1, class Graph2, class IsoMap, class VertexIndexMap>
    void operator()(Graph1& g1, Graph2* g2, IsoMap map, VertexIndexMap index1,
                    VertexIndexMap index2, bool& result) const
    {
        result = isomorphism(g1, *g2, isomorphism_map(map).
                             vertex_index1_map(index1).
                             vertex_index2_map(index2));
    }
};

struct directed_graph_view_pointers:
    mpl::transform<graph_tool::detail::always_directed,
                   mpl::quote1<add_pointer> >::type {};

struct undirected_graph_view_pointers:
    mpl::transform<graph_tool::detail::never_directed,
                   mpl::quote1<add_pointer> >::type {};

typedef property_map_types::apply<integer_types,
                                  GraphInterface::vertex_index_map_t,
                                  mpl::bool_<false> >::type
    vertex_props_t;

bool check_isomorphism(GraphInterface& gi1, GraphInterface& gi2,
                       boost::any iso_map)
{
    bool result;

    if (gi1.GetDirected() != gi2.GetDirected())
        return false;
    if (gi1.GetDirected())
    {
        run_action<graph_tool::detail::always_directed>()
            (gi1, bind<void>(check_iso(),
                             _1, _2, _3, gi1.GetVertexIndex(),
                             gi2.GetVertexIndex(), var(result)),
             directed_graph_view_pointers(), vertex_props_t())
            (gi2.GetGraphView(), iso_map);
    }
    else
    {
        run_action<graph_tool::detail::never_directed>()
            (gi1, bind<void>(check_iso(),
                             _1, _2, _3, gi1.GetVertexIndex(),
                             gi2.GetVertexIndex(), var(result)),
             undirected_graph_view_pointers(), vertex_props_t())
            (gi2.GetGraphView(), iso_map);
    }

    return result;
}
