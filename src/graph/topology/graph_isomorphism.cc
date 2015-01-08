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

#include <boost/graph/isomorphism.hpp>

using namespace graph_tool;
using namespace boost;

struct check_iso
{

    template <class Graph1, class Graph2, class IsoMap, class InvMap,
              class VertexIndexMap>
    void operator()(Graph1& g1, Graph2* g2, InvMap cinv_map1, InvMap cinv_map2,
                    int64_t max_inv, IsoMap map, VertexIndexMap index1,
                    VertexIndexMap index2, bool& result) const
    {
        auto inv_map1 = cinv_map1.get_unchecked(num_vertices(g1));
        auto inv_map2 = cinv_map2.get_unchecked(num_vertices(*g2));

        vinv_t<decltype(inv_map1)> vinv1(inv_map1, max_inv);
        vinv_t<decltype(inv_map2)> vinv2(inv_map2, max_inv);

        result = isomorphism(g1, *g2,
                             isomorphism_map(map.get_unchecked(num_vertices(g1))).
                             vertex_invariant1(vinv1).
                             vertex_invariant2(vinv2).
                             vertex_index1_map(index1).
                             vertex_index2_map(index2));
    }

    template <class Prop>
    struct vinv_t
    {
        vinv_t(Prop& prop, int64_t max)
            : _prop(prop), _max(max) {}
        Prop& _prop;
        int64_t _max;

        template <class Vertex>
        int64_t operator()(Vertex v) const
        {
            return _prop[v];
        };

        int64_t max() const { return _max; }

        typedef int64_t result_type;
        typedef size_t argument_type;
    };
};

struct directed_graph_view_pointers:
    mpl::transform<graph_tool::detail::always_directed,
                   mpl::quote1<std::add_pointer> >::type {};

struct undirected_graph_view_pointers:
    mpl::transform<graph_tool::detail::never_directed,
                   mpl::quote1<std::add_pointer> >::type {};

typedef property_map_types::apply<integer_types,
                                  GraphInterface::vertex_index_map_t,
                                  mpl::bool_<false> >::type
    vertex_props_t;

bool check_isomorphism(GraphInterface& gi1, GraphInterface& gi2,
                       boost::any ainv_map1, boost::any ainv_map2,
                       int64_t max_inv, boost::any aiso_map)
{
    bool result;

    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        iso_map_t;
    auto iso_map = any_cast<iso_map_t>(aiso_map);

    typedef property_map_type::apply<int64_t,
                                     GraphInterface::vertex_index_map_t>::type
        inv_map_t;
    auto inv_map1 = any_cast<inv_map_t>(ainv_map1);
    auto inv_map2 = any_cast<inv_map_t>(ainv_map2);

    if (gi1.GetDirected() != gi2.GetDirected())
        return false;
    if (gi1.GetDirected())
    {
        run_action<graph_tool::detail::always_directed>()
            (gi1, std::bind(check_iso(),
                            placeholders::_1, placeholders::_2,
                            inv_map1, inv_map2, max_inv, iso_map,
                            gi1.GetVertexIndex(),
                            gi2.GetVertexIndex(), std::ref(result)),
             directed_graph_view_pointers())
            (gi2.GetGraphView());
    }
    else
    {
        run_action<graph_tool::detail::never_directed>()
            (gi1, std::bind(check_iso(),
                            placeholders::_1, placeholders::_2,
                            inv_map1, inv_map2, max_inv, iso_map,
                            gi1.GetVertexIndex(),
                            gi2.GetVertexIndex(), std::ref(result)),
             undirected_graph_view_pointers())
            (gi2.GetGraphView());
    }

    return result;
}
