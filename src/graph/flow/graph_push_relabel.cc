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

#include "graph_filtering.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"
#include "graph.hh"

#include "graph_augment.hh"

#include <boost/mpl/if.hpp>
#include <boost/mpl/or.hpp>
namespace std
{
using namespace boost;

// we need a min() function with arguments of different types

template <class T1, class T2>
typename boost::mpl::if_<
    typename boost::mpl::or_<typename std::is_floating_point<T1>::type,
                             typename std::is_floating_point<T2>::type>::type,
    double, int>::type
min(const T1& v1, const T2& v2)
{
    if (v1 <= T1(v2))
        return v1;
    else
        return v2;
}
}

#include <boost/graph/push_relabel_max_flow.hpp>

#include <boost/bind.hpp>

using namespace graph_tool;
using namespace boost;



struct get_push_relabel_max_flow
{
    template <class Graph, class VertexIndex, class EdgeIndex, class CapacityMap,
              class ResidualMap>
    void operator()(Graph& g, VertexIndex vertex_index, EdgeIndex edge_index,
                    size_t max_e, size_t src, size_t sink, CapacityMap cm,
                    ResidualMap res) const
    {
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;
        checked_vector_property_map<bool,EdgeIndex>
            augmented(edge_index);
        unchecked_vector_property_map<edge_t,EdgeIndex>
            reverse_map(edge_index, max_e);

        typedef typename std::remove_const<Graph>::type GT;
        GT& u = const_cast<GT&>(g);

        augment_graph(u, augmented, cm,
                      reverse_map.get_checked(), res);

        boost::push_relabel_max_flow(g, vertex(src, g), vertex(sink, g),
                                     get_unchecked(cm),
                                     res.get_unchecked(),
                                     reverse_map, vertex_index);
        deaugment_graph(u, augmented);
    }
};


void push_relabel_max_flow(GraphInterface& gi, size_t src, size_t sink,
                           boost::any capacity, boost::any res)
{
    run_action<graph_tool::detail::always_directed, boost::mpl::true_>()
        (gi, std::bind(get_push_relabel_max_flow(),
                       placeholders::_1, gi.GetVertexIndex(), gi.GetEdgeIndex(),
                       gi.GetMaxEdgeIndex(),
                       src, sink,  placeholders::_2,  placeholders::_3),
         writable_edge_scalar_properties(), writable_edge_scalar_properties())
        (capacity,res);
}
