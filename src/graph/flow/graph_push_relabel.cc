// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007-2010 Tiago de Paula Peixoto <tiago@forked.de>
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
typename mpl::if_<
    typename mpl::or_<typename is_floating_point<T1>::type,
                      typename is_floating_point<T2>::type>::type,
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
    template <class Graph, class EdgeIndex, class CapacityMap,
              class ResidualMap>
    void operator()(Graph& g, EdgeIndex edge_index, size_t max_e, size_t src,
                    size_t sink, CapacityMap cm, ResidualMap res) const
    {
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;
        unchecked_vector_property_map<bool,EdgeIndex>
            augmented(edge_index, max_e);
        unchecked_vector_property_map<edge_t,EdgeIndex>
            reverse_map(edge_index, max_e);
        augment_graph(g, augmented.get_checked(), get_checked(cm),
                      reverse_map.get_checked(), res.get_checked());
        boost::push_relabel_max_flow(g._g, vertex(src, g), vertex(sink, g),
                                     capacity_map(cm).
                                     reverse_edge_map(reverse_map).
                                     residual_capacity_map(res));
        deaugment_graph(g, augmented.get_checked());
    }
};


void push_relabel_max_flow(GraphInterface& gi, size_t src, size_t sink,
                           boost::any capacity, boost::any res)
{
    run_action<graph_tool::detail::always_directed, mpl::true_>()
        (gi, bind<void>(get_push_relabel_max_flow(),
                        _1, gi.GetEdgeIndex(),  gi.GetMaxEdgeIndex(),
                        src, sink, _2, _3),
         edge_scalar_properties(), writable_edge_scalar_properties())
        (capacity,res);
}
