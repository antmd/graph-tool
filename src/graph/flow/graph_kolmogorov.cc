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

#include <boost/bind.hpp>

using namespace graph_tool;
using namespace boost;

#if (BOOST_VERSION >= 104400)
# include <boost/graph/boykov_kolmogorov_max_flow.hpp>
# define KOLMOGOROV_MAX_FLOW boost::boykov_kolmogorov_max_flow
#else
# include <boost/graph/kolmogorov_max_flow.hpp>
# define KOLMOGOROV_MAX_FLOW boost::kolmogorov_max_flow
#endif

struct get_kolmogorov_max_flow
{
    template <class Graph, class EdgeIndex, class VertexIndex,
              class CapacityMap, class ResidualMap>
    void operator()(Graph& g, EdgeIndex edge_index, size_t,
                    VertexIndex vertex_index, size_t src, size_t sink,
                    CapacityMap cm, ResidualMap res) const
    {
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;
        checked_vector_property_map<uint8_t,EdgeIndex>
            augmented(edge_index);
        checked_vector_property_map<edge_t,EdgeIndex>
            reverse_map(edge_index);
        unchecked_vector_property_map<edge_t,VertexIndex>
            pred_map(vertex_index, num_vertices(g));
        unchecked_vector_property_map<size_t,VertexIndex>
            color_map(vertex_index, num_vertices(g));
        unchecked_vector_property_map<size_t,VertexIndex>
            dist_map(vertex_index, num_vertices(g));

        typedef typename std::remove_const<Graph>::type GT;
        GT& u = const_cast<GT&>(g);

        augment_graph(u, augmented, cm, reverse_map, res, true);

        KOLMOGOROV_MAX_FLOW(g, cm, res, reverse_map, pred_map, color_map,
                            dist_map, vertex_index, vertex(src, g),
                            vertex(sink, g));

        deaugment_graph(u, augmented);

    }
};


void kolmogorov_max_flow(GraphInterface& gi, size_t src, size_t sink,
                         boost::any capacity, boost::any res)
{
    run_action<graph_tool::detail::always_directed, boost::mpl::true_>()
        (gi, std::bind(get_kolmogorov_max_flow(),
                       placeholders::_1, gi.GetEdgeIndex(), gi.GetMaxEdgeIndex(),
                       gi.GetVertexIndex(), src, sink,  placeholders::_2,
                       placeholders::_3),
         writable_edge_scalar_properties(), writable_edge_scalar_properties())
        (capacity,res);
}
