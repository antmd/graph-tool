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
#include "graph_properties.hh"

#include <boost/mpl/contains.hpp>
#include <boost/python/extract.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

//
// Property map copying
// ====================

template <class IteratorSel, class PropertyMaps>
struct copy_property
{
    template <class GraphTgt, class GraphSrc, class PropertyTgt>
    void operator()(const GraphTgt& tgt, const GraphSrc* src,
                    PropertyTgt dst_map, boost::any prop_src) const
    {
        try
        {
            typedef typename property_traits<PropertyTgt>::value_type val_tgt;
            typedef typename IteratorSel::template get_descriptor<GraphSrc>::type src_d;

            DynamicPropertyMapWrap<val_tgt, src_d> src_map(prop_src, PropertyMaps());

            typename IteratorSel::template apply<GraphSrc>::type vs, vs_end;
            typename IteratorSel::template apply<GraphTgt>::type vt, vt_end;
            tie(vt, vt_end) = IteratorSel::range(tgt);
            for (tie(vs, vs_end) = IteratorSel::range(*src); vs != vs_end; ++vs)
            {
                if (vt == vt_end)
                    throw ValueException("Error copying properties: "
                                         "graphs not compatible");
                dst_map[*vt] = get(src_map, *vs);
                ++vt;
            }
        }
        catch (bad_lexical_cast&)
        {
            throw ValueException("property values are not convertible");
        }
    }
};

struct edge_selector
{
    template <class Graph>
    struct apply
    {
        typedef typename graph_traits<Graph>::edge_iterator type;
    };

    template <class Graph>
    struct get_descriptor
    {
        typedef typename graph_traits<Graph>::edge_descriptor type;
    };

    template <class Graph>
    static pair<typename apply<Graph>::type,
                typename apply<Graph>::type>
    range(Graph& g)
    {
        return edges(g);
    }
};

struct vertex_selector
{
    template <class Graph>
    struct apply
    {
        typedef typename graph_traits<Graph>::vertex_iterator type;
    };

    template <class Graph>
    struct get_descriptor
    {
        typedef typename graph_traits<Graph>::vertex_descriptor type;
    };

    template <class Graph>
    static pair<typename apply<Graph>::type,
                typename apply<Graph>::type>
    range(Graph& g)
    {
        return vertices(g);
    }
};

struct graph_views:
    boost::mpl::transform<graph_tool::detail::always_directed_never_reversed,
                          boost::mpl::quote1<std::add_pointer> >::type {};

void GraphInterface::CopyVertexProperty(const GraphInterface& src,
                                        boost::any prop_src,
                                        boost::any prop_tgt)
{
    run_action<>()
        (*this, std::bind(copy_property<vertex_selector,vertex_properties>(),
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3, prop_src),
         graph_views(), writable_vertex_properties())
        (src.GetGraphView(), prop_tgt);
}

void GraphInterface::CopyEdgeProperty(const GraphInterface& src,
                                      boost::any prop_src,
                                      boost::any prop_tgt)
{
    run_action<>()
        (*this, std::bind(copy_property<edge_selector,edge_properties>(),
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3, prop_src),
         graph_views(), writable_edge_properties())
        (src.GetGraphView(), prop_tgt);
}
