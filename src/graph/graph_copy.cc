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

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_properties.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

struct graph_copy
{
    template <class GraphDst, class GraphSrc, class DstVertexIndexMap,
              class SrcVertexIndexMap,  class DstEdgeIndexMap,
              class SrcEdgeIndexMap>
    void operator()(GraphDst* dstp, GraphSrc* srcp,
                    DstVertexIndexMap dst_vertex_index,
                    SrcVertexIndexMap src_vertex_index,
                    DstEdgeIndexMap dst_edge_index,
                    SrcEdgeIndexMap src_edge_index) const
    {
        GraphDst& dst = *dstp;
        GraphSrc& src = *srcp;

        vector<size_t> index_map(num_vertices(src));
        typename graph_traits<GraphSrc>::vertex_iterator v, v_end;
        for (tie(v, v_end) = vertices(src); v != v_end; ++v)
        {
            if (src_vertex_index[*v] >= index_map.size())
                index_map.resize(src_vertex_index[*v]+1);
            typename graph_traits<GraphDst>::vertex_descriptor new_v =
                add_vertex(dst);
            index_map[src_vertex_index[*v]] = dst_vertex_index[new_v];
        }

        typename graph_traits<GraphSrc>::edge_iterator e, e_end;
        for (tie(e, e_end) = edges(src); e != e_end; ++e)
        {
            size_t s = index_map[src_vertex_index[source(*e, src)]];
            size_t t = index_map[src_vertex_index[target(*e, src)]];
            typename graph_traits<GraphDst>::edge_descriptor new_e =
                add_edge(vertex(s,dst), vertex(t,dst), dst).first;
            dst_edge_index[new_e] = src_edge_index[new_e];
        }
    }
};

struct copy_vertex_property
{
    template <class GraphDst, class GraphSrc, class VertexProperty>
    void operator()(GraphDst* dstp, GraphSrc* srcp,
                    VertexProperty dst_map,
                    boost::any src_property) const
    {
        GraphDst& dst = *dstp;
        GraphSrc& src = *srcp;
        VertexProperty src_map = any_cast<VertexProperty>(src_property);
        typename graph_traits<GraphSrc>::vertex_iterator v, v_end;
        size_t count = 0;
        for (tie(v, v_end) = vertices(src); v != v_end; ++v)
        {
            dst_map[vertex(count,dst)] = src_map[*v];
            ++count;
        }
    }
};

struct copy_edge_property
{
    template <class GraphSrc, class EdgeProperty>
    void operator()(GraphSrc* srcp,
                    EdgeProperty dst_map,
                    boost::any src_property) const
    {
        GraphSrc& src = *srcp;
        EdgeProperty src_map = any_cast<EdgeProperty>(src_property);
        typename graph_traits<GraphSrc>::edge_iterator e, e_end;
        for (tie(e, e_end) = edges(src); e != e_end; ++e)
        {
            put(dst_map, *e, get(src_map, *e));
        }
    }
};


// copy constructor
GraphInterface::GraphInterface(const GraphInterface& gi)
    :_mg(),
     _reversed(false),
     _directed(gi._directed),
     _vertex_index(get(vertex_index,_mg)),
     _edge_index(get(edge_index_t(),_mg)),
     _vertex_filter_map(_vertex_index),
     _vertex_filter_invert(false),
     _vertex_filter_active(false),
     _edge_filter_map(_edge_index),
     _edge_filter_invert(false),
     _edge_filter_active(false)
{
    typedef mpl::transform<detail::all_graph_views,
                           mpl::quote1<add_pointer> >::type graph_views;

    run_action<detail::never_filtered>()(*this, bind<void>(graph_copy(), _1, _2,
                                                           _vertex_index,
                                                           gi._vertex_index,
                                                           _edge_index,
                                                           gi._edge_index),
                                         graph_views())(gi.GetGraphView());

    typedef property_map_types::apply<value_types, vertex_index_map_t,
                                      mpl::bool_<false> >::type
        vertex_properties;

    typedef property_map_types::apply<value_types, edge_index_map_t,
                                      mpl::bool_<false> >::type
        edge_properties;

    for (typeof(gi._properties.begin()) iter = gi._properties.begin();
         iter != gi._properties.end(); ++iter)
    {
        string name = iter->first;
        string type = get_type_name<>()(iter->second->value());
        if (iter->second->key() == typeid(vertex_t))
        {
            AddVertexProperty(name, type);
            run_action<detail::never_filtered>()
                (*this, bind<void>(copy_vertex_property(), _1, _2, _3,
                                   prop(name, gi._vertex_index,
                                        gi._properties)),
                 graph_views(), vertex_properties())
                (gi.GetGraphView(), prop(name, _vertex_index, _properties));
        }
        else if (iter->second->key() == typeid(edge_t))
        {
            AddEdgeProperty(name, type);
            run_action<detail::never_filtered>()
                (*this, bind<void>(copy_edge_property(), _1, _2,
                                   prop(name, gi._edge_index, gi._properties)),
                 edge_properties())
                (prop(name, _edge_index, _properties));
        }
        else
        {
            AddGraphProperty(name, type);
            put(name, _properties, graph_property_tag(),
                iter->second->get(graph_property_tag()));
        }
    }
}
