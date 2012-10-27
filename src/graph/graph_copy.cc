// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007-2012 Tiago de Paula Peixoto <tiago@skewed.de>
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

#include <boost/mpl/contains.hpp>
#include <boost/python/extract.hpp>

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_properties.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

struct do_graph_copy
{
    template <class GraphDst, class GraphSrc, class DstVertexIndexMap,
              class SrcVertexIndexMap,  class DstEdgeIndexMap,
              class SrcEdgeIndexMap>
    void operator()(const GraphSrc& src, GraphDst& dst,
                    DstVertexIndexMap src_vertex_index,
                    SrcVertexIndexMap dst_vertex_index,
                    DstEdgeIndexMap src_edge_index,
                    SrcEdgeIndexMap dst_edge_index) const
    {
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
            dst_edge_index[new_e] = src_edge_index[*e];
        }
    }
};

// copy constructor
GraphInterface::GraphInterface(const GraphInterface& gi, bool keep_ref)
    :_state(keep_ref ? gi._state : shared_ptr<state_t>(new state_t())),
     _vertex_index(get(vertex_index, _state->_mg)),
     _edge_index(get(edge_index_t(), _state->_mg)),
     _reversed(gi._reversed),
     _directed(gi._directed),
     _vertex_filter_map(_vertex_index),
     _vertex_filter_invert(false),
     _vertex_filter_active(false),
     _edge_filter_map(_edge_index),
     _edge_filter_invert(false),
     _edge_filter_active(false)
{
    if (keep_ref)
        return;
    _state->_nedges = gi._state->_nedges;
    _state->_max_edge_index = gi._state->_max_edge_index;
    _state->_free_indexes = gi._state->_free_indexes;

    run_action<>()
        (const_cast<GraphInterface&>(gi),
         bind<void>(do_graph_copy(), _1, ref(_state->_mg),
                    gi._vertex_index, _vertex_index, 
                    gi._edge_index, _edge_index))();
    // filters will be copied in python
}

