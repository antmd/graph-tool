// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2013 Tiago de Paula Peixoto <tiago@skewed.de>
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

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <iostream>

using namespace std;
using namespace boost;
using namespace graph_tool;


// this is the constructor for the graph interface
GraphInterface::GraphInterface()
    :_state(new state_t()),
     _vertex_index(get(vertex_index, _state->_mg)),
     _edge_index(get(edge_index_t(), _state->_mg)),
     _reversed(false),
     _directed(true),
     _graph_index(0),
     _vertex_filter_map(_vertex_index),
     _vertex_filter_invert(false),
     _vertex_filter_active(false),
     _edge_filter_map(_edge_index),
     _edge_filter_invert(false),
     _edge_filter_active(false)
{
    _state->_nedges = _state->_max_edge_index = 0;
}

// the destructor
GraphInterface::~GraphInterface()
{
}

// this will get the number of vertices, either the "soft" O(1) way, or the hard
// O(V) way, which is necessary if the graph is filtered
size_t GraphInterface::GetNumberOfVertices()
{
    size_t n = 0;
    if (IsVertexFilterActive())
        run_action<>()(*this, lambda::var(n) =
                       lambda::bind<size_t>(HardNumVertices(),lambda::_1))();
    else
        run_action<>()(*this, lambda::var(n) =
                       lambda::bind<size_t>(SoftNumVertices(),lambda::_1))();
    return n;
}

// this will get the number of edges, either the "soft" O(E) way, or the hard
// O(E) way, which is necessary if the graph is filtered. Both cases are of
// linear complexity, since num_edges() is O(E) in Boost's adjacency_list
size_t GraphInterface::GetNumberOfEdges()
{
    using namespace boost::lambda;
    size_t n = 0;
    if (IsEdgeFilterActive() || IsVertexFilterActive())
        run_action<>()(*this, lambda::var(n) =
                       lambda::bind<size_t>(HardNumEdges(),lambda::_1))();
    else
        n = _state->_nedges;
    return n;
}

void GraphInterface::Clear()
{
    _state->_mg.clear();
    _state->_nedges = 0;
    _state->_max_edge_index = 0;
    _state->_free_indexes.clear();
}

void GraphInterface::ClearEdges()
{
    graph_traits<multigraph_t>::vertex_iterator v, v_end;
    for (tie(v, v_end) = vertices(_state->_mg); v != v_end; ++v)
        clear_vertex(*v, _state->_mg);
    _state->_nedges = 0;
    _state->_max_edge_index = 0;
    _state->_free_indexes.clear();
}

