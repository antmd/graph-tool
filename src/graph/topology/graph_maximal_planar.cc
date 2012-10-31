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

#include "graph_filtering.hh"
#include "graph.hh"
#include "graph_properties.hh"

#include <boost/graph/make_maximal_planar.hpp>
#include <boost/graph/make_biconnected_planar.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

template <class EdgeMap>
struct mark_planar_edge
{
    mark_planar_edge(EdgeMap map, bool force): _map(map), _force(force) {}
    EdgeMap _map;
    bool _force;

    template <typename Graph, typename Vertex>
    void visit_vertex_pair(Vertex u, Vertex v, Graph& g)
    {
        if (!is_adjacent(u, v, g))
            add_edge(u, v, g);
    }

    template <typename Graph, typename Vertex, class EdgePredicate, class VertexPredicate>
    void visit_vertex_pair(Vertex u, Vertex v, UndirectedAdaptor<filtered_graph<Graph,
                                                                                EdgePredicate,
                                                                                VertexPredicate> >& g)
    {
        if (_force && !is_adjacent(u, v, g))
        {
            add_edge(u, v, g);
            return;
        }

        typedef typename graph_traits<Graph>::edge_descriptor edge_t;
        std::pair<edge_t, bool> e = edge(u, v, UndirectedAdaptor<Graph>(g.OriginalGraph().m_g));
        if (e.second)
            _map[e.first] = true;
    }

};

struct do_maximal_planar
{
    template <class Graph, class VertexIndex, class EdgeIndex, class EdgeMap>
    void operator()(Graph& g, VertexIndex vertex_index, EdgeIndex edge_index,
                    EdgeMap emap, bool augment) const
    {

        unchecked_vector_property_map
            <vector<typename graph_traits<Graph>::edge_descriptor>, VertexIndex>
            embedding(vertex_index, num_vertices(g));
        bool is_planar = boyer_myrvold_planarity_test
            (boyer_myrvold_params::graph = g,
             boyer_myrvold_params::edge_index_map = edge_index,
             boyer_myrvold_params::embedding = embedding);

        if (!is_planar)
            throw GraphException("Graph is not planar!");

        mark_planar_edge<EdgeMap> vis(emap, augment);
        make_biconnected_planar(g, embedding, edge_index, vis);
        boyer_myrvold_planarity_test
            (boyer_myrvold_params::graph = g,
             boyer_myrvold_params::edge_index_map = edge_index,
             boyer_myrvold_params::embedding = embedding);
        make_maximal_planar(g, embedding, vertex_index, edge_index, vis);
    }

};


void maximal_planar(GraphInterface& gi, boost::any edge_map, bool augment)
{
    if (augment)
    {
        run_action<graph_tool::detail::never_directed, mpl::true_>()
            (gi, bind<void>(do_maximal_planar(), _1, gi.GetVertexIndex(),
                            gi.GetEdgeIndex(), false, true))();
    }
    else
    {
        run_action<graph_tool::detail::never_directed, mpl::true_>()
            (gi, bind<void>(do_maximal_planar(), _1, gi.GetVertexIndex(),
                            gi.GetEdgeIndex(), _2, false),
             edge_scalar_properties())
            (edge_map);
    }
}
