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

#include "graph_filtering.hh"
#include "graph.hh"
#include "graph_properties.hh"

#include <tr1/unordered_set>
#include <iostream>
#include <iomanip>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/device/file.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

// retrieves the line graph

struct get_line_graph
{
    template <class Graph, class VertexIndex, class LineGraph,
              class EdgeIndexMap, class LGVertexIndex>
    void operator()(const Graph& g, VertexIndex vertex_index,
                    LineGraph& line_graph, EdgeIndexMap edge_index,
                    LGVertexIndex vmap) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        typedef boost::property<edge_index_t, size_t> EdgeProperty;

        typedef typename property_map<LineGraph, vertex_index_t>::type
            line_vertex_index_map_t;
        line_vertex_index_map_t line_vertex_index(get(boost::vertex_index,
                                                      line_graph));

        typedef typename graph_traits<LineGraph>::vertex_descriptor lg_vertex_t;
        typedef HashedDescriptorMap<EdgeIndexMap,lg_vertex_t>
            edge_to_vertex_map_t;
        edge_to_vertex_map_t edge_to_vertex_map(edge_index);

        typename LGVertexIndex::checked_t vertex_map = vmap.get_checked();

        typename graph_traits<Graph>::edge_iterator e, e_end;
        for (tie(e, e_end) = edges(g); e != e_end; ++e)
        {
            typename graph_traits<LineGraph>::vertex_descriptor v;
            v = add_vertex(line_graph);
            edge_to_vertex_map[*e] = v;
            vertex_map[v] = edge_index[*e];
        }

        typedef typename property_map<LineGraph,edge_index_t>::type
            line_edge_index_map_t;
        line_edge_index_map_t line_edge_index(get(edge_index_t(), line_graph));

        size_t e_index = 0;
        if (boost::is_directed(g))
        {
            typename graph_traits<Graph>::vertex_iterator v, v_end;
            for (tie(v, v_end) = vertices(g); v != v_end; ++v)
            {
                typename graph_traits<Graph>::out_edge_iterator e1, e2, e1_end,
                    e2_end;;
                for (tie(e1, e1_end) = out_edges(*v, g); e1 != e1_end; ++e1)
                {
                    for (tie(e2, e2_end) = out_edges(target(*e1, g), g);
                         e2 != e2_end; ++e2)
                    {
                        typename graph_traits<LineGraph>::edge_descriptor
                            new_edge;
                        new_edge = add_edge(edge_to_vertex_map[*e1],
                                            edge_to_vertex_map[*e2],
                                            line_graph).first;
                        line_edge_index[new_edge] = e_index++;
                    }
                }
            }
        }
        else
        {
            typename graph_traits<Graph>::vertex_iterator v, v_end;
            for (tie(v, v_end) = vertices(g); v != v_end; ++v)
            {
                typename graph_traits<Graph>::out_edge_iterator e1, e2, e_end;
                for (tie(e1, e_end) = out_edges(*v, g); e1 != e_end; ++e1)
                    for (e2 = e1; e2 != e_end; ++e2)
                        if (*e1 != *e2)
                        {
                            typename graph_traits<LineGraph>::edge_descriptor
                                new_edge;
                            new_edge = add_edge(edge_to_vertex_map[*e1],
                                                edge_to_vertex_map[*e2],
                                                line_graph).first;
                            line_edge_index[new_edge] = e_index++;
                        }
            }
        }
    }
};

void line_graph(GraphInterface& gi, GraphInterface& lgi,
                boost::any edge_index)
{
    typedef property_map_types::apply<mpl::vector<int64_t>,
                                      GraphInterface::vertex_index_map_t,
                                       mpl::false_>::type
        vertex_properties;

     run_action<>()(gi, bind<void>(get_line_graph(), _1,
                                   gi.GetVertexIndex(),
                                   ref(lgi.GetGraph()), lgi.GetEdgeIndex(),
                                   _2),
                    vertex_properties())(edge_index);
     lgi.ReIndexEdges();
}
