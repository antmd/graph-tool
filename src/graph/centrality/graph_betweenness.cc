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

#include <boost/python.hpp>
#include <boost/graph/betweenness_centrality.hpp>

#include "graph.hh"
#include "graph_selectors.hh"
#include "graph_util.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

template <class Graph, class EdgeBetweenness, class VertexBetweenness>
void normalize_betweenness(const Graph& g,
                           EdgeBetweenness edge_betweenness,
                           VertexBetweenness vertex_betweenness,
                           size_t n)
{
    double vfactor = (n > 2) ? 1.0/((n-1)*(n-2)) : 1.0;
    double efactor = (n > 1) ? 1.0/(n*(n-1)) : 1.0;
    if (std::is_convertible<typename graph_traits<Graph>::directed_category,
                            undirected_tag>::value)
    {
        vfactor *= 2;
        efactor *= 2;
    }

    int i, N = num_vertices(g);
    #pragma omp parallel for default(shared) private(i)   \
        schedule(runtime) if (N > 100)
    for (i = 0; i < N; ++i)
    {
        typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
        if (v == graph_traits<Graph>::null_vertex())
            continue;
        put(vertex_betweenness, v, vfactor * get(vertex_betweenness, v));
    }

    typename graph_traits<Graph>::edge_iterator e, e_end;
    for (tie(e, e_end) = edges(g); e != e_end; ++e)
    {
        put(edge_betweenness, *e, efactor * get(edge_betweenness, *e));
    }
}

struct get_betweenness
{
    typedef void result_type;
    template <class Graph, class EdgeBetweenness, class VertexBetweenness>
    void operator()(Graph& g,
                    GraphInterface::vertex_index_map_t index_map,
                    EdgeBetweenness edge_betweenness,
                    VertexBetweenness vertex_betweenness,
                    bool normalize, size_t n) const
    {
        vector<vector<typename graph_traits<Graph>::edge_descriptor> >
            incoming_map(num_vertices(g));
        vector<size_t> distance_map(num_vertices(g));
        vector<typename property_traits<VertexBetweenness>::value_type>
            dependency_map(num_vertices(g));
        vector<size_t> path_count_map(num_vertices(g));
        brandes_betweenness_centrality
            (g, vertex_betweenness, edge_betweenness,
             make_iterator_property_map(incoming_map.begin(), index_map),
             make_iterator_property_map(distance_map.begin(), index_map),
             make_iterator_property_map(dependency_map.begin(), index_map),
             make_iterator_property_map(path_count_map.begin(), index_map),
             index_map);
        if (normalize)
            normalize_betweenness(g, edge_betweenness, vertex_betweenness, n);
    }
};

struct get_weighted_betweenness
{
    typedef void result_type;
    template <class Graph, class EdgeBetweenness, class VertexBetweenness,
              class VertexIndexMap>
        void operator()(Graph& g, VertexIndexMap vertex_index,
                        EdgeBetweenness edge_betweenness,
                        VertexBetweenness vertex_betweenness,
                        boost::any weight_map, bool normalize,
                        size_t n, size_t max_eindex) const
    {
        vector<vector<typename graph_traits<Graph>::edge_descriptor> >
            incoming_map(num_vertices(g));
        vector<typename property_traits<EdgeBetweenness>::value_type>
            distance_map(num_vertices(g));
        vector<typename property_traits<VertexBetweenness>::value_type>
            dependency_map(num_vertices(g));
        vector<size_t> path_count_map(num_vertices(g));

        typename EdgeBetweenness::checked_t weight =
            any_cast<typename EdgeBetweenness::checked_t>(weight_map);

        brandes_betweenness_centrality
            (g, vertex_betweenness, edge_betweenness,
             make_iterator_property_map(incoming_map.begin(), vertex_index),
             make_iterator_property_map(distance_map.begin(), vertex_index),
             make_iterator_property_map(dependency_map.begin(), vertex_index),
             make_iterator_property_map(path_count_map.begin(), vertex_index),
             vertex_index, weight.get_unchecked(max_eindex+1));
        if (normalize)
            normalize_betweenness(g, edge_betweenness, vertex_betweenness, n);
    }
};

void betweenness(GraphInterface& g, boost::any weight,
                 boost::any edge_betweenness,
                 boost::any vertex_betweenness,
                 bool normalize)
{
    if (!belongs<edge_floating_properties>()(edge_betweenness))
        throw ValueException("edge property must be of floating point value"
                             " type");

    if (!belongs<vertex_floating_properties>()(vertex_betweenness))
        throw ValueException("vertex property must be of floating point value"
                             " type");

    if (!weight.empty())
    {
        run_action<>()
            (g, std::bind<>(get_weighted_betweenness(),
                            std::placeholders::_1, g.GetVertexIndex(),
                            std::placeholders::_2,
                            std::placeholders::_3, weight, normalize,
                            g.GetNumberOfVertices(), g.GetMaxEdgeIndex()),
             edge_floating_properties(),
             vertex_floating_properties())
            (edge_betweenness, vertex_betweenness);
    }
    else
    {
        run_action<>()
            (g, std::bind<void>(get_betweenness(), std::placeholders::_1,
                                g.GetVertexIndex(), std::placeholders::_2,
                                std::placeholders::_3, normalize,
                                g.GetNumberOfVertices()),
             edge_floating_properties(),
             vertex_floating_properties())
            (edge_betweenness, vertex_betweenness);
    }
}

struct get_central_point_dominance
{
    template <class Graph, class VertexBetweenness>
    void operator()(Graph& g, VertexBetweenness vertex_betweenness, double& c)
        const
    {
        c = double(central_point_dominance(g, vertex_betweenness));
    }
};

double central_point(GraphInterface& g,
                     boost::any vertex_betweenness)
{
    double c = 0.0;
    run_action<graph_tool::detail::never_reversed>()
        (g, std::bind<>(get_central_point_dominance(), std::placeholders::_1,
                        std::placeholders::_2, std::ref(c)),
         vertex_scalar_properties()) (vertex_betweenness);
    return c;
}

void export_betweenness()
{
    using namespace boost::python;
    def("get_betweenness", &betweenness);
    def("get_central_point_dominance", &central_point);
}
