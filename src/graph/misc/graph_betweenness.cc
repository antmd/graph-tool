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

#include <boost/lambda/bind.hpp>
#include <boost/graph/betweenness_centrality.hpp>

#include "graph.hh"
#include "graph_selectors.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

template <class Graph, class EdgeBetweenness, class VertexBetweenness>
void normalize_betweenness(const Graph& g,
                           EdgeBetweenness edge_betweenness,
                           VertexBetweenness vertex_betweenness)
{
    size_t n = HardNumVertices()(&g);
    double factor = 2.0/(n*n - 3*n + 2);

    typename graph_traits<Graph>::vertex_iterator v, v_end;
    for (tie(v, v_end) = vertices(g); v != v_end; ++v)
    {
        put(vertex_betweenness, *v, factor * get(vertex_betweenness, *v));
    }

    typename graph_traits<Graph>::edge_iterator e, e_end;
    for (tie(e, e_end) = edges(g); e != e_end; ++e)
    {
        put(edge_betweenness, *e, factor * get(edge_betweenness, *e));
    }
}

struct get_betweenness
{
    template <class Graph, class VertexIndexMap, class EdgeBetweenness,
              class VertexBetweenness>
    void operator()(const Graph* gp, VertexIndexMap index_map,
                    EdgeBetweenness edge_betweenness,
                    VertexBetweenness vertex_betweenness) const
    {
        vector<vector<typename graph_traits<Graph>::edge_descriptor> >
            incoming_map(num_vertices(*gp));
        vector<size_t> distance_map(num_vertices(*gp));
        vector<typename property_traits<VertexBetweenness>::value_type>
            dependency_map(num_vertices(*gp));
        vector<size_t> path_count_map(num_vertices(*gp));

        brandes_betweenness_centrality
            (*gp, vertex_betweenness,edge_betweenness,
             make_iterator_property_map(incoming_map.begin(), index_map),
             make_iterator_property_map(distance_map.begin(), index_map),
             make_iterator_property_map(dependency_map.begin(), index_map),
             make_iterator_property_map(path_count_map.begin(), index_map),
             index_map);
        normalize_betweenness(*gp, edge_betweenness, vertex_betweenness);
    }
};

template <class VertexIndexMap>
struct get_weighted_betweenness
{
    get_weighted_betweenness(VertexIndexMap index): _vertex_index(index) {}

    template <class Graph, class EdgeBetweenness, class VertexBetweenness,
              class WeightMap>
    void operator()(const Graph* gp, EdgeBetweenness edge_betweenness,
                    VertexBetweenness vertex_betweenness,
                    WeightMap weight_map) const
    {
        vector<vector<typename graph_traits<Graph>::edge_descriptor> >
            incoming_map(num_vertices(*gp));
        vector<typename property_traits<WeightMap>::value_type>
            distance_map(num_vertices(*gp));
        vector<typename property_traits<VertexBetweenness>::value_type>
            dependency_map(num_vertices(*gp));
        vector<size_t> path_count_map(num_vertices(*gp));

        brandes_betweenness_centrality
            (*gp, vertex_betweenness,edge_betweenness,
             make_iterator_property_map(incoming_map.begin(), _vertex_index),
             make_iterator_property_map(distance_map.begin(), _vertex_index),
             make_iterator_property_map(dependency_map.begin(), _vertex_index),
             make_iterator_property_map(path_count_map.begin(), _vertex_index),
             _vertex_index, weight_map);
        normalize_betweenness(*gp, edge_betweenness, vertex_betweenness);
    }

    VertexIndexMap _vertex_index;
};

void GraphInterface::GetBetweenness(string weight, string edge_betweenness,
                                    string vertex_betweenness)
{
    boost::any edge_prop, vertex_prop;

    try
    {
        find_property_map(_properties, edge_betweenness, typeid(edge_t));
        edge_prop = prop(edge_betweenness, _edge_index, _properties);
        if (!belongs<edge_floating_properties>()(edge_prop))
            throw GraphException("edge property " + edge_betweenness +
                                 " is not of floating type");
    }
    catch (property_not_found)
    {
        typedef vector_property_map<double, edge_index_map_t>
            edge_betweenness_map_t;
        edge_betweenness_map_t edge_betweenness_map(_edge_index);
        edge_prop = edge_betweenness_map;
    }

    try
    {
        find_property_map(_properties, vertex_betweenness, typeid(vertex_t));
        vertex_prop = prop(vertex_betweenness, _vertex_index, _properties);
        if (!belongs<vertex_floating_properties>()(vertex_prop))
            throw GraphException("vertex property " + vertex_betweenness +
                                 " is not of floating type");
    }
    catch (property_not_found)
    {
        typedef vector_property_map<double, vertex_index_map_t>
            vertex_betweenness_map_t;
        vertex_betweenness_map_t vertex_betweenness_map(_vertex_index);
        vertex_prop = vertex_betweenness_map;
    }

    if (weight != "")
    {
        try
        {
            run_action<>()
                (*this,
                 get_weighted_betweenness<vertex_index_map_t>(_vertex_index),
                 edge_floating_properties(), vertex_floating_properties(),
                 edge_scalar_properties())(edge_prop, vertex_prop,
                                           prop(weight, _edge_index,
                                                _properties));
        }
        catch (property_not_found& e)
        {
            throw GraphException("edge property " + weight + " not found");
        }
    }
    else
    {
        run_action<>()(*this, bind<void>(get_betweenness(), _1,
                                         _vertex_index, _2, _3),
                       edge_floating_properties(), vertex_floating_properties())
            (edge_prop, vertex_prop);
    }
}

struct get_central_point_dominance
{
    template <class Graph, class VertexBetweenness>
    void operator()(Graph* g, VertexBetweenness vertex_betweenness, double& c)
        const
    {
        c = central_point_dominance(*g, vertex_betweenness);
    }
};

double GraphInterface::GetCentralPointDominance(string vertex_betweenness)
{
    try
    {
        double c = 0.0;

        bool directed = this->GetDirected();
        bool reversed = this->GetReversed();
        this->SetReversed(false);
        this->SetDirected(true);
        run_action<detail::never_reversed>()
            (*this, bind<void>(get_central_point_dominance(), _1, _2,
                               var(c)), vertex_scalar_properties())
            (prop(vertex_betweenness, _vertex_index, _properties));
        this->SetReversed(reversed);
        this->SetDirected(directed);
        return c;
    }
    catch (property_not_found)
    {
        throw GraphException("vertex property " + vertex_betweenness
                             + " not found");
    }
}
