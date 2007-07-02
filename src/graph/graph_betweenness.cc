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
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/mpl/for_each.hpp>

#include "graph.hh"
#include "histogram.hh"
#include "graph_filtering.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"
#include "graph_adaptor.hh"

#include <boost/graph/betweenness_centrality.hpp>


using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

template <class Graph, class VertexBetweenness, class EdgeBetweenness>
void normalize_betweenness(const Graph& g, VertexBetweenness vertex_betweenness, EdgeBetweenness edge_betweenness)
{
    size_t n = HardNumVertices()(g);
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
    template <class Graph, class VertexIndexMap, class EdgeBetweenness, class VertexBetweenness>
    void operator()(Graph& g, VertexIndexMap index_map, EdgeBetweenness edge_betweenness, VertexBetweenness vertex_betweenness) const
    {        
        brandes_betweenness_centrality(g, vertex_index_map(index_map).edge_centrality_map(edge_betweenness).centrality_map(vertex_betweenness));
        normalize_betweenness(g, vertex_betweenness, edge_betweenness);
    }
};

template <class Graph, class VertexIndexMap,  class EdgeIndexMap, class EdgeBetweenness, class VertexBetweenness>
struct get_weighted_betweenness
{
    get_weighted_betweenness(Graph& g, VertexIndexMap vertex_index, EdgeIndexMap edge_index, string& weight, dynamic_properties& dp,  
                             EdgeBetweenness edge_betweenness, VertexBetweenness vertex_betweenness, bool& found)
        :_g(g), _vertex_index(vertex_index), _edge_index(edge_index), _weight(weight), _dp(dp), _edge_betweenness(edge_betweenness), _vertex_betweenness(vertex_betweenness), _found(found) {}

    template <class WeightType>
    void operator()(WeightType)
    {
        try
        {
            typedef vector_property_map<WeightType, EdgeIndexMap> weight_map_t;
            weight_map_t weight_map(_edge_index);
            weight_map = get_static_property_map<weight_map_t>(find_property_map(_dp, _weight, typeid(graph_traits<typename GraphInterface::multigraph_t>::edge_descriptor)));
            brandes_betweenness_centrality(_g, vertex_index_map(_vertex_index).weight_map(weight_map).edge_centrality_map(_edge_betweenness).centrality_map(_vertex_betweenness));
            normalize_betweenness(_g, _vertex_betweenness, _edge_betweenness);
            _found = true;
        }
        catch (property_not_found) {}
        catch (bad_cast) {}        
    }
    
    Graph& _g;
    VertexIndexMap _vertex_index;
    EdgeIndexMap _edge_index;
    string& _weight;
    dynamic_properties& _dp;
    EdgeBetweenness _edge_betweenness;
    VertexBetweenness _vertex_betweenness;
    bool& _found;
};

struct choose_weight_map
{
    template <class Graph, class VertexIndexMap,  class EdgeIndexMap, class EdgeBetweenness, class VertexBetweenness>
    void operator()(Graph& g, VertexIndexMap vertex_index, EdgeIndexMap edge_index, string& weight, dynamic_properties& dp,  
                      EdgeBetweenness edge_betweenness, VertexBetweenness vertex_betweenness, bool& found) const
    {
        mpl::for_each<scalar_types>(get_weighted_betweenness<Graph,VertexIndexMap,EdgeIndexMap,EdgeBetweenness,VertexBetweenness>
                                    (g, vertex_index, edge_index, weight, dp, edge_betweenness, vertex_betweenness, found));
    }
};

void GraphInterface::GetBetweenness(string weight, string edge_betweenness, string vertex_betweenness)
{
    typedef vector_property_map<double, vertex_index_map_t> vertex_betweenness_map_t;
    vertex_betweenness_map_t vertex_betweenness_map(_vertex_index);

    typedef vector_property_map<double, edge_index_map_t> edge_betweenness_map_t;
    edge_betweenness_map_t edge_betweenness_map(_edge_index);

    if (vertex_betweenness != "" && edge_betweenness != "")
    {
        if (weight != "")
        {
            bool found = false;
            check_filter(*this, bind<void>(choose_weight_map(), _1, var(_vertex_index), var(_edge_index), var(weight), var(_properties), var(edge_betweenness_map), var(vertex_betweenness_map), var(found)), reverse_check(), directed_check());
            if (!found)
                throw GraphException("error getting scalar property: " + weight);
        }
        else
        {
            check_filter(*this, bind<void>(get_betweenness(), _1, var(_vertex_index),  var(edge_betweenness_map), var(vertex_betweenness_map)), reverse_check(), directed_check());
        }        
    }
    else if (vertex_betweenness != "")
    {
        if (weight != "")
        {
            bool found = false;
            check_filter(*this, bind<void>(choose_weight_map(), _1, var(_vertex_index), var(_edge_index), var(weight), var(_properties), var(edge_betweenness_map), dummy_property_map(), var(found)), reverse_check(), directed_check());
            if (!found)
                throw GraphException("error getting scalar property: " + weight);
        }
        else
        {
            check_filter(*this, bind<void>(get_betweenness(), _1, var(_vertex_index), var(edge_betweenness_map), dummy_property_map()), reverse_check(), directed_check());
        }
    }
    else
    {
        if (weight != "")
        {
            bool found = false;
            check_filter(*this, bind<void>(choose_weight_map(), _1, var(_vertex_index), var(_edge_index), var(weight), var(_properties), dummy_property_map(), var(vertex_betweenness_map), var(found)), reverse_check(), directed_check());
            if (!found)
                throw GraphException("error getting scalar property: " + weight);
        }
        else
        {
            check_filter(*this, bind<void>(get_betweenness(), _1, var(_vertex_index), dummy_property_map(), var(vertex_betweenness_map)), reverse_check(), directed_check());
        }        

    }

    if (vertex_betweenness != "")
    {
        try
        {
            find_property_map(_properties, vertex_betweenness, typeid(graph_traits<multigraph_t>::vertex_descriptor));
            RemoveVertexProperty(vertex_betweenness);
        }
        catch (property_not_found) {}

        _properties.property(vertex_betweenness, vertex_betweenness_map);
    }

    if (edge_betweenness != "")
    {
        try
        {
            find_property_map(_properties, edge_betweenness, typeid(graph_traits<multigraph_t>::edge_descriptor));
            RemoveEdgeProperty(edge_betweenness);
        }
        catch (property_not_found) {}

        _properties.property(edge_betweenness, edge_betweenness_map);
    }

}

struct get_central_point_dominance
{    
    template <class Graph, class VertexBetweenness>
    void operator()(Graph& g, VertexBetweenness vertex_betweenness, double& c) const
    {        
        c = central_point_dominance(g, vertex_betweenness);
    }
};


double GraphInterface::GetCentralPointDominance(string vertex_betweenness)
{
    try
    {
        typedef DynamicPropertyMapWrap<double, graph_traits<multigraph_t>::vertex_descriptor> betweenness_map_t;
        betweenness_map_t betweenness(find_property_map(_properties, vertex_betweenness, typeid(graph_traits<multigraph_t>::vertex_descriptor)));
        double c = 0.0;
        
        bool reversed = this->GetReversed();
        bool directed = this->GetDirected();
        check_filter(*this, bind<void>(get_central_point_dominance(), _1, var(betweenness), var(c)), never_reversed(), always_directed());
        this->SetReversed(reversed);
        this->SetDirected(directed);
        return c;
    }
    catch (property_not_found) 
    {
        throw GraphException("vertex property " + vertex_betweenness + " not found");
    }
}
