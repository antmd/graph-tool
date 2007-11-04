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


#include <algorithm>
#include <tr1/unordered_set>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/graph/gursoy_atun_layout.hpp>
#include <boost/graph/kamada_kawai_spring_layout.hpp>
#include <boost/graph/fruchterman_reingold.hpp>
#include <boost/graph/random_layout.hpp>
#include <boost/random.hpp>

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_properties.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;


struct compute_gursoy
{
    template <class Graph, class PosMap, class WeightMap, class IndexMap>
    void operator()(Graph &g, size_t iter, size_t seed, PosMap pos,
                    WeightMap weight, string topology, IndexMap index_map) const
    {
        mt19937 rng(static_cast<mt19937::result_type>(seed));
        size_t n = HardNumVertices()(g);

        vector_property_map<square_topology<mt19937>::point_type, IndexMap>
            position_map(index_map);
        if (iter == 0)
            iter = n;
        if (topology == "square")
        {
            square_topology<mt19937> top(rng, sqrt(n));
            gursoy_atun_layout (g, top,  position_map, iter, sqrt(double(n)),
                                1.0, 0.8, 0.2, index_map, weight);
        }
        else if (topology == "circle")
        {
            circle_topology<mt19937> top(rng, n/2);
            gursoy_atun_layout (g, top,  position_map, iter, sqrt(double(n)),
                                1.0, 0.8, 0.2, index_map, weight);
        }
        else if (topology == "heart")
        {
            heart_topology<mt19937> top(rng);
            vector_property_map<heart_topology<mt19937>::point_type, IndexMap>
                pos_map(index_map);
            gursoy_atun_layout (g, top, pos_map, iter, sqrt(double(n)), 1.0,
                                0.8, 0.2, index_map, weight);
            int i, N = num_vertices(g);
            #pragma omp parallel for default(shared) private(i) \
                schedule(dynamic)
            for (i = 0; i < N; ++i)
            {
                typename graph_traits<Graph>::vertex_descriptor v =
                    vertex(i, g);
                if (v == graph_traits<Graph>::null_vertex())
                    continue;
                position_map[v][0] = pos_map[v][0];
                position_map[v][1] = pos_map[v][1];
            }
        }
        else
        {
            throw GraphException("invalid topology for graph layout: " +
                                 topology);
        }

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            pos[v].x = position_map[v][0];
            pos[v].y = position_map[v][1];
        }
    }
};


void GraphInterface::ComputeGraphLayoutGursoy(string prop, string weight,
                                              string topology, size_t iter,
                                              size_t seed)
{
    // vertex postion map
    typedef vector_property_map<pos_t, vertex_index_map_t> pos_map_t;
    pos_map_t pos_map(num_vertices(_mg), _vertex_index);

    if (weight == "")
    {
        check_filter(*this,bind<void>(compute_gursoy(),_1,iter,seed, pos_map,
                                      dummy_property_map(), topology,
                                      _vertex_index),
                     reverse_check(), directed_check());
    }
    else
    {
        try
        {
            dynamic_property_map& weight_prop =
                find_property_map(_properties, weight, typeid(edge_t));
            try
            {
                vector_property_map<double, edge_index_map_t> weight_map;
                weight_map =
                    get_static_property_map
                        <vector_property_map<double, edge_index_map_t> >
                            (weight_prop);
                check_filter(*this,bind<void>(compute_gursoy(),_1,iter,seed,
                                              pos_map, weight_map, topology,
                                              _vertex_index),
                             reverse_check(),directed_check());
            }
            catch (bad_cast)
            {
                DynamicPropertyMapWrap<double, edge_t> weight_map(weight_prop);
                check_filter(*this,bind<void>(compute_gursoy(),_1,iter,seed,
                                              pos_map, weight_map, topology,
                                              _vertex_index),reverse_check(),
                             directed_check());
            }
        }
        catch (property_not_found& e)
        {
            throw GraphException("error getting scalar property: " +
                                 string(e.what()));
        }
    }

    try
    {
        find_property_map(_properties, prop, typeid(vertex_t));
        RemoveVertexProperty(prop);
    }
    catch (property_not_found) {}
    _properties.property(prop, pos_map);
}

struct compute_spring_block
{
    template <class Graph, class PosMap, class IndexMap, class WeightMap>
    void operator()(Graph &g, string type, size_t iter, size_t seed, PosMap pos,
                    bool progressive, WeightMap weight, IndexMap index_map) const
    {
        mt19937 rng(static_cast<mt19937::result_type>(seed));
        size_t n = HardNumVertices()(g);

        if (iter == 0)
            iter = 100;
        double radius = n;
        if (!progressive)
            random_graph_layout(g, pos,
                                -radius/2, radius/2,
                                -radius/2, radius/2, rng);
        if (type == "fg-grid")
        {
            fruchterman_reingold_force_directed_layout
                (g, pos, radius, radius,
                 cooling(linear_cooling<double>(iter)).
                 vertex_index_map(index_map));
        }
        else if (type == "fg-all-pairs")
        {
            fruchterman_reingold_force_directed_layout
                (g, pos, radius, radius,
                 cooling(linear_cooling<double>(iter)).
                 vertex_index_map(index_map).force_pairs(all_force_pairs()));
        }
        else if (type == "kw")
        {
            typedef typename graph_traits<Graph>::directed_category
                directed_category;
            bool retval;
            retval = compute_kamada_kawai(g, iter, n, pos, weight,
                                          typename is_convertible
                                                      <directed_category,
                                                       undirected_tag>::type());
            if (!retval)
                throw GraphException("the Kamada-Kawai layout algorithm only "
                                     "works for connected graphs!");
        }
        else
        {
            throw GraphException("invalid type of spring-block graph layout: " +
                                 type);
        }
    }

    template <class Graph, class PosMap, class WeightMap>
    bool compute_kamada_kawai(Graph &g, size_t iter, size_t n, PosMap pos,
                              WeightMap weight, boost::true_type) const
    {
        return kamada_kawai_spring_layout(g, pos, weight, side_length(n));
    }

    template <class Graph, class PosMap, class WeightMap>
    bool compute_kamada_kawai(Graph &g,  size_t iter, size_t n, PosMap pos,
                              WeightMap weight, boost::false_type) const
    {
        UndirectedAdaptor<Graph> ug(g);
        return kamada_kawai_spring_layout(ug, pos, weight, side_length(n));
    }
};

void GraphInterface::ComputeGraphLayoutSpringBlock(string prop, string weight,
                                                   string type, size_t iter,
                                                   size_t seed)
{
    // vertex postion map
    typedef vector_property_map<pos_t, vertex_index_map_t> pos_map_t;
    pos_map_t pos_map(num_vertices(_mg), vertex_index);
    bool progressive = false;
    try
    {
        dynamic_property_map& pos = find_property_map(_properties, weight,
                                                      typeid(edge_t));
        pos_map =
            get_static_property_map
                <vector_property_map<pos_t,vertex_index_map_t> >(pos);
        progressive = true;
    }
    catch (bad_cast) {}
    catch (property_not_found) {}

    if (weight == "")
    {
        check_filter(*this,bind<void>(compute_spring_block(), _1, type, iter,
                                      seed, pos_map, progressive,
                                      ConstantPropertyMap<double,edge_t>(1.0),
                                      _vertex_index),
                     reverse_check(), directed_check());
    }
    else
    {
        try
        {
            dynamic_property_map& weight_prop =
                find_property_map(_properties, weight, typeid(edge_t));
            try
            {
                vector_property_map<double, edge_index_map_t> weight_map;
                weight_map =
                    get_static_property_map
                        <vector_property_map<double, edge_index_map_t> >
                            (weight_prop);
                check_filter(*this,bind<void>(compute_spring_block(), _1, type,
                                              iter, seed, pos_map, progressive,
                                              weight_map, _vertex_index),
                             reverse_check(), directed_check());
            }
            catch (bad_cast)
            {
                DynamicPropertyMapWrap<double,edge_t> weight_map(weight_prop);
                check_filter(*this,bind<void>(compute_spring_block(), _1, type,
                                              iter, seed, pos_map, progressive,
                                              weight_map, _vertex_index),
                             reverse_check(), directed_check());
            }
        }
        catch (property_not_found& e)
        {
            throw GraphException("error getting scalar property: " +
                                 string(e.what()));
        }
    }

    try
    {
        find_property_map(_properties, prop, typeid(vertex_t));
        RemoveVertexProperty(prop);
    }
    catch (property_not_found) {}
    _properties.property(prop, pos_map);
}
