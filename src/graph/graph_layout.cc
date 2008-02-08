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

#include <boost/graph/gursoy_atun_layout.hpp>
#include <boost/graph/kamada_kawai_spring_layout.hpp>
#include <boost/graph/fruchterman_reingold.hpp>
#include <boost/graph/random_layout.hpp>
#include <boost/random.hpp>
#include <boost/lambda/lambda.hpp>

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

struct compute_gursoy
{
    template <class Graph, class PosMap, class WeightMap, class IndexMap>
    void operator()(const Graph* gp, size_t iter, size_t seed, PosMap pos,
                    WeightMap weight, string topology, IndexMap index_map) const
    {
        const Graph& g = *gp;
        mt19937 rng(static_cast<mt19937::result_type>(seed));
        size_t n = HardNumVertices()(&g);

        typename get_weight_map<WeightMap>::type weight_map(weight);

        if (iter == 0)
            iter = n;
        if (topology == "square")
        {
            square_topology<mt19937> top(rng, sqrt(n));
            gursoy_atun_layout(g, top,  pos, iter, sqrt(double(n)),
                               1.0, 0.8, 0.2, index_map, weight_map);
        }
        else if (topology == "circle")
        {
            circle_topology<mt19937> top(rng, n/2);
            gursoy_atun_layout(g, top,  pos, iter, sqrt(double(n)),
                               1.0, 0.8, 0.2, index_map, weight_map);
        }
        else if (topology == "heart")
        {
            heart_topology<mt19937> top(rng);
            vector_property_map<heart_topology<mt19937>::point_type,
                                IndexMap> pos_map (index_map);
            gursoy_atun_layout(g, top, pos_map, iter, sqrt(double(n)),
                              1.0, 0.8, 0.2, index_map, weight_map);
            typename graph_traits<Graph>::vertex_iterator v, v_end;
            for (tie(v, v_end) = vertices(g); v != v_end; ++v)
            {
                pos[*v][0] = pos_map[*v][0];
                pos[*v][1] = pos_map[*v][1];
            }
        }
        else
        {
            throw GraphException("invalid topology for graph layout: " +
                                 topology);
        }
    }

    template <class PropertyMap>
    struct get_weight_map
    {
        typedef typename property_traits<PropertyMap>::value_type
            value_type;
        typedef typename mpl::if_<
            typename mpl::or_<is_floating_point<value_type>,
                              is_same<PropertyMap,dummy_property_map> >::type,
            PropertyMap,
            ConvertedPropertyMap<PropertyMap,double> >::type type;
    };
};

struct copy_points
{
    template <class Graph, class VectorProperty, class PointsMap>
    void operator()(const Graph* g, VectorProperty vec_prop, PointsMap pos_map)
        const
    {
        typename graph_traits<Graph>::vertex_iterator v, v_end;
        for (tie(v, v_end) = vertices(*g); v != v_end; ++v)
        {
            vec_prop[*v].resize(2);
            vec_prop[*v][0] = pos_map[*v][0];
            vec_prop[*v][1] = pos_map[*v][1];
        }
    }
};

struct copy_weights
{
    template <class Graph, class WeightSrcMap, class WeightTgtMap>
    void operator()(const Graph* g, WeightSrcMap src_map, WeightTgtMap tgt_map)
        const
    {
        do_copy(g, src_map, tgt_map, is_same<WeightSrcMap, WeightTgtMap>());
    }

    template <class Graph, class WeightSrcMap, class WeightTgtMap>
    void do_copy(const Graph* g, WeightSrcMap src_map, WeightTgtMap tgt_map,
                 true_type) const
    {
        tgt_map = src_map;
    }

    template <class Graph, class WeightSrcMap, class WeightTgtMap>
    void do_copy(const Graph* g, WeightSrcMap src_map, WeightTgtMap tgt_map,
                 false_type) const
    {
        typename graph_traits<Graph>::edge_iterator e, e_end;
        for (tie(e, e_end) = edges(*g); e != e_end; ++e)
        {
            tgt_map[*e] = src_map[*e];
        }
    }
};


struct compute_spring_block
{
    template <class Graph, class PosMap, class IndexMap, class WeightMap>
    void operator()(const Graph* gp, string type, size_t iter, size_t seed,
                    PosMap pos, bool progressive, WeightMap weight,
                    IndexMap index_map)
        const
    {
        const Graph& g = *gp;
        mt19937 rng(static_cast<mt19937::result_type>(seed));
        size_t n = HardNumVertices()(&g);

        if (iter == 0)
            iter = 100;
        double radius = n;
        if (!progressive)
            random_graph_layout(g, pos, -radius/2, radius/2,
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


struct copy_positions
{
    template <class Graph, class VectorProperty, class PosMap>
    void operator()(const Graph* g, VectorProperty vec_prop, PosMap pos_map)
        const
    {
        typename graph_traits<Graph>::vertex_iterator v, v_end;
        for (tie(v, v_end) = vertices(*g); v != v_end; ++v)
        {
            vec_prop[*v].resize(2);
            vec_prop[*v][0] = pos_map[*v].x;
            vec_prop[*v][1] = pos_map[*v].y;
        }
    }
};

void GraphInterface::ComputeGraphLayoutGursoy(string pos_prop, string weight,
                                              string topology, size_t iter,
                                              size_t seed)
{
    typedef vector_property_map<double, edge_index_map_t> temp_weight_t;
    boost::any pos_map, weight_map;
    try
    {
        pos_map = prop(pos_prop, _vertex_index, _properties);
        if (!belongs<vertex_floating_vector_properties>()(pos_map))
            throw GraphException("property " + pos_prop +
                                 "is not of floating point vector type");
        if (weight == "")
        {
            weight_map = dummy_property_map();
        }
        else
        {
            temp_weight_t temp_weight(_edge_index);
            run_action<>()(*this, bind<void>(copy_weights(), _1, _2,
                                             temp_weight),
                           edge_scalar_properties())
                (prop(weight, _edge_index, _properties));
            weight_map = temp_weight;
        }
    }
    catch (property_not_found& e)
    {
        throw GraphException("error getting property: " + string(e.what()));
    }

    vector_property_map<convex_topology<2>::point_type, vertex_index_map_t>
        pos(_vertex_index);
    run_action<>()(*this, bind<void>(compute_gursoy(), _1, iter, seed, pos,
                                     _2, topology, _vertex_index),
                   mpl::vector<temp_weight_t,dummy_property_map>())
        (weight_map);
    run_action<>()(*this, bind<void>(copy_points(), _1, _2, pos),
                   vertex_floating_vector_properties())(pos_map);
}

struct pos_t
{
    double x;
    double y;
};

void GraphInterface::ComputeGraphLayoutSpringBlock(string pos_prop,
                                                   string weight, string type,
                                                   size_t iter,
                                                   bool progressive,
                                                   size_t seed)
{
    typedef vector_property_map<double, edge_index_map_t> temp_weight_t;
    boost::any pos_map, weight_map;
    try
    {
        pos_map = prop(pos_prop, _vertex_index, _properties);
        if (!belongs<vertex_floating_vector_properties>()(pos_map))
            throw GraphException("property " + pos_prop +
                                 "is not of floating point vector type");

        if (weight == "")
        {
            weight_map = ConstantPropertyMap<double,edge_t>(1.0);
        }
        else
        {
            temp_weight_t temp_weight(_edge_index);
            run_action<>()(*this, bind<void>(copy_weights(), _1, _2,
                                             temp_weight),
                           edge_scalar_properties())
                (prop(weight, _edge_index, _properties));
            weight_map = temp_weight;
        }
    }
    catch (property_not_found& e)
    {
        throw GraphException("error getting property: " + string(e.what()));
    }

    typedef mpl::vector<temp_weight_t, ConstantPropertyMap<double,edge_t> >
        weight_maps;

    vector_property_map<pos_t, vertex_index_map_t> pos(_vertex_index);
    run_action<>()(*this, bind<void>(compute_spring_block(), _1, type, iter,
                                     seed, pos, progressive, _2, _vertex_index),
                   weight_maps())
        (weight_map);
    run_action<>()(*this, bind<void>(copy_positions(), _1, _2, pos),
                   vertex_floating_vector_properties())(pos_map);
}
