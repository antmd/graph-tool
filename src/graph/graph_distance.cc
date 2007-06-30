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

#include <algorithm>
#include <tr1/unordered_set>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/lambda/bind.hpp>
#include <iostream>

#include "graph.hh"
#include "histogram.hh"
#include "graph_filtering.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

//==============================================================================
// GetDistanceHistogram()
// retrieves the vertex-vertex distance histogram
//==============================================================================
struct no_weightS {};

struct get_distance_histogram
{

    template <class Graph, class IndexMap, class WeightMap, class Hist>
    void operator()(const Graph &g, IndexMap index_map, WeightMap weights, Hist& hist) const
    {        
        // select get_vertex_dists based on the existence of weights
        typedef typename mpl::if_<is_same<WeightMap, no_weightS>,
                                       get_dists_bfs,
                                  get_dists_djk>::type get_vertex_dists_t;

        get_vertex_dists_t get_vertex_dists;
        int i, N = num_vertices(g);

        #pragma omp parallel for default(shared) private(i) 
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            typedef tr1::unordered_map<typename graph_traits<Graph>::vertex_descriptor,double,DescriptorHash<IndexMap> > dmap_t;
            dmap_t dmap(0, DescriptorHash<IndexMap>(index_map));
            InitializedPropertyMap<dmap_t> dist_map(dmap, numeric_limits<double>::max());

            dist_map[v] = 0.0;
            get_vertex_dists(g, v, index_map, dist_map, weights);

            typename graph_traits<Graph>::vertex_iterator v2, v_end;
            for (tie(v2, v_end) = vertices(g); v2 != v_end; ++v2)
                if (*v2 != v && dist_map[*v2] != numeric_limits<double>::max())
                {
                    double dist = dist_map[*v2];
                    #pragma omp atomic
                    hist[dist]++;
                }
        }        
    }

    // weighted version. Use dijkstra_shortest_paths()
    struct get_dists_djk
    {
        template <class Graph, class Vertex, class IndexMap, class DistanceMap, class WeightMap>
        void operator()(const Graph& g, Vertex s, IndexMap index_map, DistanceMap dist_map, WeightMap weights) const
        {
            dijkstra_shortest_paths(g, s, vertex_index_map(index_map).weight_map(weights).distance_map(dist_map));  
        }
    };

    // unweighted version. Use BFS.
    struct get_dists_bfs
    {
        template <class Graph, class Vertex, class IndexMap, class DistanceMap>
        void operator()(const Graph& g, Vertex s, IndexMap index_map, DistanceMap dist_map, no_weightS) const
        {
            typedef tr1::unordered_map<typename graph_traits<Graph>::vertex_descriptor,default_color_type,DescriptorHash<IndexMap> > cmap_t;
            cmap_t cmap(0, DescriptorHash<IndexMap>(index_map));
            InitializedPropertyMap<cmap_t> color_map(cmap, color_traits<default_color_type>::white());
                        
            breadth_first_visit(g, s, visitor(make_bfs_visitor(record_distances(dist_map, on_tree_edge()))).color_map(color_map)); 
        }
    };

    
};

GraphInterface::hist_t GraphInterface::GetDistanceHistogram(string weight) const
{
    hist_t hist;

    if (weight == "")
    {
        check_filter(*this, bind<void>(get_distance_histogram(), _1, _vertex_index, no_weightS(), var(hist)),
                     reverse_check(), directed_check()); 
    }
    else
    {
        try 
        {
            dynamic_property_map& weight_prop = find_property_map(_properties, weight, typeid(graph_traits<multigraph_t>::edge_descriptor));
            try 
            {
                vector_property_map<double, edge_index_map_t> weight_map;
                weight_map = get_static_property_map<vector_property_map<double, edge_index_map_t> >(weight_prop);
                check_filter(*this, bind<void>(get_distance_histogram(), _1, _vertex_index, weight_map, var(hist)),
                             reverse_check(), directed_check()); 
            }
            catch (bad_cast)
            {
                DynamicPropertyMapWrap<double, graph_traits<multigraph_t>::edge_descriptor> weight_map(weight_prop);
                check_filter(*this, bind<void>(get_distance_histogram(), _1, _vertex_index, weight_map, var(hist)),
                             reverse_check(), directed_check()); 
            }
        }
        catch (property_not_found& e)
        {
            throw GraphException("error getting scalar property: " + string(e.what()));
        }
    }
    return hist;
}
