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

// based on code written by Alexandre Hannud Abdo <abdo@member.fsf.org>

#include <algorithm>
#include <tr1/unordered_set>
#include <boost/lambda/bind.hpp>
#include <boost/graph/breadth_first_search.hpp>

#include "graph.hh"
#include "histogram.hh"
#include "graph_filtering.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;


// filters out a single vertex
template <class Vertex>
struct single_vertex_filter 
{
    single_vertex_filter() {}
    single_vertex_filter(Vertex v):_v(v) {}

    bool operator()(Vertex v) const { return v != _v; }

    Vertex _v;
};

class bfs_stop_exception {};

// this will abort the BFS search when no longer useful

template <class TargetSet, class DistanceMap>
struct bfs_max_depth_watcher 
{
    typedef on_tree_edge event_filter;

    bfs_max_depth_watcher(TargetSet& targets, size_t max_depth, DistanceMap distance)
        : _targets(targets), _max_depth(max_depth), _distance(distance) {}
    
    template <class Graph>
    void operator()(typename graph_traits<Graph>::edge_descriptor e, const Graph& g) 
    {
        typename graph_traits<Graph>::vertex_descriptor v = target(e,g);
        if (get(_distance, v) > _max_depth)
            throw bfs_stop_exception();
        if (_targets.find(v) != _targets.end())
            _targets.erase(v);
        if (_targets.empty())
            throw bfs_stop_exception();
    }
    
    TargetSet& _targets;
    size_t _max_depth;
    DistanceMap _distance;
};


struct get_extended_clustering
{
    template <class Graph, class IndexMap, class ClusteringMap>
    void operator()(Graph& g, IndexMap vertex_index, vector<ClusteringMap>& cmaps) const
    {        

        int i, N = num_vertices(g);

        #pragma omp parallel for default(shared) private(i) schedule(dynamic) 
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
        
            // We must disconsider paths through the original vertex
            typedef single_vertex_filter<typename graph_traits<Graph>::vertex_descriptor> filter_t;
            typedef filtered_graph<Graph, keep_all, filter_t> fg_t;
            fg_t fg(g, keep_all(), filter_t(v));

            typedef DescriptorHash<IndexMap> hasher_t;
            typedef tr1::unordered_set<typename graph_traits<Graph>::vertex_descriptor,hasher_t> neighbour_set_t;
            neighbour_set_t neighbours(0, hasher_t(vertex_index));
            neighbour_set_t neighbours2(0, hasher_t(vertex_index));
            neighbour_set_t targets(0, hasher_t(vertex_index));
            
            // collect the targets
            typename graph_traits<Graph>::adjacency_iterator a, a_end;
            for(tie(a, a_end) = adjacent_vertices(v, g); a != a_end; ++a)
                if (*a != v) // no self-loops
                    targets.insert(*a);
            size_t k = targets.size();

            // And now we setup and start the BFS bonanza
            for(tie(a, a_end) = adjacent_vertices(v, g); a != a_end; ++a)
            {
                if (*a == v) // no self-loops
                    continue;
                if (neighbours.find(*a) != neighbours.end()) // avoid parallel edges
                    continue;
                neighbours.insert(*a);

                typedef tr1::unordered_map<typename graph_traits<Graph>::vertex_descriptor,size_t,DescriptorHash<IndexMap> > dmap_t;
                dmap_t dmap(0, DescriptorHash<IndexMap>(vertex_index));
                InitializedPropertyMap<dmap_t> distance_map(dmap, numeric_limits<size_t>::max());

                typedef tr1::unordered_map<typename graph_traits<Graph>::vertex_descriptor,default_color_type,DescriptorHash<IndexMap> > cmap_t;
                cmap_t cmap(0, DescriptorHash<IndexMap>(vertex_index));
                InitializedPropertyMap<cmap_t> color_map(cmap, color_traits<default_color_type>::white());
                
                try
                {
                    distance_map[*a] = 0;
                    neighbour_set_t specific_targets = targets;
                    specific_targets.erase(*a);
                    bfs_max_depth_watcher<neighbour_set_t,InitializedPropertyMap<dmap_t> > watcher(specific_targets, cmaps.size(), distance_map);
                    breadth_first_visit(fg, *a, visitor(make_bfs_visitor(make_pair(record_distances(distance_map, boost::on_tree_edge()),watcher))).
                                        color_map(color_map));
                }
                catch(bfs_stop_exception) {}

                neighbours2.clear();
                typename graph_traits<Graph>::adjacency_iterator a2;
                for(a2 = adjacent_vertices(v, g).first ; a2 != a_end ; ++a2) 
                {
                    if (*a2 == v || *a2 == *a) // no self-loops
                        continue;
                    if (neighbours2.find(*a2) != neighbours2.end()) // avoid parallel edges
                        continue;
                    neighbours2.insert(*a2);
                    if (distance_map[*a2] <= cmaps.size())
                        cmaps[distance_map[*a2]-1][v] += 1.0/(k*(k-1));
                }
            }
        }
    }
};


void GraphInterface::SetExtendedClusteringToProperty(string property_prefix, size_t max_depth)
{
    typedef vector_property_map<double, vertex_index_map_t> cmap_t;
    vector<cmap_t> cmaps(max_depth);
    for (size_t i = 0; i < cmaps.size(); ++i)
        cmaps[i] = cmap_t(_vertex_index);

    bool directed = _directed;
    _directed = false;
    check_filter(*this, bind<void>(get_extended_clustering(), _1, _vertex_index, var(cmaps)), reverse_check(), always_undirected()); 
    _directed = directed;

    for (size_t i = 0; i < cmaps.size(); ++i)
    {
        string name = property_prefix + lexical_cast<string>(i+1);
        try
        {
            find_property_map(_properties, name, typeid(graph_traits<multigraph_t>::vertex_descriptor));
            RemoveVertexProperty(name);
        }
        catch (property_not_found) {}

        _properties.property(name, cmaps[i]);
    }
}
