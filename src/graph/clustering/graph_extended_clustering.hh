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

// based on code written by Alexandre Hannud Abdo <abdo@member.fsf.org>

#ifndef GRAPH_EXTENDED_CLUSTERING_HH
#define GRAPH_EXTENDED_CLUSTERING_HH

#include <unordered_set>

#include <boost/graph/breadth_first_search.hpp>

namespace graph_tool
{

using namespace std;
using namespace boost;

// graph filter to remove a single vertex

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

    bfs_max_depth_watcher(TargetSet& targets, size_t max_depth,
                          DistanceMap distance)
        : _targets(targets), _max_depth(max_depth), _distance(distance) {}

    template <class Graph>
    void operator()(typename graph_traits<Graph>::edge_descriptor e,
                    const Graph& g)
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

// abstract target collecting so algorithm works for bidirectional and
// undirected graphs

template<class Graph, class Vertex, class Targets, class DirectedCategory>
void collect_targets(Vertex v, Graph& g, Targets& t, DirectedCategory)
{
    typename graph_traits<Graph>::in_edge_iterator ei, ei_end;
    typename graph_traits<Graph>::vertex_descriptor u;
    for(tie(ei, ei_end) = in_edges(v, g); ei != ei_end; ++ei)
    {
        u = source(*ei, g);
        if (u == v) // no self-loops
            continue;
        if (t.find(u) != t.end()) // avoid parallel edges
            continue;
        t.insert(u);
    }
}

template<class Graph, class Vertex, class Targets>
void collect_targets(Vertex v, Graph& g, Targets& t, undirected_tag)
{
    typename graph_traits<Graph>::out_edge_iterator ei, ei_end;
    typename graph_traits<Graph>::vertex_descriptor u;
    for(tie(ei, ei_end) = out_edges(v, g); ei != ei_end; ++ei)
    {
        u = target(*ei, g);
        if (u == v) // no self-loops
            continue;
        if (t.find(u) != t.end()) // avoid parallel edges
            continue;
        t.insert(u);
    }
}

// get_extended_clustering

struct get_extended_clustering
{
    template <class Graph, class IndexMap, class ClusteringMap>
    void operator()(const Graph& g, IndexMap vertex_index,
                    vector<ClusteringMap> cmaps) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

        int i, N = num_vertices(g);

        #pragma omp parallel for default(shared) private(i) schedule(runtime) if (N > 100)
        for (i = 0; i < N; ++i)
        {
            vertex_t v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            // We must disconsider paths through the original vertex
            typedef single_vertex_filter<vertex_t> filter_t;
            typedef filtered_graph<Graph, keep_all, filter_t> fg_t;
            fg_t fg(g, keep_all(), filter_t(v));

            typedef DescriptorHash<IndexMap> hasher_t;
            typedef std::unordered_set<vertex_t,hasher_t> neighbour_set_t;
            neighbour_set_t neighbours(0, hasher_t(vertex_index));
            neighbour_set_t targets(0, hasher_t(vertex_index));
            typename neighbour_set_t::iterator ni, ti;

            // collect targets, neighbours and calculate normalization factor
            collect_targets(v, g, targets,
                            typename graph_traits<Graph>::directed_category());
            size_t k_in = targets.size(), k_out, k_inter=0, z;
            typename graph_traits<Graph>::adjacency_iterator a, a_end;
            for (tie(a, a_end) = adjacent_vertices(v, g); a != a_end; ++a)
            {
                if (*a == v) // no self-loops
                    continue;
                if (neighbours.find(*a) != neighbours.end()) // avoid parallel
                    continue;                                // edges

                neighbours.insert(*a);
                if (targets.find(*a) != targets.end())
                    ++k_inter;
            }
            k_out = neighbours.size();
            z = (k_in*k_out) - k_inter;

            // And now we setup and start the BFS bonanza
            for (ni = neighbours.begin(); ni != neighbours.end(); ++ni)
            {
                typedef std::unordered_map<vertex_t,size_t,
                                           DescriptorHash<IndexMap> > dmap_t;
                dmap_t dmap(0, DescriptorHash<IndexMap>(vertex_index));
                InitializedPropertyMap<dmap_t>
                    distance_map(dmap, numeric_limits<size_t>::max());

                typedef std::unordered_map<vertex_t,default_color_type,
                                           DescriptorHash<IndexMap> > cmap_t;
                cmap_t cmap(0, DescriptorHash<IndexMap>(vertex_index));
                InitializedPropertyMap<cmap_t>
                    color_map(cmap, color_traits<default_color_type>::white());

                try
                {
                    distance_map[*ni] = 0;
                    neighbour_set_t specific_targets = targets;
                    specific_targets.erase(*ni);
                    bfs_max_depth_watcher<neighbour_set_t,
                                          InitializedPropertyMap<dmap_t> >
                        watcher(specific_targets, cmaps.size(), distance_map);
                    breadth_first_visit(fg, *ni,
                                        visitor
                                        (make_bfs_visitor
                                         (make_pair(record_distances
                                                    (distance_map,
                                                     boost::on_tree_edge()),
                                                    watcher))).
                                        color_map(color_map));
                }
                catch (bfs_stop_exception) {}

                for (ti = targets.begin(); ti != targets.end(); ++ti)
                {
                    if (*ti == *ni) // no self-loops
                        continue;
                    if (distance_map[*ti] <= cmaps.size())
                        cmaps[distance_map[*ti]-1][v] += 1.0/z;
                }
            }
        }
    }
};

} // graph_tool namespace

#endif // GRAPH_EXTENDED_CLUSTERING_HH
