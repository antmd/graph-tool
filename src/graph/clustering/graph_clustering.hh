// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007-2010 Tiago de Paula Peixoto <tiago@forked.de>
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

#ifndef GRAPH_CLUSTERING_HH
#define GRAPH_CLUSTERING_HH

#include <tr1/unordered_set>
#include <boost/mpl/if.hpp>

namespace graph_tool
{
using namespace boost;

// calculates the number of triangles to which v belongs
template <class Graph>
pair<int,int>
get_triangles(typename graph_traits<Graph>::vertex_descriptor v, const Graph &g)
{
    tr1::unordered_set<typename graph_traits<Graph>::vertex_descriptor>
        neighbour_set1, neighbour_set2, neighbour_set3;

    size_t triangles = 0, k = 0;

    typename graph_traits<Graph>::adjacency_iterator n1_begin, n1_end, n1;
    tie(n1_begin, n1_end) = adjacent_vertices(v, g);
    for (n1 = n1_begin; n1 != n1_end; ++n1)
    {
        if (*n1 == v) // no self-loops
            continue;
        if (neighbour_set1.find(*n1) != neighbour_set1.end())
            continue;
        else
            neighbour_set1.insert(*n1);

        typename graph_traits<Graph>::adjacency_iterator n2_begin, n2_end, n2;
        tie(n2_begin, n2_end) = adjacent_vertices(*n1, g);
        for (n2 = n2_begin; n2 != n2_end; ++n2)
        {
            if (*n2 == *n1) // no self-loops
                continue;
            if (neighbour_set2.find(*n2) != neighbour_set2.end())
                continue;
            else
                neighbour_set2.insert(*n2);

            typename graph_traits<Graph>::adjacency_iterator
                n3_begin, n3_end, n3;
            tie(n3_begin, n3_end) = adjacent_vertices(*n2, g);
            for (n3 = n3_begin; n3 != n3_end; ++n3)
            {
                if (*n3 == *n2) // no self-loops
                    continue;
                if (neighbour_set3.find(*n3) != neighbour_set3.end())
                    continue;
                else
                    neighbour_set3.insert(*n3);

                if (*n3 == v) //found a triangle
                    triangles++;
            }
            neighbour_set3.clear();
        }
        neighbour_set2.clear();
        k++;
    }
    neighbour_set1.clear();
    return make_pair(triangles/2,(k*(k-1))/2);
}


// retrieves the global clustering coefficient
struct get_global_clustering
{
    template <class Graph>
    void operator()(const Graph& g, double& c, double& c_err) const
    {
        size_t triangles = 0, n = 0;
        pair<size_t, size_t> temp;

        int i, N = num_vertices(g);

        #pragma omp parallel for default(shared) private(i,temp) \
            schedule(dynamic) reduction(+:triangles)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            temp = get_triangles(v, g);
            triangles += temp.first;
            n += temp.second;
        }
        c = double(triangles)/n;

        // "jackknife" variance
        c_err = 0.0;

        #pragma omp parallel for default(shared) private(i,temp) \
            schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            temp = get_triangles(v, g);
            double cl = double(triangles - temp.first)/(n - temp.second);

            #pragma omp atomic
            c_err += (c - cl)*(c - cl);
        }
        c_err = sqrt(c_err);
    }
};

// sets the local clustering coefficient to a property
struct set_clustering_to_property
{
    template <class Graph, class ClustMap>
    void operator()(const Graph& g, ClustMap clust_map) const
    {
        typename get_undirected_graph<Graph>::type ug(g);
        int i, N = num_vertices(g);

        #pragma omp parallel for default(shared) private(i) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            pair<size_t,size_t> triangles = get_triangles(v,ug); // get from ug
            double clustering = (triangles.second > 0) ?
                double(triangles.first)/triangles.second :
                0.0;

            #pragma omp critical
            {
                clust_map[v] = clustering;
            }
        }
    }

    template <class Graph>
    struct get_undirected_graph
    {
        typedef typename mpl::if_
            < is_convertible<typename graph_traits<Graph>::directed_category,
                             directed_tag>,
              const UndirectedAdaptor<Graph>,
              const Graph& >::type type;
    };
};

} //graph-tool namespace

#endif // GRAPH_CLUSTERING_HH
