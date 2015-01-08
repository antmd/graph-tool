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
// you should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#ifndef GRAPH_CLUSTERING_HH
#define GRAPH_CLUSTERING_HH

#include "config.h"

#include <unordered_set>
#include <boost/mpl/if.hpp>

#ifdef HAVE_SPARSEHASH
#include SPARSEHASH_INCLUDE(dense_hash_set)
#endif

#ifndef __clang__
#include <ext/numeric>
using __gnu_cxx::power;
#else
template <class Value>
Value power(Value value, int n)
{
    return pow(value, n);
}
#endif

namespace graph_tool
{
using namespace boost;

#ifdef HAVE_SPARSEHASH
using google::dense_hash_set;
#else
using std::unordered_set;
#endif

// calculates the number of triangles to which v belongs
template <class Graph>
pair<int,int>
get_triangles(typename graph_traits<Graph>::vertex_descriptor v, const Graph &g)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

#ifdef HAVE_SPARSEHASH
    typedef dense_hash_set<vertex_t> set_t;
#else
    typedef unordered_set<vertex_t> set_t;
#endif

    set_t neighbour_set;

#ifdef HAVE_SPARSEHASH
     neighbour_set.set_empty_key(numeric_limits<vertex_t>::max());
     neighbour_set.resize(out_degree(v, g));
#endif

    size_t triangles = 0;

    typename graph_traits<Graph>::adjacency_iterator n, n_end;
    for (tie(n, n_end) = adjacent_vertices(v, g); n != n_end; ++n)
    {
        if (*n == v) // no self-loops
            continue;
        neighbour_set.insert(*n);
    }

    for (tie(n, n_end) = adjacent_vertices(v, g); n != n_end; ++n)
    {
        typename graph_traits<Graph>::adjacency_iterator n2, n2_end;
        for (tie(n2, n2_end) = adjacent_vertices(*n, g); n2 != n2_end; ++n2)
        {
            if (*n2 == *n) // no self-loops
                continue;
            if (neighbour_set.find(*n2) != neighbour_set.end())
                ++triangles;
        }
    }

    size_t k = out_degree(v, g);
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
            schedule(runtime) if (N > 100) reduction(+:triangles, n)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            temp = get_triangles(v, g);
            triangles += temp.first;
            n += temp.second;
        }
        c = double(triangles) / n;

        // "jackknife" variance
        c_err = 0.0;
        double cerr = 0.0;

        #pragma omp parallel for default(shared) private(i,temp) \
            schedule(runtime) if (N > 100) reduction(+:cerr)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            temp = get_triangles(v, g);
            double cl = double(triangles - temp.first) / (n - temp.second);

            cerr += power(c - cl, 2);
        }
        c_err = sqrt(cerr);
    }
};

// sets the local clustering coefficient to a property
struct set_clustering_to_property
{
    template <class Graph, class ClustMap>
    void operator()(const Graph& g, ClustMap clust_map) const
    {
        typedef typename property_traits<ClustMap>::value_type c_type;
        typename get_undirected_graph<Graph>::type ug(g);
        int i, N = num_vertices(g);

        #pragma omp parallel for default(shared) private(i) schedule(runtime) if (N > 100)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            pair<size_t,size_t> triangles = get_triangles(v,ug); // get from ug
            double clustering = (triangles.second > 0) ?
                double(triangles.first)/triangles.second :
                0.0;

            clust_map[v] = c_type(clustering);
        }
    }

    template <class Graph>
    struct get_undirected_graph
    {
        typedef typename mpl::if_
           <std::is_convertible<typename graph_traits<Graph>::directed_category,
                                directed_tag>,
            const UndirectedAdaptor<Graph>,
            const Graph& >::type type;
    };
};

} //graph-tool namespace

#endif // GRAPH_CLUSTERING_HH
