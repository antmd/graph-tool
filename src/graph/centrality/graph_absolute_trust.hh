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

#ifndef GRAPH_ABSOLUTE_TRUST_HH
#define GRAPH_ABSOLUTE_TRUST_HH

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_util.hh"

#include <tr1/unordered_set>
#include <tr1/tuple>
#include <algorithm>

#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <iostream>

namespace graph_tool
{
using namespace std;
using namespace boost;


struct stop_search {}; // exception to be thrown to stop search

template <class SourceMap, class WeightMap>
class source_counter:
    public boost::dijkstra_visitor<null_visitor>
{
public:
    source_counter(SourceMap source_map, WeightMap weight_map, size_t n_sources)
        : _source_map(source_map), _weight_map(weight_map),
          _n_sources(n_sources) {}

    template <class Graph>
    void examine_vertex(typename graph_traits<Graph>::vertex_descriptor u,
                        Graph& g)
    {
        // stop if all sources are found
        if (_source_map[u])
        {
            _n_sources--;
            if (_n_sources == 0)
                throw stop_search();
        }
    }

private:
    SourceMap _source_map;
    WeightMap _weight_map;
    size_t _n_sources;
};

struct dist_compare
{
    template <class Type1, class Type2>
    bool operator()(const Type1& d1, const Type2& d2) const
    {
        return d1 > d2; // we want trust paths with _maximum_ "distance"
    }
};

struct dist_combine
{
    template <class DistType, class WeightType>
    DistType operator()(const DistType& d, const WeightType& w) const
    {
        return d*w;
    }
};

// predicate to filter a single vertex from the graph
struct filter_vertex_pred
{
    filter_vertex_pred() {}
    filter_vertex_pred(size_t v): _v(v) {}
    template <class Vertex>
    bool operator()(const Vertex& v) const { return v != _v; }
    size_t _v;
};

struct get_absolute_trust
{
    template <class Graph, class VertexIndex, class EdgeIndex, class TrustMap,
              class InferredTrustMap>
    void operator()(Graph& g, VertexIndex vertex_index, EdgeIndex edge_index,
                    size_t max_edge_index, int64_t source, int64_t target,
                    TrustMap c, InferredTrustMap t) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;
        typedef typename property_traits<TrustMap>::value_type c_type;
        typedef typename property_traits<InferredTrustMap>::value_type t_type;

        vertex_t v = vertex(source, g);

        int i, N = (target == -1) ? num_vertices(g) : target + 1;
        #pragma omp parallel for default(shared) private(i) schedule(dynamic)
        for (i = (target == -1) ? 0 : target; i < N; ++i)
        {

            vertex_t w = vertex(i, g);

            // mark the sources
            typedef unchecked_vector_property_map<uint8_t, VertexIndex>
                source_map_t;
            source_map_t source_map(vertex_index, num_vertices(g));

            typename in_edge_iteratorS<Graph>::type e, e_end;
            for (tie(e, e_end) = in_edge_iteratorS<Graph>::get_edges(w, g);
                 e != e_end; ++e)
                source_map[boost::source(*e,g)] = true;

            // filter vertex w out of the graph
            filtered_graph<Graph, boost::keep_all, filter_vertex_pred>
                fg(g, boost::keep_all(), filter_vertex_pred(w));

            // distance map (target weight map)
            typedef unchecked_vector_property_map<t_type, VertexIndex>
                dist_map_t;
            dist_map_t dist_map(vertex_index, num_vertices(g));

            // compute the targets weights
            try
            {
                size_t k = in_degreeS()(w,g);
                source_counter<source_map_t,dist_map_t>
                    visitor(source_map, dist_map, k);
                dijkstra_shortest_paths(fg, v, weight_map(c).
                                        vertex_index_map(vertex_index).
                                        distance_map(dist_map).
                                        distance_compare(dist_compare()).
                                        distance_combine(dist_combine()).
                                        distance_inf(t_type(0)).
                                        distance_zero(t_type(1)).
                                        visitor(visitor));
            }
            catch (const stop_search&) {}

            // compute the target's trust
            t_type sum_w = 0, avg = 0;
            for (tie(e, e_end) = in_edge_iteratorS<Graph>::get_edges(w, g);
                 e != e_end; ++e)
            {
                t_type weight = dist_map[boost::source(*e,g)];
                sum_w += weight;
                avg += c[*e]*weight*weight;
            }
            if (sum_w > 0)
                t[w] = avg/sum_w;
        }
    }
};

}

#endif
