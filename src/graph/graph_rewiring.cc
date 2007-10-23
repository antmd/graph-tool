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
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/random.hpp>
#include <boost/functional/hash.hpp>

#include "graph.hh"
#include "graph_filtering.hh"

#include "shared_map.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

typedef boost::mt19937 rng_t;

struct get_in_degree
{
    template <class Graph>
    size_t operator()(typename graph_traits<Graph>::vertex_descriptor v, const Graph& g)
    {
        return get_degree(v, g, typename is_convertible<typename graph_traits<Graph>::directed_category,directed_tag>::type());
    }

    template <class Graph>
    size_t get_degree(typename graph_traits<Graph>::vertex_descriptor v, const Graph& g, true_type)
    {
        return in_degree(v,g);
    }

    template <class Graph>
    size_t get_degree(typename graph_traits<Graph>::vertex_descriptor v, const Graph& g, false_type)
    {
        return out_degree(v,g);
    }
};

template <class Graph>
struct get_in_edges
{
    typedef typename mpl::if_<typename is_convertible<typename graph_traits<Graph>::directed_category,directed_tag>::type,
                              typename graph_traits<Graph>::in_edge_iterator,
                              typename graph_traits<Graph>::out_edge_iterator>::type iterator;

    pair<iterator,iterator> operator()(typename graph_traits<Graph>::vertex_descriptor v, const Graph& g)
    {
        return get_edges(v, g, typename is_convertible<typename graph_traits<Graph>::directed_category,directed_tag>::type());
    }

    pair<iterator,iterator> get_edges(typename graph_traits<Graph>::vertex_descriptor v, const Graph& g, true_type)
    {
        return in_edges(v,g);
    }

    pair<iterator,iterator> get_edges(typename graph_traits<Graph>::vertex_descriptor v, const Graph& g, false_type)
    {
        return out_edges(v,g);
    }
};

template <class Graph>
bool is_adjacent(typename graph_traits<Graph>::vertex_descriptor u, typename graph_traits<Graph>::vertex_descriptor v, Graph& g )
{
    typename graph_traits<Graph>::out_edge_iterator e, e_end;
    for (tie(e, e_end) = out_edges(u, g); e != e_end; ++e)
    {
        if (target(*e,g) == v)
            return true;
    }
    return false;
}

template <template <class Graph> class RewireStrategy>
struct graph_rewire
{
    template <class Graph, class EdgeIndexMap>
    void operator()(const Graph& g, EdgeIndexMap edge_index, size_t seed, bool self_loops, bool parallel_edges) const
    {
        rng_t rng(static_cast<rng_t::result_type>(seed));

        RewireStrategy<Graph> rewire(g, rng);
       
        if (!self_loops)
        {
            // check the existence of self-loops
            bool has_self_loops = false;
            int i, N = num_vertices(g);
            #pragma omp parallel for default(shared) private(i) schedule(dynamic)
            for (i = 0; i < N; ++i)
            {
                typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
                if (v == graph_traits<Graph>::null_vertex())
                    continue;
                if (is_adjacent(v, v, g))
                    self_loops = true;
            }
            if (has_self_loops)
                throw GraphException("Self-loop detected. Can't rewire graph without self-loops if it already contains self-loops!");
        }

        if (!parallel_edges)
        {
            // check the existence of parallel edges 
            bool has_parallel_edges = false;
            int i, N = num_vertices(g);
            #pragma omp parallel for default(shared) private(i) schedule(dynamic)
            for (i = 0; i < N; ++i)
            {
                typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
                if (v == graph_traits<Graph>::null_vertex())
                    continue;

                tr1::unordered_set<typename graph_traits<Graph>::vertex_descriptor> targets;
                typename graph_traits<Graph>::out_edge_iterator e, e_end;
                for (tie(e, e_end) = out_edges(v, g); e != e_end; ++e)
                {
                    if (targets.find(target(*e, g)) != targets.end())
                        has_parallel_edges = true;
                    else
                        targets.insert(target(*e, g));
                }
            }

            if (has_parallel_edges)
                throw GraphException("Parallel edge detected. Can't rewire graph without parallel edges if it already contains parallel edges!");
        }

        vector<typename graph_traits<Graph>::edge_descriptor> edges(num_edges(g));
        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            typename graph_traits<Graph>::out_edge_iterator e, e_end;
            for (tie(e, e_end) = out_edges(v, g); e != e_end; ++e)
                edges[edge_index[*e]] = *e;
        }

        vector<pair<size_t, size_t> > swap_list;
        swap_list.reserve(num_edges(g));

        // for each vertex, rewire the targets of its out edges
        N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            
            tr1::unordered_set<typename graph_traits<Graph>::vertex_descriptor> target_set;

            //rewire out edges
            typename graph_traits<Graph>::out_edge_iterator e, e_end;
            for (tie(e, e_end) = out_edges(v, g); e != e_end; ++e)
            {              
                typename graph_traits<Graph>::edge_descriptor swap_partner;

                //rewire edge
                swap_partner = rewire(*e, target_set, self_loops, parallel_edges);

                if (swap_partner == *e)
                    continue;

                target_set.insert(target(swap_partner, g));
                #pragma omp critical
                {
                    swap_list.push_back(make_pair(edge_index[*e], edge_index[swap_partner]));
                }
            }
        }
        
        // do the actual graph modification
        N = swap_list.size();
        #pragma omp parallel for default(shared) private(i) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            size_t e1, e2;
            tie(e1, e2) = swap_list[i];

            typename graph_traits<Graph>::vertex_descriptor s1, s2, new_target1, new_target2;
            typename graph_traits<Graph>::edge_descriptor  new_edge1, new_edge2 ;            
            
            s1 = source(edges[e1], g);
            new_target1 = target(edges[e2], g);
            s2 = source(edges[e2], g);
            new_target2 = target(edges[e1], g);

            #pragma omp critical
            {
                new_edge1 = add_edge(s1, new_target1, const_cast<Graph&>(g)).first;
                edge_index[new_edge1] = e1;
                new_edge2 = add_edge(s2, new_target2, const_cast<Graph&>(g)).first;
                edge_index[new_edge2] = e2;

                remove_edge(edges[e1], const_cast<Graph&>(g));
                edges[e1] = new_edge1;                    
                remove_edge(edges[e2], const_cast<Graph&>(g));
                edges[e2] = new_edge2;
            }
        }
    }
};

template <class Graph>
class RandomRewireStrategy
{
public:
    RandomRewireStrategy (const Graph& g, rng_t& rng): _g(g), _rng(rng)
    {
        int i, N = num_vertices(_g);
        #pragma omp parallel for default(shared) private(i) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, _g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            size_t k = get_in_degree()(v, _g);
            if (k == 0)
                continue;
            #pragma omp critical
            {
                if (!_vertices.empty())
                    _vertices[_vertices.rbegin()->first + k] = v;
                else
                    _vertices[k] = v;
            }
        }
    }
    
    typename graph_traits<Graph>::edge_descriptor operator()(const typename graph_traits<Graph>::edge_descriptor& e, 
                                                             tr1::unordered_set<typename graph_traits<Graph>::vertex_descriptor>& target_set,
                                                             bool self_loops, bool parallel_edges)
    {
        uniform_int<size_t> vertex_sample(0, _vertices.rbegin()->first-1);
        size_t v;
        typename graph_traits<Graph>::vertex_descriptor vertex;
        
        do
        {
            #pragma omp critical
            {
                v = vertex_sample(_rng);
            }
            vertex = _vertices.upper_bound(v)->second;
        }
        while ( (!self_loops && (vertex == source(e, _g)) ) ||
                (!parallel_edges && (target_set.find(vertex) != target_set.end())));

        vector<typename graph_traits<Graph>::edge_descriptor> candidates;
        typename get_in_edges<Graph>::iterator w, w_end;
        for (tie(w, w_end) = get_in_edges<Graph>()(vertex, _g); w != w_end; ++w)
        {
            if (parallel_edges || !is_adjacent(source(*w, _g), vertex, _g))
                candidates.push_back(*w);
        }

        if (candidates.empty())
            return e;

        uniform_int<size_t> edge_sample(0, candidates.size() - 1);
        size_t edge;
        #pragma omp critical
        {
           edge = edge_sample(_rng);
        }
        return candidates[edge];
    }

private:
    const Graph& _g;
    rng_t& _rng;
    typedef map<size_t, typename graph_traits<Graph>::vertex_descriptor> vertices_t;
    vertices_t _vertices;
};

template <class Graph>
class CorrelatedRewireStrategy
{
public:
    CorrelatedRewireStrategy (const Graph& g, rng_t& rng): _g(g), _rng(rng), _deg_sum(0)
    {
        int i, N = num_vertices(_g);
        #pragma omp parallel for default(shared) private(i) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {            
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, _g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            size_t j = get_in_degree()(v, _g);
            size_t k = out_degree(v, _g);

            #pragma omp critical
            {
                typename deg_vertices_t::value_type::second_type& vertices = _deg_vertices[make_pair(j,k)];
                if (!vertices.empty())
                    vertices[vertices.rbegin()->first + j] = v;
                else
                    vertices[j] = v;
            }
        }
    }

    typename graph_traits<Graph>::edge_descriptor operator()(const typename graph_traits<Graph>::edge_descriptor& e, 
                                                             tr1::unordered_set<typename graph_traits<Graph>::vertex_descriptor>& target_set,
                                                             bool self_loops, bool parallel_edges)
    {
        typename graph_traits<Graph>::vertex_descriptor t = target(e, _g); 
        typename deg_vertices_t::value_type::second_type& vertices = _deg_vertices[make_pair(get_in_degree()(t, _g), out_degree(t, _g))];

        uniform_int<size_t> vertex_sample(0, vertices.rbegin()->first - 1);
        size_t v;
        typename graph_traits<Graph>::vertex_descriptor vertex;
        
        do
        {
            #pragma omp critical
            {
                v = vertex_sample(_rng);
            }
            vertex = vertices.upper_bound(v)->second;
        }
        while ( (!self_loops && (vertex == source(e, _g)) ) ||
                (!parallel_edges && (target_set.find(vertex) != target_set.end())));

        vector<typename graph_traits<Graph>::edge_descriptor> candidates;
        typename get_in_edges<Graph>::iterator w, w_end;
        for (tie(w, w_end) = get_in_edges<Graph>()(vertex, _g); w != w_end; ++w)
        {
            if (parallel_edges || !is_adjacent(source(*w, _g), vertex, _g))
                candidates.push_back(*w);
        }

        if (candidates.empty())
            return e;

        uniform_int<size_t> edge_sample(0, candidates.size() - 1);
        size_t edge;
        #pragma omp critical
        {
            edge = edge_sample(_rng);
        }
        return candidates[edge];
    }

private:
    const Graph& _g;
    rng_t& _rng;
    typedef tr1::unordered_map<pair<size_t, size_t>, map<size_t, typename graph_traits<Graph>::vertex_descriptor>, hash<pair<size_t, size_t> > > deg_vertices_t;
    deg_vertices_t _deg_vertices;
    size_t _deg_sum;
};

//==============================================================================
// RandomRewire
//==============================================================================
void GraphInterface::RandomRewire(rewire_strat_t strat, bool self_loops, bool parallel_edges, size_t seed)
{
    bool reversed = GetReversed();
    SetReversed(false);

    if (strat == UNCORRELATED_STRAT)
        check_filter(*this, bind<void>(graph_rewire<RandomRewireStrategy>(), _1, _edge_index, seed, self_loops, parallel_edges),
                     never_reversed(), directed_check());
    else
        check_filter(*this, bind<void>(graph_rewire<CorrelatedRewireStrategy>(), _1, _edge_index, seed, self_loops, parallel_edges),
                     never_reversed(), directed_check());

    SetReversed(reversed);
}
