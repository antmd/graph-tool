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

template <class EdgeIterator>
typename iterator_traits<EdgeIterator>::value_type sample_edge(const EdgeIterator& e_begin, const EdgeIterator& e_end, rng_t& rng)
{
    vector<typename iterator_traits<EdgeIterator>::value_type> candidates;
    EdgeIterator e;
    for (e = e_begin; e != e_end; ++e)
        candidates.push_back(*e);

    uniform_int<size_t> edge_sample(0, candidates.size() - 1);
    size_t edge = edge_sample(rng);
    return candidates[edge];
}

template <class Graph, class EdgeIndexMap>
pair<typename graph_traits<Graph>::edge_descriptor,typename graph_traits<Graph>::edge_descriptor>
swap_edge_target(const typename graph_traits<Graph>::edge_descriptor& e1, const typename graph_traits<Graph>::edge_descriptor& e2, EdgeIndexMap edge_index, Graph& g)
{
    if (e1 == e2)
        return make_pair(e1, e2);
    typename graph_traits<Graph>::edge_descriptor ne1, ne2;
    ne1 = add_edge(source(e1, g), target(e2, g), g).first;
    edge_index[ne1] = edge_index[e1];
    ne2 = add_edge(source(e2, g), target(e1, g), g).first;
    edge_index[ne2] = edge_index[e2];
    remove_edge(e1, g);
    remove_edge(e2, g);
    return make_pair(ne1, ne2);
}

template <class Graph, class EdgeIndexMap>
pair<typename graph_traits<Graph>::edge_descriptor,typename graph_traits<Graph>::edge_descriptor>
swap_edge_source(const typename graph_traits<Graph>::edge_descriptor& e1, const typename graph_traits<Graph>::edge_descriptor& e2, EdgeIndexMap edge_index, Graph& g)
{
    if (e1 == e2)
        return make_pair(e1, e2);
    typename graph_traits<Graph>::edge_descriptor ne1, ne2;
    ne1 = add_edge(source(e2, g), target(e1, g), g).first;
    edge_index[ne1] = edge_index[e1];
    ne2 = add_edge(source(e1 , g), target(e2, g), g).first;
    edge_index[ne2] = edge_index[e2];
    remove_edge(e1, g);
    remove_edge(e2, g);
    return make_pair(ne1, ne2);
}


template <template <class Graph> class RewireStrategy>
struct graph_rewire
{
    template <class Graph, class EdgeIndexMap>
    void operator()(const Graph& og, EdgeIndexMap edge_index, size_t seed, bool self_loops, bool parallel_edges) const
    {
        Graph& g = const_cast<Graph&>(og);

        rng_t rng(static_cast<rng_t::result_type>(seed));

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

        RewireStrategy<Graph> rewire(g, rng);

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

        // for each edge simultaneously rewire its source and target
        for (i = 0; i < int(edges.size()); ++i)
        {
            typename graph_traits<Graph>::edge_descriptor e = edges[i];
            typename graph_traits<Graph>::vertex_descriptor s,t;

            tie(s, t) = rewire(e, self_loops);
            cout << source(e,g) << "->" << s << ",\t " << target(e,g) << "->" << t << endl;

            typename graph_traits<Graph>::edge_descriptor se = e, te = e;

            if (t != target(e, g))
            {
                typename get_in_edges<Graph>::iterator ie_begin, ie_end;
                tie(ie_begin, ie_end) = get_in_edges<Graph>()(t, g);
                if (ie_begin != ie_end)
                    te = sample_edge(ie_begin, ie_end, rng);
                if (te != e)
                {
                    assert(t == target(te,g));
                    cout << "(" << source(e,g) << "," << target(e,g) << ") -> (" << source(te,g) << "," << target(te,g) << ")"<< endl;
                    tie(e,te) = swap_edge_target(e, te, edge_index, g);
                    edges[edge_index[e]] = e;
                    edges[edge_index[te]] = te;
                }
            }

            if (s != source(e, g))
            {
                typename graph_traits<Graph>::out_edge_iterator oe, oe_begin, oe_end;
                tie(oe_begin, oe_end) = out_edges(s, g);
                for (oe = oe_begin; oe != oe_end; ++oe)
                    assert(s == source(*oe,g));
                if (oe_begin != oe_end)
                    se = sample_edge(oe_begin, oe_end, rng);
                if (se != e)
                {
                    assert(s == source(se,g));
                    cout << "(" << source(e,g) << "," << target(e,g) << ") -> (" << source(se,g) << "," << target(se,g) << ")"<< endl;
                    tie(e, se) = swap_edge_source(e, se, edge_index, g);
                    edges[edge_index[e]] = e;
                    edges[edge_index[se]] = se;
                }
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
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, _g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            size_t k = get_in_degree()(v, _g);
            if (k == 0)
                continue;
            if (!_vertices.empty())
                _vertices[_vertices.rbegin()->first + k] = v;
            else
                _vertices[k] = v;
        }
    }

    pair<typename graph_traits<Graph>::vertex_descriptor, typename graph_traits<Graph>::vertex_descriptor>
    operator()(const typename graph_traits<Graph>::edge_descriptor& e, bool self_loops)
    {
        uniform_int<size_t> vertex_sample(0, _vertices.rbegin()->first-1);
        size_t v1,v2;
        typename graph_traits<Graph>::vertex_descriptor s,t;
        do
        {
            v1 = vertex_sample(_rng);
            v2 = vertex_sample(_rng);
            s = _vertices.upper_bound(v1)->second;
            t = _vertices.upper_bound(v2)->second;
        }
        while (!self_loops && (s == t));
        return make_pair(s, t);
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
    CorrelatedRewireStrategy (const Graph& g, rng_t& rng): _g(g), _rng(rng)
    {
        int i, N = num_vertices(_g);
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, _g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            size_t j = get_in_degree()(v, _g);
            size_t k = out_degree(v, _g);

            if (j > 0)
            {
                typename deg_vertices_t::value_type::second_type& target_vertices = _deg_target_vertices[make_pair(j,k)];
                if (!target_vertices.empty())
                    target_vertices[target_vertices.rbegin()->first + j] = v;
                else
                    target_vertices[j] = v;
            }
            
            if (k > 0)
            {
                typename deg_vertices_t::value_type::second_type& source_vertices = _deg_source_vertices[make_pair(j,k)];
                if (!source_vertices.empty())
                    source_vertices[source_vertices.rbegin()->first + k] = v;
                else
                    source_vertices[k] = v;
            }
        }
    }

    pair<typename graph_traits<Graph>::vertex_descriptor, typename graph_traits<Graph>::vertex_descriptor>
    operator()(const typename graph_traits<Graph>::edge_descriptor& e, bool self_loops)
    {
        typename deg_vertices_t::value_type::second_type& source_vertices = _deg_source_vertices[make_pair(get_in_degree()(source(e, _g), _g), out_degree(source(e, _g), _g))];
        typename deg_vertices_t::value_type::second_type& target_vertices = _deg_target_vertices[make_pair(get_in_degree()(target(e, _g), _g), out_degree(target(e, _g), _g))];

        //sample the new target
        uniform_int<size_t> source_sample(0, source_vertices.rbegin()->first - 1);
        uniform_int<size_t> target_sample(0, target_vertices.rbegin()->first - 1);

        size_t v1,v2;
        typename graph_traits<Graph>::vertex_descriptor s,t;
        do
        {
            #pragma omp critical
            {
                v1 = source_sample(_rng);
                v2 = target_sample(_rng);
            }
            s = source_vertices.upper_bound(v1)->second;
            t = target_vertices.upper_bound(v2)->second;
        }
        while ( !self_loops && (s == t) );

        return make_pair(s, t);
    }

private:
    const Graph& _g;
    rng_t& _rng;
    typedef tr1::unordered_map<pair<size_t, size_t>, map<size_t, typename graph_traits<Graph>::vertex_descriptor>, hash<pair<size_t, size_t> > > deg_vertices_t;
    deg_vertices_t _deg_source_vertices;
    deg_vertices_t _deg_target_vertices;
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
