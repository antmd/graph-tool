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

#ifndef GRAPH_REWIRING_HH
#define GRAPH_REWIRING_HH

#include <tr1/unordered_set>
#include <tr1/random>
#include <boost/functional/hash.hpp>
#include <iostream>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/multi_index/identity.hpp>

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_util.hh"

namespace graph_tool
{
using namespace std;
using namespace boost;
using namespace boost::multi_index;

// returns true if vertices u and v are adjacent. This is O(k(u)).
template <class Graph>
bool is_adjacent(typename graph_traits<Graph>::vertex_descriptor u,
                 typename graph_traits<Graph>::vertex_descriptor v,
                 const Graph& g )
{
    typename graph_traits<Graph>::out_edge_iterator e, e_end;
    for (tie(e, e_end) = out_edges(u, g); e != e_end; ++e)
    {
        if (target(*e,g) == v)
            return true;
    }
    return false;
}

template <class Graph>
typename graph_traits<Graph>::vertex_descriptor
source(const pair<typename graph_traits<Graph>::edge_descriptor, bool>& e,
       const Graph& g)
{
    if (e.second)
        return target(e.first, g);
    else
        return source(e.first, g);
}

template <class Graph>
typename graph_traits<Graph>::vertex_descriptor
target(const pair<typename graph_traits<Graph>::edge_descriptor, bool>& e,
       const Graph& g)
{
    if (e.second)
        return source(e.first, g);
    else
        return target(e.first, g);
}

// This will iterate over a random permutation of a random access sequence, by
// swapping the values of the sequence as it iterates
template <class RandomAccessIterator, class RNG,
          class RandomDist = tr1::uniform_int<size_t> >
class random_permutation_iterator : public
    std::iterator<input_iterator_tag, typename RandomAccessIterator::value_type>
{
public:
    random_permutation_iterator(RandomAccessIterator begin,
                                RandomAccessIterator end, RNG& rng)
        : _i(begin), _end(end), _rng(&rng)
    {
        if(_i != _end)
        {
            RandomDist random(0,  _end - _i - 1);
            std::iter_swap(_i, _i + random(*_rng));
        }
    }

    typename RandomAccessIterator::value_type operator*()
    {
        return *_i;
    }

    random_permutation_iterator& operator++()
    {
        ++_i;
        if(_i != _end)
        {
            RandomDist random(0,  _end - _i - 1);
            std::iter_swap(_i, _i + random(*_rng));
        }
        return *this;
    }

    bool operator==(const random_permutation_iterator& ri)
    {
        return _i == ri._i;
    }

    bool operator!=(const random_permutation_iterator& ri)
    {
        return _i != ri._i;
    }

    size_t operator-(const random_permutation_iterator& ri)
    {
        return _i - ri._i;
    }

private:
    RandomAccessIterator _i, _end;
    RNG* _rng;
};

// this functor will swap the source of the edge e with the source of edge se
// and the target of edge e with the target of te
struct swap_edge
{
    template <class Graph, class NewEdgeMap>
    static bool
    parallel_check_source
    (const typename graph_traits<Graph>::edge_descriptor& e,
     const pair<typename graph_traits<Graph>::edge_descriptor, bool>& se,
     NewEdgeMap edge_is_new, const Graph &g)
    {
        // We want to check that if we swap the source of 'e' with the source of
        // 'se', as such
        //
        //  (s)    -e--> (t)          (ns)   -e--> (t)
        //  (ns)   -se-> (se_t)   =>  (s)    -se-> (se_t)
        //
        // no parallel edges are introduced. We must considered only "new
        // edges", i.e., edges which were already sampled and swapped. "Old
        // edges" will have their chance of being swapped, and then they'll be
        // checked for parallelism.

        typename graph_traits<Graph>::vertex_descriptor
            s = source(e, g),          // current source
            t = target(e, g),          // current target
            ns = source(se, g),        // new source
            se_t = target(se, g);      // source edge target

        if (is_adjacent_in_new(ns,  t, edge_is_new, g))
            return true; // e would clash with an existing (new) edge
        if (edge_is_new[se.first] && is_adjacent_in_new(s, se_t, edge_is_new,g))
            return true; // se would clash with an existing (new) edge
        return false; // the coast is clear - hooray!
    }

    template <class Graph, class NewEdgeMap>
    static bool
    parallel_check_target
    (const typename graph_traits<Graph>::edge_descriptor& e,
     const pair<typename graph_traits<Graph>::edge_descriptor, bool>& te,
     NewEdgeMap edge_is_new, const Graph &g)
    {
        // We want to check that if we swap the target of 'e' with the target of
        // 'te', as such
        //
        //  (s)    -e--> (t)          (s)    -e--> (nt)
        //  (te_s) -te-> (nt)   =>    (te_s) -te-> (t)
        //
        // no parallel edges are introduced. We must considered only "new
        // edges", i.e., edges which were already sampled and swapped. "Old
        // edges" will have their chance of being swapped, and then they'll be
        // checked for parallelism.

        typename graph_traits<Graph>::vertex_descriptor
            s = source(e, g),          // current source
            t = target(e, g),          // current target
            nt = target(te, g),        // new target
            te_s = source(te, g);      // target edge source

        if (is_adjacent_in_new(s,  nt, edge_is_new, g))
            return true; // e would clash with an existing (new) edge
        if (edge_is_new[te.first] && is_adjacent_in_new(te_s, t, edge_is_new,g))
            return true; // te would clash with an existing (new) edge
        return false; // the coast is clear - hooray!
    }

    // returns true if vertices u and v are adjacent in the new graph. This is
    // O(k(u)).
    template <class Graph, class EdgeIsNew>
    static bool is_adjacent_in_new
        (typename graph_traits<Graph>::vertex_descriptor u,
         typename graph_traits<Graph>::vertex_descriptor v,
         EdgeIsNew edge_is_new, const Graph& g)
    {
        typename graph_traits<Graph>::out_edge_iterator e, e_end;
        for (tie(e, e_end) = out_edges(u, g); e != e_end; ++e)
        {
            if (edge_is_new[*e] && target(*e,g) == v)
                return true;
        }
        return false;
    }

    template <class Graph, class EdgeIndexMap, class EdgesType>
    void operator()(const typename graph_traits<Graph>::edge_descriptor& e,
                    const pair<typename graph_traits<Graph>::edge_descriptor,
                               bool>& se,
                    const pair<typename graph_traits<Graph>::edge_descriptor,
                               bool>& te,
                    EdgesType& edges, EdgeIndexMap edge_index, Graph& g)
    {
        // swap the source of the edge 'e' with the source of edge 'se' and the
        // target of edge 'e' with the target of 'te', as such:
        //
        //  (s)    -e--> (t)          (ns)   -e--> (nt)
        //  (ns)   -se-> (se_t)   =>  (s)    -se-> (se_t)
        //  (te_s) -te-> (nt)         (te_s) -te-> (t),

        // new edges which will replace the old ones
        typename graph_traits<Graph>::edge_descriptor ne, nse, nte;

        // split cases where different combinations of the three edges are
        // the same
        if(se.first != te.first)
        {
            ne = add_edge(source(se, g), target(te, g), g).first;
            if(e != se.first)
            {
                if (!se.second)
                    nse = add_edge(source(e, g), target(se, g), g).first;
                else // keep invertedness (only for undirected graphs)
                    nse = add_edge(target(se, g), source(e, g), g).first;
                edge_index[nse] = edge_index[se.first];
                remove_edge(se.first, g);
                edges[edge_index[nse]] = nse;
            }
            if(e != te.first)
            {
                if (!te.second)
                    nte = add_edge(source(te, g), target(e, g), g).first;
                else // keep invertedness (only for undirected graphs)
                    nte = add_edge(target(e, g), source(te, g), g).first;
                edge_index[nte] = edge_index[te.first];
                remove_edge(te.first, g);
                edges[edge_index[nte]] = nte;
            }
            edge_index[ne] = edge_index[e];
            remove_edge(e, g);
            edges[edge_index[ne]] = ne;
        }
        else
        {
            if(e != se.first)
            {
                // se and te are the same. swapping indexes only.
                swap(edge_index[se.first], edge_index[e]);
                edges[edge_index[se.first]] = se.first;
                edges[edge_index[e]] = e;
            }
        }
    }
};

// used for verbose display
void print_progress(size_t current, size_t total, stringstream& str);

// main rewire loop
template <template <class Graph, class EdgeIndexMap> class RewireStrategy>
struct graph_rewire
{
    template <class Graph, class EdgeIndexMap>
    void operator()(Graph& g, EdgeIndexMap edge_index, rng_t& rng,
                    bool self_loops, bool parallel_edges) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;

        if (!self_loops)
        {
            // check the existence of self-loops
            bool has_self_loops = false;
            int i, N = num_vertices(g);
            #pragma omp parallel for default(shared) private(i) \
                schedule(dynamic)
            for (i = 0; i < N; ++i)
            {
                vertex_t v = vertex(i, g);
                if (v == graph_traits<Graph>::null_vertex())
                    continue;
                if (is_adjacent(v, v, g))
                    has_self_loops = true;
            }
            if (has_self_loops)
                throw ValueException("Self-loop detected. Can't rewire graph "
                                     "without self-loops if it already contains"
                                     " self-loops!");
        }

        if (!parallel_edges)
        {
            // check the existence of parallel edges
            bool has_parallel_edges = false;
            int i, N = num_vertices(g);
            #pragma omp parallel for default(shared) private(i) \
                schedule(dynamic)
            for (i = 0; i < N; ++i)
            {
                vertex_t v = vertex(i, g);
                if (v == graph_traits<Graph>::null_vertex())
                    continue;

                tr1::unordered_set<vertex_t> targets;
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
                throw ValueException("Parallel edge detected. Can't rewire "
                                     "graph without parallel edges if it "
                                     "already contains parallel edges!");
        }

        vector<edge_t> edges(num_edges(g));
        vector<size_t> edges_sequence;
        typename graph_traits<Graph>::edge_iterator e, e_end;
        for (tie(e, e_end) = boost::edges(g); e != e_end; ++e)
        {
            if (edge_index[*e] >= edges.size())
                edges.resize(edge_index[*e] + 1);
            edges[edge_index[*e]] = *e;
            edges_sequence.push_back(edge_index[*e]);
        }

        RewireStrategy<Graph, EdgeIndexMap> rewire(g, edge_index, edges, rng);

        typedef random_permutation_iterator<typename vector<size_t>::iterator,
                                            rng_t>
            random_edge_iter;

        random_edge_iter
            ei(edges_sequence.begin(), edges_sequence.end(), rng),
            ei_end(edges_sequence.end(), edges_sequence.end(), rng);

        stringstream str;
        //size_t count = 0;
        // for each edge simultaneously rewire its source or target
        for (; ei != ei_end; ++ei)
        {
            //print_progress(count++, edges_sequence.size(), str);
            rewire(edges[*ei], self_loops, parallel_edges);
        }
    }
};


// this will rewire the edges so that the resulting graph will be entirely
// random (i.e. Erdos-Renyi)
template <class Graph, class EdgeIndexMap>
class ErdosRewireStrategy
{
public:
    typedef typename graph_traits<Graph>::edge_descriptor edge_t;
    typedef typename EdgeIndexMap::value_type index_t;

    ErdosRewireStrategy(Graph& g, EdgeIndexMap edge_index,
                        vector<edge_t>& edges, rng_t& rng)
        : _g(g), _edge_index(edge_index), _edges(edges),
          _vertices(HardNumVertices()(g)), _rng(rng)
    {
        typeof(_vertices.begin()) viter = _vertices.begin();
        typename graph_traits<Graph>::vertex_iterator v, v_end;
        for (tie(v, v_end) = vertices(_g); v != v_end; ++v)
            *(viter++) = *v;
    }

    void operator()(const edge_t& e, bool self_loops, bool parallel_edges)
    {
        //try randomly drawn pairs of vertices until one satisfies all the
        //consistency checks
        typedef random_permutation_iterator
            <typename graph_traits<Graph>::vertex_iterator, rng_t>
            random_vertex_iter;

        tr1::uniform_int<size_t> sample(0, _vertices.size());
        typename graph_traits<Graph>::vertex_descriptor s, t;
        while (true)
        {
            s = sample(_rng);
            t = sample(_rng);

            if(s == t && !self_loops) // reject self-loops if not allowed
                continue;
            if (!parallel_edges &&
                swap_edge::is_adjacent_in_new(s, t, _edge_is_new, _g))
                continue;  // reject parallel edges if not allowed
            break;
        }
        edge_t ne = add_edge(s, t, _g).first;
        _edges[_edge_index[e]] = ne;
        remove_edge(e, _g);
        if (_edge_index[ne] >= _edges.size())
            _edges.resize(_edge_index[ne] + 1);
        _edges[_edge_index[ne]] = ne;
        _edge_is_new[ne] = true;
    }

private:
    Graph& _g;
    EdgeIndexMap _edge_index;
    vector<edge_t>& _edges;
    vector<typename graph_traits<Graph>::vertex_descriptor> _vertices;
    checked_vector_property_map<bool, EdgeIndexMap> _edge_is_new;
    rng_t& _rng;
};


// this is the mother class for edge-based rewire strategies
// it contains the common loop for finding edges to swap, so different
// strategies need only to specify where to sample the edges from.
template <class Graph, class EdgeIndexMap, class RewireStrategy>
class RewireStrategyBase
{
public:
    typedef typename graph_traits<Graph>::edge_descriptor edge_t;
    typedef typename graph_traits<Graph>::edge_descriptor vertex_t;

    typedef typename EdgeIndexMap::value_type index_t;

    RewireStrategyBase(Graph& g, EdgeIndexMap edge_index, vector<edge_t>& edges,
                       rng_t& rng)
        : _g(g), _edge_index(edge_index), _edges(edges),
          _edge_is_new(edge_index, 1), _rng(rng)
    {
        typename graph_traits<Graph>::edge_iterator e_i, e_i_end;
        for (tie(e_i, e_i_end) = boost::edges(g); e_i != e_i_end; ++e_i)
            _edge_is_new.get_checked()[*e_i] == true;
    }

    void operator()(const edge_t& e, bool self_loops, bool parallel_edges)
    {
        RewireStrategy& self = *static_cast<RewireStrategy*>(this);
        tr1::bernoulli_distribution coin;

        //try randomly drawn pairs of edges until one satisfies all the
        //consistency checks
        pair<edge_t,bool> es, et;

        _edge_is_new[e] = false;
        while (true)
        {
            es.first = et.first = e;
            es.second = et.second = false;

            if (coin(_rng))
            {
                // rewire source
                es = self.get_source_edge(e);

                if(!self_loops) // reject self-loops if not allowed
                {
                    if((source(e, _g) == target(es, _g)) ||
                       (target(e, _g) == source(es, _g)))
                        continue;
                }
                if (!parallel_edges) // reject parallel edges if not allowed
                {
                    if (swap_edge::parallel_check_source(e, es, _edge_is_new,
                                                         _g))
                        continue;
                }
            }
            else
            {
                // rewire target
                et = self.get_target_edge(e);

                if (!self_loops) // reject self-loops if not allowed
                {
                    if((source(e, _g) == target(et, _g)) ||
                       (target(e, _g) == source(et, _g)))
                        continue;
                }
                if (!parallel_edges) // reject parallel edges if not allowed
                {
                    if (swap_edge::parallel_check_target(e, et, _edge_is_new,
                                                         _g))
                        continue;
                }
            }
            break;
        }

        swap_edge()(e, es, et, _edges, _edge_index, _g);
        _edge_is_new[e] = true;

        const edge_t* elist[3] = {&e, &es.first, &et.first};
        for (size_t i = 0; i < 3; ++i)
        {
            self.update_edge(*elist[i], false);
            self.update_edge(*elist[i], true);
        }

    }

protected:
    Graph& _g;
    EdgeIndexMap _edge_index;
    vector<edge_t>& _edges;
    unchecked_vector_property_map<bool, EdgeIndexMap> _edge_is_new;
    rng_t& _rng;
};

// this will rewire the edges so that the combined (in, out) degree distribution
// will be the same, but all the rest is random
template <class Graph, class EdgeIndexMap>
class RandomRewireStrategy:
    public RewireStrategyBase<Graph, EdgeIndexMap,
                              RandomRewireStrategy<Graph, EdgeIndexMap> >
{
public:
    typedef RewireStrategyBase<Graph, EdgeIndexMap,
                               RandomRewireStrategy<Graph, EdgeIndexMap> >
        base_t;

    typedef Graph graph_t;
    typedef EdgeIndexMap edge_index_t;

    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    typedef typename graph_traits<Graph>::edge_descriptor edge_t;
    typedef typename EdgeIndexMap::value_type index_t;

    struct hash_index {};
    struct random_index {};

    typedef multi_index_container<
        pair<index_t,bool>,
        indexed_by<
            hashed_unique<tag<hash_index>, identity<pair<index_t,bool> > >,
            random_access<tag<random_index> > > >
        edge_set_t;

    RandomRewireStrategy(Graph& g, EdgeIndexMap edge_index,
                         vector<edge_t>& edges, rng_t& rng)
        : base_t(g, edge_index, edges, rng), _g(g)
    {
        size_t E = HardNumEdges()(_g);

        typename graph_traits<Graph>::vertex_iterator v_i, v_i_end;
        for (tie(v_i, v_i_end) = vertices(g); v_i != v_i_end; ++v_i)
        {
            if (out_degree(*v_i, _g) > E/3)
            {
                _hub_non_adjacent[*v_i] = edge_set_t();
                if(!is_directed::apply<Graph>::type::value)
                    _hub_non_incident[*v_i] = edge_set_t();
            }

            if (in_degreeS()(*v_i, _g) > E/3)
                _hub_non_incident[*v_i] = edge_set_t();
        }

        typename graph_traits<Graph>::edge_iterator e_i, e_i_end;
        for (tie(e_i, e_i_end) = boost::edges(g); e_i != e_i_end; ++e_i)
        {
            // For undirected graphs, there is no difference between source and
            // target, and each edge will appear _twice_ in the list below,
            // once for each different ordering of source and target.

            _all_edges.push_back(make_pair(edge_index[*e_i], false));
            if (!is_directed::apply<Graph>::type::value)
                _all_edges.push_back(make_pair(edge_index[*e_i], true));
            update_edge(*e_i, true);
        }

        for (typeof(_hub_non_adjacent.begin()) iter = _hub_non_adjacent.begin();
             iter != _hub_non_adjacent.end(); ++iter)
            cout << iter->first << " " << iter->second.size() << endl;
        for (typeof(_hub_non_incident.begin()) iter = _hub_non_incident.begin();
             iter != _hub_non_incident.end(); ++iter)
            cout << iter->first << " " << iter->second.size() << endl;
    }

    pair<edge_t,bool> get_edge(const edge_t& e, bool src)
    {
        tr1::unordered_map<size_t, edge_set_t>& hub =
            (src || !is_directed::apply<Graph>::type::value) ?
            _hub_non_adjacent : _hub_non_incident;

        typeof(hub.begin()) iter =
            hub.find(src ? target(e,_g) : source(e,_g));
        if (iter != hub.end())
        {
            edge_set_t& eset = iter->second;

            tr1::uniform_int<> sample(0, eset.size());
            size_t i = sample(base_t::_rng);
            if (i == eset.size()) // no rewire option
            {
                return make_pair(e, false);
            }
            else
            {
                pair<index_t,bool> ep =
                    *(eset.get<random_index>().begin() + i);
                return  make_pair(base_t::_edges[ep.first], ep.second);
            }
        }

        tr1::uniform_int<> sample(0, _all_edges.size()-1);
        pair<index_t,bool> ep = _all_edges[sample(base_t::_rng)];
        return make_pair(base_t::_edges[ep.first], ep.second);
    }


    pair<edge_t,bool> get_source_edge(const edge_t& e)
    {
        return get_edge(e, true);
    }

    pair<edge_t,bool> get_target_edge(const edge_t& e)
    {
        return get_edge(e, false);
    }

    void update_edge_list(index_t v, const pair<edge_t, bool>& e,
                          edge_set_t& edge_set, bool src, bool insert)
    {
        pair<index_t,bool> ep;
        if (!base_t::_edge_is_new[e.first] || !insert ||
            (src && source(e.first, _g) != v) ||
            (!src && target(e.first, _g) != v))
        {
            ep = make_pair(_edge_index[e.first], e.second);
            if (insert)
                edge_set.get<hash_index>().insert(ep);
            else
                edge_set.get<hash_index>().erase(ep);
        }
    }

    void update_edge(const edge_t& e, bool insert)
    {
        for (typeof(_hub_non_adjacent.begin()) iter = _hub_non_adjacent.begin();
             iter != _hub_non_adjacent.end(); ++iter)
        {
            update_edge_list(iter->first, make_pair(e, false), iter->second,
                             true, insert);
            if (!is_directed::apply<Graph>::type::value)
                update_edge_list(iter->first, make_pair(e, true), iter->second,
                                 true, insert);
        }

        for (typeof(_hub_non_incident.begin()) iter = _hub_non_incident.begin();
             iter != _hub_non_incident.end(); ++iter)
        {
            update_edge_list(iter->first, make_pair(e, false), iter->second,
                             false, insert);
            if (!is_directed::apply<Graph>::type::value)
                update_edge_list(iter->first, make_pair(e, true), iter->second,
                                 false, insert);
        }
    }

private:
    Graph& _g;
    EdgeIndexMap _edge_index;
    vector<pair<index_t,bool> > _all_edges;

    tr1::unordered_map<size_t, edge_set_t>  _hub_non_incident;
    tr1::unordered_map<size_t, edge_set_t>  _hub_non_adjacent;
};


// this will rewire the edges so that the (in,out) degree distributions and the
// (in,out)->(in,out) correlations will be the same, but all the rest is random
template <class Graph, class EdgeIndexMap>
class CorrelatedRewireStrategy:
    public RewireStrategyBase<Graph, EdgeIndexMap,
                              CorrelatedRewireStrategy<Graph, EdgeIndexMap> >
{
public:
    typedef RewireStrategyBase<Graph, EdgeIndexMap,
                               CorrelatedRewireStrategy<Graph, EdgeIndexMap> >
        base_t;

    typedef Graph graph_t;
    typedef EdgeIndexMap edge_index_t;

    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    typedef typename graph_traits<Graph>::edge_descriptor edge_t;
    typedef typename EdgeIndexMap::value_type index_t;

    CorrelatedRewireStrategy (Graph& g, EdgeIndexMap edge_index,
                              vector<edge_t>& edges, rng_t& rng)
        : base_t(g, edge_index, edges, rng), _g(g)
    {
        typename graph_traits<Graph>::edge_iterator e_i, e_i_end;
        for (tie(e_i, e_i_end) = boost::edges(g); e_i != e_i_end; ++e_i)
        {
            // For undirected graphs, there is no difference between source and
            // target, and each edge will appear _twice_ in the lists below,
            // once for each different ordering of source and target.

            _edges_by_source
                [make_pair(in_degreeS()(source(*e_i, _g), _g),
                           out_degree(source(*e_i, _g), _g))]
                .push_back(make_pair(edge_index[*e_i], false));

            _edges_by_target
                [make_pair(in_degreeS()(target(*e_i, _g), _g),
                           out_degree(target(*e_i, _g), _g))]
                .push_back(make_pair(edge_index[*e_i], false));

            if (!is_directed::apply<Graph>::type::value)
            {
                _edges_by_source
                    [make_pair(in_degreeS()(target(*e_i, _g), _g),
                               out_degree(target(*e_i, _g), _g))]
                    .push_back(make_pair(edge_index[*e_i], true));
                _edges_by_target
                    [make_pair(in_degreeS()(source(*e_i, _g), _g),
                               out_degree(source(*e_i, _g), _g))]
                    .push_back(make_pair(edge_index[*e_i], true));
            }
        }
    }

    pair<edge_t,bool> get_edge(const edge_t& e, bool src)
    {
        pair<size_t, size_t> deg =
            make_pair(in_degreeS()(src ? source(e, _g) : target(e, _g), _g),
                      out_degree(src ? source(e, _g) : target(e, _g), _g));
        edges_by_end_deg_t& edges = src ? _edges_by_source : _edges_by_target;
        typename edges_by_end_deg_t::mapped_type& elist = edges[deg];
        tr1::uniform_int<> sample(0, elist.size()-1);
        pair<index_t,bool> ep = elist[sample(base_t::_rng)];
        return make_pair(base_t::_edges[ep.first], ep.second);
    }

    pair<edge_t,bool> get_source_edge(const edge_t& e)
    {
        return get_edge(e, true);
    }

    pair<edge_t,bool> get_target_edge(const edge_t& e)
    {
        return get_edge(e, false);
    }

    void update_edge(const edge_t& e, bool insert) {}

private:
    typedef tr1::unordered_map<pair<size_t, size_t>,
                               vector<pair<index_t, bool> >,
                               hash<pair<size_t, size_t> > > edges_by_end_deg_t;
    edges_by_end_deg_t _edges_by_source, _edges_by_target;

protected:
    const Graph& _g;
};

} // graph_tool namespace

#endif // GRAPH_REWIRING_HH
