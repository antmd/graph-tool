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
#include "sampler.hh"

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
source(const pair<size_t, bool>& e,
       const vector<typename graph_traits<Graph>::edge_descriptor>& edges,
       const Graph& g)
{
    if (e.second)
        return target(edges[e.first], g);
    else
        return source(edges[e.first], g);
}

template <class Graph>
typename graph_traits<Graph>::vertex_descriptor
target(const pair<size_t, bool>& e,
       const vector<typename graph_traits<Graph>::edge_descriptor>& edges,
       const Graph& g)
{
    if (e.second)
        return source(edges[e.first], g);
    else
        return target(edges[e.first], g);
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
    template <class Graph>
    static bool
    parallel_check_source (size_t e, pair<size_t, bool> se,
                           vector<typename graph_traits<Graph>::edge_descriptor>& edges,
                           const Graph &g)
    {
        // We want to check that if we swap the source of 'e' with the source of
        // 'se', as such
        //
        //  (s)    -e--> (t)          (ns)   -e--> (t)
        //  (ns)   -se-> (se_t)   =>  (s)    -se-> (se_t)
        //
        // no parallel edges are introduced.

        typename graph_traits<Graph>::vertex_descriptor
            s = source(edges[e], g),          // current source
            t = target(edges[e], g),          // current target
            ns = source(se, edges, g),        // new source
            se_t = target(se, edges, g);      // source edge target

        if (is_adjacent(ns,  t, g))
            return true; // e would clash with an existing (new) edge
        if (is_adjacent(s, se_t, g))
            return true; // se would clash with an existing (new) edge
        return false; // the coast is clear - hooray!
    }

    template <class Graph>
    static bool
    parallel_check_target (size_t e, const pair<size_t, bool>& te,
                           vector<typename graph_traits<Graph>::edge_descriptor>& edges,
                           const Graph &g)
    {
        // We want to check that if we swap the target of 'e' with the target of
        // 'te', as such
        //
        //  (s)    -e--> (t)          (s)    -e--> (nt)
        //  (te_s) -te-> (nt)   =>    (te_s) -te-> (t)
        //
        // no parallel edges are introduced.

        typename graph_traits<Graph>::vertex_descriptor
            s = source(edges[e], g),          // current source
            t = target(edges[e], g),          // current target
            nt = target(te, edges, g),        // new target
            te_s = source(te, edges, g);      // target edge source

        if (is_adjacent(s,  nt, g))
            return true; // e would clash with an existing (new) edge
        if (is_adjacent(te_s, t, g))
            return true; // te would clash with an existing (new) edge
        return false; // the coast is clear - hooray!
    }

    template <class Graph, class EdgeIndexMap>
    static void swap_source (size_t e, pair<size_t, bool>& se,
                             vector<typename graph_traits<Graph>::edge_descriptor>& edges,
                             EdgeIndexMap edge_index, Graph& g)
    {
        // swap the source of the edge 'e' with the source of edge 'se', as
        // such:
        //
        //  (s)    -e--> (t)          (ns)   -e--> (t)
        //  (ns)   -se-> (se_t)   =>  (s)    -se-> (se_t)

        if (e == se.first)
            return;

        // new edges which will replace the old ones
        typename graph_traits<Graph>::edge_descriptor ne, nse;
        ne = add_edge(source(se, edges, g), target(edges[e], g), g).first;
        if (!se.second)
            nse = add_edge(source(edges[e], g), target(se, edges, g), g).first;
        else // keep invertedness (only for undirected graphs)
            nse = add_edge(target(se, edges, g), source(edges[e], g),  g).first;
        edge_index[nse] = se.first;
        remove_edge(edges[se.first], g);
        edges[se.first] = nse;

        edge_index[ne] = e;
        remove_edge(edges[e], g);
        edges[e] = ne;
    }

    template <class Graph, class EdgeIndexMap>
    static void swap_target
        (size_t e, const pair<size_t, bool>& te,
         vector<typename graph_traits<Graph>::edge_descriptor>& edges,
         EdgeIndexMap edge_index, Graph& g)
    {
        // swap the source of the edge 'e' with the source of edge 'se', as
        // such:
        //
        //  (s)    -e--> (t)          (s)    -e--> (nt)
        //  (te_s) -te-> (nt)   =>    (te_s) -te-> (t)

        if (e == te.first)
            return;

        // new edges which will replace the old ones
        typename graph_traits<Graph>::edge_descriptor ne, nte;
        ne = add_edge(source(edges[e], g), target(te, edges, g), g).first;
        if (!te.second)
            nte = add_edge(source(te, edges, g), target(edges[e], g), g).first;
        else // keep invertedness (only for undirected graphs)
            nte = add_edge(target(edges[e], g), source(te, edges, g),  g).first;

        edge_index[nte] = te.first;
        remove_edge(edges[te.first], g);
        edges[te.first] = nte;

        edge_index[ne] = e;
        remove_edge(edges[e], g);
        edges[e] = ne;
    }

};

// used for verbose display
void print_progress(size_t current, size_t total, stringstream& str);

// main rewire loop
template <template <class Graph, class EdgeIndexMap, class CorrProb>
          class RewireStrategy>
struct graph_rewire
{
    template <class Graph, class EdgeIndexMap, class CorrProb>
    void operator()(Graph& g, EdgeIndexMap edge_index, CorrProb corr_prob,
                    rng_t& rng, bool self_loops, bool parallel_edges,
                    bool verbose) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;

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

        RewireStrategy<Graph, EdgeIndexMap, CorrProb>
            rewire(g, edge_index, edges, corr_prob, rng);

        typedef random_permutation_iterator<typename vector<size_t>::iterator,
                                            rng_t>
            random_edge_iter;

        random_edge_iter
            ei(edges_sequence.begin(), edges_sequence.end(), rng),
            ei_end(edges_sequence.end(), edges_sequence.end(), rng);

        if (verbose)
            cout << "rewiring edges: ";

        stringstream str;
        size_t count = 0;
        // for each edge rewire its source or target
        for (; ei != ei_end; ++ei)
        {
            if (verbose)
                print_progress(count++, edges_sequence.size(), str);
            rewire(*ei, self_loops, parallel_edges);
        }
        if (verbose)
            cout << endl;
    }
};


// this will rewire the edges so that the resulting graph will be entirely
// random (i.e. Erdos-Renyi)
template <class Graph, class EdgeIndexMap, class CorrProb>
class ErdosRewireStrategy
{
public:
    typedef typename graph_traits<Graph>::edge_descriptor edge_t;
    typedef typename EdgeIndexMap::value_type index_t;

    ErdosRewireStrategy(Graph& g, EdgeIndexMap edge_index,
                        vector<edge_t>& edges, CorrProb corr, rng_t& rng)
        : _g(g), _edge_index(edge_index), _edges(edges),
          _vertices(HardNumVertices()(g)), _rng(rng)
    {
        typeof(_vertices.begin()) viter = _vertices.begin();
        typename graph_traits<Graph>::vertex_iterator v, v_end;
        for (tie(v, v_end) = vertices(_g); v != v_end; ++v)
            *(viter++) = *v;
    }

    void operator()(size_t e, bool self_loops, bool parallel_edges)
    {
        //try randomly drawn pairs of vertices until one satisfies all the
        //consistency checks
        typedef random_permutation_iterator
            <typename graph_traits<Graph>::vertex_iterator, rng_t>
            random_vertex_iter;

        tr1::uniform_int<size_t> sample(0, _vertices.size()-1);
        typename graph_traits<Graph>::vertex_descriptor s, t;
        while (true)
        {
            s = sample(_rng);
            t = sample(_rng);

            if(s == t && !self_loops) // reject self-loops if not allowed
                continue;
            if (!parallel_edges && is_adjacent(s, t, _g))
                continue;  // reject parallel edges if not allowed
            break;
        }
        edge_t ne = add_edge(s, t, _g).first;
        _edge_index[ne] = e;
        remove_edge(_edges[e], _g);
        _edges[e] = ne;
    }

private:
    Graph& _g;
    EdgeIndexMap _edge_index;
    vector<edge_t>& _edges;
    vector<typename graph_traits<Graph>::vertex_descriptor> _vertices;
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
        : _g(g), _edge_index(edge_index), _edges(edges),  _rng(rng) {}

    void operator()(size_t e, bool self_loops, bool parallel_edges)
    {
        RewireStrategy& self = *static_cast<RewireStrategy*>(this);
        tr1::bernoulli_distribution coin;

        //try randomly drawn pairs of edges until one satisfies all the
        //consistency checks
        pair<size_t, bool> es, et;

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
                    if((source(_edges[e], _g) == target(es, _edges, _g)) ||
                       (target(_edges[e], _g) == source(es, _edges, _g)))
                        continue;
                }

                // reject parallel edges if not allowed
                if (!parallel_edges && (es.first != e))
                {
                    if (swap_edge::parallel_check_source(e, es, _edges, _g))
                        continue;
                }

                if (e != es.first)
                {
                    self.update_edge(e, false, true);
                    self.update_edge(es.first, false, false);

                    swap_edge::swap_source(e, es, _edges, _edge_index, _g);

                    self.update_edge(e, true, true);
                    self.update_edge(es.first, true, false);
                }
            }
            else
            {
                // rewire target
                et = self.get_target_edge(e);

                if (!self_loops) // reject self-loops if not allowed
                {
                    if((source(_edges[e], _g) == target(et, _edges, _g)) ||
                       (target(_edges[e], _g) == source(et, _edges, _g)))
                        continue;
                }

                // reject parallel edges if not allowed
                if (!parallel_edges && (et.first != e))
                {
                    if (swap_edge::parallel_check_target(e, et, _edges, _g))
                        continue;
                }

                if (e != et.first)
                {
                    self.update_edge(e, false, true);
                    self.update_edge(et.first, false, false);

                    swap_edge::swap_target(e, et, _edges, _edge_index, _g);

                    self.update_edge(e, true, true);
                    self.update_edge(et.first, true, false);
                }
            }
            break;
        }
    }

protected:
    Graph& _g;
    EdgeIndexMap _edge_index;
    vector<edge_t>& _edges;
    rng_t& _rng;
};

// this will rewire the edges so that the combined (in, out) degree distribution
// will be the same, but all the rest is random
template <class Graph, class EdgeIndexMap, class CorrProb>
class RandomRewireStrategy:
    public RewireStrategyBase<Graph, EdgeIndexMap,
                              RandomRewireStrategy<Graph, EdgeIndexMap,
                                                   CorrProb> >
{
public:
    typedef RewireStrategyBase<Graph, EdgeIndexMap,
                               RandomRewireStrategy<Graph, EdgeIndexMap,
                                                    CorrProb> >
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
                         vector<edge_t>& edges, CorrProb corr, rng_t& rng)
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
            update_edge(edge_index[*e_i], true);
        }
    }

    pair<size_t,bool> get_edge(size_t e, bool src)
    {
        tr1::unordered_map<size_t, edge_set_t>& hub =
            (src || !is_directed::apply<Graph>::type::value) ?
            _hub_non_adjacent : _hub_non_incident;

        typeof(hub.begin()) iter =
            hub.find(src ? target(base_t::_edges[e], _g) :
                     source(base_t::_edges[e], _g));
        if (iter != hub.end())
        {
            edge_set_t& eset = iter->second;

            tr1::uniform_int<> sample(0, eset.size());
            size_t i = sample(base_t::_rng);
            if (i == eset.size()) // no rewire option
                return make_pair(e, false);
            else
                return *(eset.get<random_index>().begin() + i);
        }

        tr1::uniform_int<> sample(0, _all_edges.size()-1);
        return _all_edges[sample(base_t::_rng)];
    }


    pair<size_t,bool> get_source_edge(size_t e)
    {
        return get_edge(e, true);
    }

    pair<size_t,bool> get_target_edge(size_t e)
    {
        return get_edge(e, false);
    }

    void update_edge_list(index_t v, const pair<size_t, bool>& e,
                          edge_set_t& edge_set, bool src, bool insert)
    {
        pair<index_t,bool> ep;
        if ((src && source(e, base_t::_edges, _g) != v) ||
            (!src && target(e, base_t::_edges, _g) != v))
        {
            if (insert)
                edge_set.get<hash_index>().insert(e);
            else
                edge_set.get<hash_index>().erase(e);
        }
    }

    void update_edge(size_t e, bool insert, bool final = false)
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
template <class Graph, class EdgeIndexMap, class CorrProb>
class CorrelatedRewireStrategy:
    public RewireStrategyBase<Graph, EdgeIndexMap,
                              CorrelatedRewireStrategy<Graph, EdgeIndexMap,
                                                       CorrProb> >
{
public:
    typedef RewireStrategyBase<Graph, EdgeIndexMap,
                               CorrelatedRewireStrategy<Graph, EdgeIndexMap,
                                                        CorrProb> >
        base_t;

    typedef Graph graph_t;
    typedef EdgeIndexMap edge_index_t;

    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    typedef typename graph_traits<Graph>::edge_descriptor edge_t;
    typedef typename EdgeIndexMap::value_type index_t;

    CorrelatedRewireStrategy (Graph& g, EdgeIndexMap edge_index,
                              vector<edge_t>& edges, CorrProb corr, rng_t& rng)
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

    pair<size_t,bool> get_edge(size_t e, bool src)
    {
        pair<size_t, size_t> deg =
            make_pair(in_degreeS()(src ? source(base_t::_edges[e], _g)
                                   : target(base_t::_edges[e], _g), _g),
                      out_degree(src ? source(base_t::_edges[e], _g) :
                                 target(base_t::_edges[e], _g), _g));
        edges_by_end_deg_t& edges = src ? _edges_by_source : _edges_by_target;
        typename edges_by_end_deg_t::mapped_type& elist = edges[deg];
        tr1::uniform_int<> sample(0, elist.size()-1);
        return elist[sample(base_t::_rng)];
    }

    pair<size_t,bool> get_source_edge(size_t e)
    {
        return get_edge(e, true);
    }

    pair<size_t,bool> get_target_edge(size_t e)
    {
        return get_edge(e, false);
    }

    void update_edge(size_t e, bool insert, bool) {}

private:
    typedef tr1::unordered_map<pair<size_t, size_t>,
                               vector<pair<index_t, bool> >,
                               hash<pair<size_t, size_t> > > edges_by_end_deg_t;
    edges_by_end_deg_t _edges_by_source, _edges_by_target;

protected:
    const Graph& _g;
};

template <class Graph, class EdgeIndexMap, class CorrProb>
class ProbabilisticRewireStrategy:
    public RewireStrategyBase<Graph, EdgeIndexMap,
                              ProbabilisticRewireStrategy<Graph, EdgeIndexMap,
                                                          CorrProb> >
{
public:
    typedef RewireStrategyBase<Graph, EdgeIndexMap,
                                ProbabilisticRewireStrategy<Graph, EdgeIndexMap,
                                                            CorrProb> >
        base_t;

    typedef Graph graph_t;
    typedef EdgeIndexMap edge_index_t;

    typedef pair<size_t,size_t> deg_t;

    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    typedef typename graph_traits<Graph>::edge_descriptor edge_t;
    typedef typename EdgeIndexMap::value_type index_t;

    ProbabilisticRewireStrategy(Graph& g, EdgeIndexMap edge_index,
                                vector<edge_t>& edges,
                                CorrProb corr_prob, rng_t& rng)
        : base_t(g, edge_index, edges, rng), _g(g)
    {
        typename graph_traits<Graph>::vertex_iterator v_i, v_i_end;
        for (tie(v_i, v_i_end) = boost::vertices(g); v_i != v_i_end; ++v_i)
            _deg_count[get_deg(*v_i, g)]++;


        typename graph_traits<Graph>::edge_iterator e_i, e_i_end;
        set<pair<size_t, size_t> > deg_set;
        for (tie(e_i, e_i_end) = boost::edges(g); e_i != e_i_end; ++e_i)
        {
            update_edge(_edge_index[*e_i], true);
            deg_set.insert(get_deg(source(*e_i, g), g));
            deg_set.insert(get_deg(target(*e_i, g), g));
        }

        for (typeof(deg_set.begin()) s_iter = deg_set.begin();
             s_iter != deg_set.end(); ++s_iter)
        {
             _sample_source_deg[*s_iter] = Sampler<deg_t>(true);
             _sample_target_deg[*s_iter] = Sampler<deg_t>(true);
        }

        for (typeof(deg_set.begin()) s_iter = deg_set.begin();
             s_iter != deg_set.end(); ++s_iter)
            for (typeof(deg_set.begin()) t_iter = deg_set.begin();
                 t_iter != deg_set.end(); ++t_iter)
        {
            double p = corr_prob(*s_iter, *t_iter);
            if (!is_directed::apply<Graph>::type::value ||
                (s_iter->second > 0 && t_iter->first > 0))
            {
                _sample_source_deg[*t_iter].Insert(*s_iter, p);
                _sample_target_deg[*s_iter].Insert(*t_iter, p);
            }
        }
    }

    deg_t get_deg(vertex_t v, Graph& g)
    {
        return make_pair(in_degreeS()(v, g), out_degree(v,g));
    }

    pair<size_t,bool> get_source_edge(size_t e)
    {
        deg_t t_deg = get_deg(target(base_t::_edges[e],_g),_g), deg;
        while (true)
        {
            deg = _sample_source_deg[t_deg](base_t::_rng);
            if (_sample_edge_by_source_deg[deg].Empty())
                _sample_source_deg[t_deg].Remove(deg);
            else
                break;
        }
        pair<size_t, bool> ep = _sample_edge_by_source_deg[deg](base_t::_rng);
        //assert(get_deg(source(ep, base_t::_edges, _g), _g) == deg);
        return ep;
    }

    pair<size_t, bool> get_target_edge(size_t e)
    {
        deg_t s_deg = get_deg(source(base_t::_edges[e],_g),_g), deg;
        while (true)
        {
            deg = _sample_target_deg[s_deg](base_t::_rng);
            if (_sample_edge_by_target_deg[deg].Empty())
                _sample_target_deg[s_deg].Remove(deg);
            else
                break;
        }
        pair<size_t, bool> ep = _sample_edge_by_target_deg[deg](base_t::_rng);
        //assert(get_deg(target(ep, base_t::_edges, _g), _g) == deg);
        return ep;
    }

    void update_edge(size_t e, bool insert, bool final = false)
    {
        deg_t s_deg = get_deg(source(base_t::_edges[e], _g), _g),
            t_deg = get_deg(target(base_t::_edges[e], _g), _g);

        pair<size_t, bool> ep = make_pair(e, false);

        if (insert && !final)
        {
            _sample_edge_by_source_deg[s_deg].Insert(ep);
            _sample_edge_by_target_deg[t_deg].Insert(ep);
            if (!is_directed::apply<Graph>::type::value)
            {
                ep.second = true;
                _sample_edge_by_source_deg[t_deg].Insert(ep);
                _sample_edge_by_target_deg[s_deg].Insert(ep);
            }
        }

        if (!insert)
        {
            _sample_edge_by_source_deg[s_deg].Remove(ep);
            _sample_edge_by_target_deg[t_deg].Remove(ep);
            if (!is_directed::apply<Graph>::type::value)
            {
                ep.second = true;
                _sample_edge_by_source_deg[t_deg].Remove(ep);
                _sample_edge_by_target_deg[s_deg].Remove(ep);
            }
        }
    }

private:
    Graph& _g;
    EdgeIndexMap _edge_index;

    typedef tr1::unordered_map<pair<deg_t, deg_t>, double,
                               hash<pair<deg_t, deg_t> > > corr_prob_t;
    typedef tr1::unordered_map<deg_t, Sampler<deg_t>, hash<deg_t> >
        deg_sampler_t;

    typedef tr1::unordered_map<deg_t, Sampler<pair<size_t,bool> >,
                               hash<deg_t> >
        edge_sampler_t;

    edge_sampler_t _sample_edge_by_source_deg;
    edge_sampler_t _sample_edge_by_target_deg;

    deg_sampler_t _sample_source_deg;
    deg_sampler_t _sample_target_deg;

    tr1::unordered_map<deg_t, size_t, hash<deg_t> > _deg_count;
};

} // graph_tool namespace

#endif // GRAPH_REWIRING_HH
