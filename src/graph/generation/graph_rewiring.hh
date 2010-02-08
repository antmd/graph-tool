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

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_util.hh"

namespace graph_tool
{
using namespace std;
using namespace boost;


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

// this functor will swap the source of the edge e with the source of edge se
// and the target of edge e with the target of te
struct swap_edge_triad
{
    template <class Graph, class NewEdgeMap>
    static bool
    parallel_check(const typename graph_traits<Graph>::edge_descriptor& e,
                   const pair<typename graph_traits<Graph>::edge_descriptor,
                              bool>& se,
                   const pair<typename graph_traits<Graph>::edge_descriptor,
                              bool>& te,
                   NewEdgeMap edge_is_new, const Graph &g)
    {
        // We want to check that if we swap the source of 'e' with the source of
        // 'se', and the target of 'te' with the target of 'e', as such
        //
        //  (s)    -e--> (t)          (ns)   -e--> (nt)
        //  (ns)   -se-> (se_t)   =>  (s)    -se-> (se_t)
        //  (te_s) -te-> (nt)         (te_s) -te-> (t),
        //
        // no parallel edges are introduced. We must considered only "new
        // edges", i.e., edges which were already sampled and swapped. "Old
        // edges" will have their chance of being swapped, and then they'll be
        // checked for parallelism.

        typename graph_traits<Graph>::vertex_descriptor
            s = source(e, g),          // current source
            t = target(e, g),          // current target
            ns = source(se, g),        // new source
            nt = target(te, g),        // new target
            te_s = source(te, g),      // target edge source
            se_t = target(se, g);      // source edge target


        if (edge_is_new[se.first] && (ns == s) && (nt == se_t))
            return true; // e is parallel to se after swap
        if (edge_is_new[te.first] && (te_s == ns) && (nt == t))
            return true; // e is parallel to te after swap
        if (edge_is_new[te.first] && edge_is_new[se.first] && (te != se) &&
             (s == te_s) && (t == se_t))
            return true; // se is parallel to te after swap
        if (is_adjacent_in_new(ns,  nt, edge_is_new, g))
            return true; // e would clash with an existing (new) edge
        if (edge_is_new[te.first] && is_adjacent_in_new(te_s, t, edge_is_new, g))
            return true; // te would clash with an existing (new) edge
        if (edge_is_new[se.first] && is_adjacent_in_new(s, se_t, edge_is_new, g))
            return true; // se would clash with an existing (new) edge
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
                nse = add_edge(source(e, g), target(se, g), g).first;
                edge_index[nse] = edge_index[se.first];
                remove_edge(se.first, g);
                edges[edge_index[nse]] = nse;
            }
            if(e != te.first)
            {
                nte = add_edge(source(te, g), target(e, g), g).first;
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

        RewireStrategy<Graph, EdgeIndexMap> rewire(g, edge_index, rng);

        vector<edge_t> edges(num_edges(g));
        vector<bool> is_edge(num_edges(g), false);
        typename graph_traits<Graph>::edge_iterator e, e_end;
        for (tie(e, e_end) = boost::edges(g); e != e_end; ++e)
        {
            if (edge_index[*e] >= edges.size())
            {
                edges.resize(edge_index[*e] + 1);
                is_edge.resize(edge_index[*e] + 1, false);
            }
            edges[edge_index[*e]] = *e;
            is_edge[edge_index[*e]] = true;
        }

        // for each edge simultaneously rewire its source and target
        for (size_t i = 0; i < edges.size(); ++i)
        {
            if (!is_edge[i])
                continue;
            typename graph_traits<Graph>::edge_descriptor e = edges[i];
            typename graph_traits<Graph>::edge_descriptor se, te;
            rewire(e, edges, is_edge, self_loops, parallel_edges);
        }
    }
};

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

// this will rewire the edges so that the resulting graph will be entirely
// random (i.e. Erdos-Renyi)
template <class Graph, class EdgeIndexMap>
class ErdosRewireStrategy
{
public:
    typedef typename graph_traits<Graph>::edge_descriptor edge_t;
    typedef typename EdgeIndexMap::value_type index_t;

    ErdosRewireStrategy(Graph& g, EdgeIndexMap edge_index, rng_t& rng)
        : _g(g), _edge_index(edge_index), _vertices(HardNumVertices()(g)),
          _rng(rng)
    {
        typeof(_vertices.begin()) viter = _vertices.begin();
        typename graph_traits<Graph>::vertex_iterator v, v_end;
        for (tie(v, v_end) = vertices(_g); v != v_end; ++v)
            *(viter++) = *v;
    }

    template<class EdgesType>
    void operator()(const edge_t& e, EdgesType& edges, vector<bool>& is_edge,
                    bool self_loops, bool parallel_edges)
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
                swap_edge_triad::is_adjacent_in_new(s, t, _edge_is_new, _g))
                continue;  // reject parallel edges if not allowed
            break;
        }
        edge_t ne = add_edge(s, t, _g).first;
        edges[_edge_index[e]] = ne;
        remove_edge(e, _g);
        if (_edge_index[ne] >= edges.size())
        {
            edges.resize(_edge_index[ne] + 1);
            is_edge.resize(_edge_index[ne] + 1, false);
        }
        edges[_edge_index[ne]] = ne;
        is_edge[_edge_index[ne]] = true;

        _edge_is_new[ne] = true;
    }

private:
    Graph& _g;
    EdgeIndexMap _edge_index;
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

    RewireStrategyBase(Graph& g, EdgeIndexMap edge_index, rng_t& rng)
        : _g(g), _edge_index(edge_index), _edge_is_new(edge_index), _rng(rng) {}

    template<class EdgesType>
    void operator()(const edge_t& e, EdgesType& edges, vector<bool>& is_edge,
                    bool self_loops, bool parallel_edges)
    {
        // where should we sample the edges from
        vector<pair<index_t,bool> >* edges_source=0, *edges_target=0;
        static_cast<RewireStrategy*>(this)->get_edges(e, edges_source,
                                                      edges_target);


        //try randomly drawn pairs of edges until one satisfies all the
        //consistency checks
        bool found = false;
        pair<edge_t,bool> es, et;
        typedef random_permutation_iterator
            <typename vector<pair<index_t,bool> >::iterator, rng_t>
            random_edge_iter;

        random_edge_iter esi(edges_source->begin(), edges_source->end(),
                             _rng),
                         esi_end(edges_source->end(), edges_source->end(),
                             _rng);
        for (; esi != esi_end && !found; ++esi)
        {
            if (!is_edge[(*esi).first])
                continue;
            es = make_pair(edges[(*esi).first], (*esi).second);
            static_cast<RewireStrategy*>(this)->check_source_edge(es, e);

            if(!self_loops) // reject self-loops if not allowed
            {
                if((source(e, _g) == target(es, _g)))
                    continue;
            }

            random_edge_iter eti(edges_target->begin(), edges_target->end(),
                                 _rng),
                             eti_end(edges_target->end(), edges_target->end(),
                                 _rng);
            for (; eti != eti_end && !found; ++eti)
            {
                if (!is_edge[(*eti).first])
                    continue;
                et = make_pair(edges[(*eti).first], (*eti).second);
                static_cast<RewireStrategy*>(this)->check_target_edge(et, e);

                if (!self_loops) // reject self-loops if not allowed
                {
                    if ((source(es, _g) == target(et, _g)) ||
                        (source(et, _g) == target(e, _g)))
                        continue;
                }
                if (!parallel_edges) // reject parallel edges if not allowed
                {
                    if (swap_edge_triad::parallel_check(e, es, et, _edge_is_new,
                                                        _g))
                        continue;
                }
                found = true;
            }
        }
        if (!found)
            throw GraphException("Couldn't find random pair of edges to swap"
                                 "... This is a bug.");
        _edge_is_new[e] = true;
        swap_edge_triad()(e, es, et, edges, _edge_index, _g);
    }

private:
    Graph& _g;
    EdgeIndexMap _edge_index;
    checked_vector_property_map<bool, EdgeIndexMap> _edge_is_new;
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

    RandomRewireStrategy(Graph& g, EdgeIndexMap edge_index,
                         rng_t& rng)
        : base_t(g, edge_index, rng)
    {
        //bernoulli_distribution coin(0.5);
        typename graph_traits<Graph>::edge_iterator e_i, e_i_end;
        for (tie(e_i, e_i_end) = edges(g); e_i != e_i_end; ++e_i)
        {
            // For undirected graphs, there is no difference between source
            // and target, and each edge will appear _twice_ on the list
            // once for each different ordering of source and target.

            _all_edges.push_back(make_pair(edge_index[*e_i], false));
            if (!is_directed::apply<Graph>::type::value)
                _all_edges.push_back(make_pair(edge_index[*e_i], true));
        }
        _all_edges2 = _all_edges;

    }

    void get_edges(const edge_t& e,
                   vector<pair<index_t,bool> >*& edges_source,
                   vector<pair<index_t,bool> >*& edges_target)
    {
        edges_source = &_all_edges;
        edges_target = &_all_edges2;
    }

    void check_source_edge(pair<edge_t,bool>& se, const edge_t& e) {}
    void check_target_edge(pair<edge_t,bool>& te, const edge_t& e) {}

private:
    vector<pair<index_t,bool> > _all_edges;
    vector<pair<index_t,bool> > _all_edges2;
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
                              rng_t& rng) : base_t(g, edge_index, rng), _g(g)
    {
        typename graph_traits<Graph>::edge_iterator e_i, e_i_end;
        for (tie(e_i, e_i_end) = edges(g); e_i != e_i_end; ++e_i)
        {
            // For undirected graphs, there is no difference between source and
            // target, and each edge will appear _twice_ on the lists below,
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

    void get_edges(const edge_t& e, vector<pair<index_t,bool> >*& edges_source,
                   vector<pair<index_t,bool> >*& edges_target)
    {
        pair<size_t, size_t> deg_source =
            make_pair(in_degreeS()(source(e, _g), _g),
                      out_degree(source(e, _g), _g));
        edges_source = &_edges_by_source[deg_source];


        pair<size_t, size_t> deg_target =
            make_pair(in_degreeS()(target(e, _g), _g),
                      out_degree(target(e, _g), _g));

        edges_target = &_edges_by_target[deg_target];
    }


    void check_source_edge(pair<edge_t,bool>& se, const edge_t& e)
    {
        check_source_edge_if_undirected
            (se, e, typename is_directed::apply<Graph>::type());
    }
    void check_target_edge(pair<edge_t,bool>& te, const edge_t& e)
    {
        check_target_edge_if_undirected
            (te, e, typename is_directed::apply<Graph>::type());
    }

    void check_source_edge_if_undirected(pair<edge_t,bool>& se, const edge_t& e,
                                         boost::true_type) {}
    void check_target_edge_if_undirected(pair<edge_t,bool>& te, const edge_t& e,
                                         boost::true_type) {}

    void check_source_edge_if_undirected(pair<edge_t,bool>& se, const edge_t& e,
                                         boost::false_type)
    {
        // check if the edge direction is correct, otherwise invert it.
        pair<size_t, size_t> deg_source1 =
            make_pair(out_degree(source(e, _g), _g),
                      out_degree(source(e, _g), _g));

        pair<size_t, size_t> deg_source2 =
            make_pair(out_degree(source(se, _g), _g),
                      out_degree(source(se, _g), _g));

        if (deg_source1 != deg_source2)
            se.first = edge_t(se.first, !se.first.IsInverted());
    }

    void check_target_edge_if_undirected(pair<edge_t,bool>& te, const edge_t& e,
                                         boost::false_type)
    {
        // check if the edge direction is correct, otherwise invert it.
        pair<size_t, size_t> deg_target1 =
            make_pair(out_degree(target(e, _g), _g),
                      out_degree(target(e, _g), _g));

        pair<size_t, size_t> deg_target2 =
            make_pair(out_degree(target(te, _g), _g),
                      out_degree(target(te, _g), _g));

        if (deg_target1 != deg_target2)
            te.first = edge_t(te.first, !te.first.IsInverted());
    }

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
