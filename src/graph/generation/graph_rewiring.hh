// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2013 Tiago de Paula Peixoto <tiago@skewed.de>
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

#if (GCC_VERSION >= 40400)
#   include <tr1/unordered_set>
#else
#   include <boost/tr1/unordered_set.hpp>
#endif
#include <boost/functional/hash.hpp>

#include <iostream>

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_util.hh"
#include "sampler.hh"

#include "random.hh"

namespace graph_tool
{
using namespace std;
using namespace boost;

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


// this functor will swap the source of the edge e with the source of edge se
// and the target of edge e with the target of te
struct swap_edge
{
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
            return true; // e would clash with an existing edge
        if (is_adjacent(te_s, t, g))
            return true; // te would clash with an existing edge
        return false; // the coast is clear - hooray!
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

        edge_index[nte] = edge_index[edges[te.first]];
        remove_edge(edges[te.first], g);
        edges[te.first] = nte;

        edge_index[ne] = edge_index[edges[e]];
        remove_edge(edges[e], g);
        edges[e] = ne;
    }

};

// used for verbose display
void print_progress(size_t i, size_t n_iter, size_t current, size_t total,
                    stringstream& str)
{
    size_t atom = (total > 200) ? total / 100 : 1;
    if ( ( (current+1) % atom == 0) || (current + 1) == total)
    {
        size_t size = str.str().length();
        for (size_t j = 0; j < str.str().length(); ++j)
            cout << "\b";
        str.str("");
        str << "(" << i + 1 << " / " << n_iter << ") "
            << current + 1 << " of " << total << " ("
            << (current + 1) * 100 / total << "%)";
        for (int j = 0; j < int(size - str.str().length()); ++j)
            str << " ";
        cout << str.str() << flush;
    }
}

//select blocks based on in/out degrees
class DegreeBlock
{
public:
    typedef pair<size_t, size_t> block_t;

    template <class Graph>
    block_t get_block(typename graph_traits<Graph>::vertex_descriptor v,
                      const Graph& g) const
    {
        return make_pair(in_degreeS()(v, g), out_degree(v, g));
    }
};

//select blocks based on property map
template <class PropertyMap>
class PropertyBlock
{
public:
    typedef typename property_traits<PropertyMap>::value_type block_t;

    PropertyBlock(PropertyMap p): _p(p) {}

    template <class Graph>
    block_t get_block(typename graph_traits<Graph>::vertex_descriptor v,
                      const Graph& g) const
    {
        return get(_p, v);
    }

private:
    PropertyMap _p;
};

// main rewire loop
template <template <class Graph, class EdgeIndexMap, class CorrProb,
                    class BlockDeg>
          class RewireStrategy>
struct graph_rewire
{

    template <class Graph, class EdgeIndexMap, class CorrProb,
              class BlockDeg>
    void operator()(Graph& g, EdgeIndexMap edge_index, CorrProb corr_prob,
                    bool self_loops, bool parallel_edges,
                    pair<size_t, bool> iter_sweep,
                    pair<bool, bool> cache_verbose,
                    size_t& pcount, rng_t& rng, BlockDeg bd)
        const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;
        bool cache = cache_verbose.first;
        bool verbose = cache_verbose.second;

        vector<edge_t> edges;
        vector<size_t> edge_pos;
        typedef random_permutation_iterator<typename vector<size_t>::iterator,
                                            rng_t>
            random_edge_iter;

        typename graph_traits<Graph>::edge_iterator e, e_end;
        for (tie(e, e_end) = boost::edges(g); e != e_end; ++e)
        {
            edges.push_back(*e);
            edge_pos.push_back(edge_pos.size());
        }

        RewireStrategy<Graph, EdgeIndexMap, CorrProb, BlockDeg>
            rewire(g, edge_index, edges, corr_prob, bd, cache, rng);

        size_t niter;
        bool no_sweep;
        tie(niter, no_sweep) = iter_sweep;
        pcount = 0;
        if (verbose)
            cout << "rewiring edges: ";
        stringstream str;
        for (size_t i = 0; i < niter; ++i)
        {
            random_edge_iter
                ei_begin(edge_pos.begin(), edge_pos.end(), rng),
                ei_end(edge_pos.end(), edge_pos.end(), rng);

            // for each edge rewire its source or target
            for (random_edge_iter ei = ei_begin; ei != ei_end; ++ei)
            {
                size_t e_pos = ei - ei_begin;
                if (verbose)
                    print_progress(i, niter, e_pos, no_sweep ? 1 : edges.size(),
                                   str);
                bool success = rewire(*ei, self_loops, parallel_edges);
                if (!success)
                    ++pcount;

                if (no_sweep)
                    break;
            }
        }
        if (verbose)
            cout << endl;
    }

    template <class Graph, class EdgeIndexMap, class CorrProb>
    void operator()(Graph& g, EdgeIndexMap edge_index, CorrProb corr_prob,
                    bool self_loops, bool parallel_edges,
                    pair<size_t, bool> iter_sweep,
                    pair<bool, bool> cache_verbose,
                    size_t& pcount, rng_t& rng)
        const
    {
        operator()(g, edge_index, corr_prob, self_loops, parallel_edges,
                   iter_sweep, cache_verbose, pcount, rng, DegreeBlock());
    }
};


// this will rewire the edges so that the resulting graph will be entirely
// random (i.e. Erdos-Renyi)
template <class Graph, class EdgeIndexMap, class CorrProb, class BlockDeg>
class ErdosRewireStrategy
{
public:
    typedef typename graph_traits<Graph>::edge_descriptor edge_t;
    typedef typename EdgeIndexMap::value_type index_t;

    ErdosRewireStrategy(Graph& g, EdgeIndexMap edge_index,
                        vector<edge_t>& edges, CorrProb, BlockDeg,
                        bool, rng_t& rng)
        : _g(g), _edge_index(edge_index), _edges(edges),
          _vertices(HardNumVertices()(g)), _rng(rng)
    {
        typeof(_vertices.begin()) viter = _vertices.begin();
        typename graph_traits<Graph>::vertex_iterator v, v_end;
        for (tie(v, v_end) = vertices(_g); v != v_end; ++v)
            *(viter++) = *v;
    }

    bool operator()(size_t ei, bool self_loops, bool parallel_edges)
    {
        //try randomly drawn pairs of vertices
        typedef random_permutation_iterator
            <typename graph_traits<Graph>::vertex_iterator, rng_t>
            random_vertex_iter;

        tr1::uniform_int<size_t> sample(0, _vertices.size() - 1);
        typename graph_traits<Graph>::vertex_descriptor s, t;
        while (true)
        {
            s = sample(_rng);
            t = sample(_rng);

            // reject self-loops if not allowed
            if(s == t && !self_loops)
                continue;
            break;
        }

        // reject parallel edges if not allowed
        if (!parallel_edges && is_adjacent(s, t, _g))
            return false;

        edge_t ne = add_edge(s, t, _g).first;
        _edge_index[ne] = _edge_index[_edges[ei]];
        remove_edge(_edges[ei], _g);
        _edges[ei] = ne;

        return true;
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

    bool operator()(size_t ei, bool self_loops, bool parallel_edges)
    {
        RewireStrategy& self = *static_cast<RewireStrategy*>(this);

        // try randomly drawn pairs of edges and check if they satisfy all the
        // consistency checks

        // rewire target
        pair<size_t, bool> et = self.get_target_edge(ei);

        if (!self_loops) // reject self-loops if not allowed
        {
            if((source(_edges[ei], _g) == target(et, _edges, _g)) ||
               (target(_edges[ei], _g) == source(et, _edges, _g)))
                return false;
        }

        // reject parallel edges if not allowed
        if (!parallel_edges && (et.first != ei))
        {
            if (swap_edge::parallel_check_target(ei, et, _edges, _g))
                return false;
        }

        if (ei != et.first)
        {
            self.update_edge(ei, false);
            self.update_edge(et.first, false);

            swap_edge::swap_target(ei, et, _edges, _edge_index, _g);

            self.update_edge(ei, true);
            self.update_edge(et.first, true);
        }
        else
        {
            return false;
        }

        return true;
    }

protected:
    Graph& _g;
    EdgeIndexMap _edge_index;
    vector<edge_t>& _edges;
    rng_t& _rng;
};

// this will rewire the edges so that the combined (in, out) degree distribution
// will be the same, but all the rest is random
template <class Graph, class EdgeIndexMap, class CorrProb, class BlockDeg>
class RandomRewireStrategy:
    public RewireStrategyBase<Graph, EdgeIndexMap,
                              RandomRewireStrategy<Graph, EdgeIndexMap,
                                                   CorrProb, BlockDeg> >
{
public:
    typedef RewireStrategyBase<Graph, EdgeIndexMap,
                               RandomRewireStrategy<Graph, EdgeIndexMap,
                                                    CorrProb, BlockDeg> >
        base_t;

    typedef Graph graph_t;
    typedef EdgeIndexMap edge_index_t;

    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    typedef typename graph_traits<Graph>::edge_descriptor edge_t;
    typedef typename EdgeIndexMap::value_type index_t;

    struct hash_index {};
    struct random_index {};

    RandomRewireStrategy(Graph& g, EdgeIndexMap edge_index,
                         vector<edge_t>& edges, CorrProb, BlockDeg,
                         bool, rng_t& rng)
        : base_t(g, edge_index, edges, rng), _g(g) {}

    pair<size_t,bool> get_target_edge(size_t e)
    {
        tr1::uniform_int<> sample(0, base_t::_edges.size() - 1);
        pair<size_t, bool> et = make_pair(sample(base_t::_rng), false);
        if (!is_directed::apply<Graph>::type::value)
        {
            tr1::bernoulli_distribution coin(0.5);
            et.second = coin(base_t::_rng);
        }
        return et;
    }

    void update_edge(size_t e, bool insert) {}

private:
    Graph& _g;
    EdgeIndexMap _edge_index;
};


// this will rewire the edges so that the (in,out) degree distributions and the
// (in,out)->(in,out) correlations will be the same, but all the rest is random
template <class Graph, class EdgeIndexMap, class CorrProb, class BlockDeg>
class CorrelatedRewireStrategy:
    public RewireStrategyBase<Graph, EdgeIndexMap,
                              CorrelatedRewireStrategy<Graph, EdgeIndexMap,
                                                       CorrProb, BlockDeg> >
{
public:
    typedef RewireStrategyBase<Graph, EdgeIndexMap,
                               CorrelatedRewireStrategy<Graph, EdgeIndexMap,
                                                        CorrProb, BlockDeg> >
        base_t;

    typedef Graph graph_t;

    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    typedef typename graph_traits<Graph>::edge_descriptor edge_t;

    CorrelatedRewireStrategy(Graph& g, EdgeIndexMap edge_index,
                             vector<edge_t>& edges, CorrProb, BlockDeg,
                             bool, rng_t& rng)
        : base_t(g, edge_index, edges, rng), _g(g)
    {
        for (size_t ei = 0; ei < base_t::_edges.size(); ++ei)
        {
            // For undirected graphs, there is no difference between source and
            // target, and each edge will appear _twice_ in the list below,
            // once for each different ordering of source and target.
            edge_t& e = base_t::_edges[ei];

            vertex_t t = target(e, _g);
            deg_t tdeg = make_pair(in_degreeS()(t, _g), out_degree(t, _g));
            _edges_by_target[tdeg].push_back(make_pair(ei, false));

            if (!is_directed::apply<Graph>::type::value)
            {
                t = source(e, _g);
                tdeg = make_pair(in_degreeS()(t, _g), out_degree(t, _g));
                _edges_by_target[tdeg].push_back(make_pair(ei, true));
            }
        }
    }

    pair<size_t,bool> get_target_edge(size_t ei)
    {
        edge_t& e = base_t::_edges[ei];
        vertex_t t = target(e, _g);
        deg_t tdeg = make_pair(in_degreeS()(t, _g), out_degree(t, _g));
        typename edges_by_end_deg_t::mapped_type& elist =
            _edges_by_target[tdeg];
        tr1::uniform_int<> sample(0, elist.size() - 1);

        return elist[sample(base_t::_rng)];
    }

    void update_edge(size_t e, bool insert) {}

private:
    typedef pair<size_t, size_t> deg_t;
    typedef tr1::unordered_map<deg_t,
                               vector<pair<size_t, bool> >,
                               hash<deg_t> >
        edges_by_end_deg_t;
    edges_by_end_deg_t _edges_by_target;

protected:
    const Graph& _g;
};


// general stochastic blockmodel
template <class Graph, class EdgeIndexMap, class CorrProb, class BlockDeg>
class ProbabilisticRewireStrategy:
    public RewireStrategyBase<Graph, EdgeIndexMap,
                              ProbabilisticRewireStrategy<Graph, EdgeIndexMap,
                                                          CorrProb, BlockDeg> >
{
public:
    typedef RewireStrategyBase<Graph, EdgeIndexMap,
                               ProbabilisticRewireStrategy<Graph, EdgeIndexMap,
                                                           CorrProb, BlockDeg> >
        base_t;

    typedef Graph graph_t;
    typedef EdgeIndexMap edge_index_t;

    typedef typename BlockDeg::block_t deg_t;

    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    typedef typename graph_traits<Graph>::edge_descriptor edge_t;
    typedef typename EdgeIndexMap::value_type index_t;

    ProbabilisticRewireStrategy(Graph& g, EdgeIndexMap edge_index,
                                vector<edge_t>& edges, CorrProb corr_prob,
                                BlockDeg blockdeg, bool cache, rng_t& rng)
        : base_t(g, edge_index, edges, rng), _g(g), _corr_prob(corr_prob),
          _blockdeg(blockdeg)
    {
        tr1::unordered_set<deg_t, hash<deg_t> > deg_set;
        for (size_t ei = 0; ei < base_t::_edges.size(); ++ei)
        {
            edge_t& e = base_t::_edges[ei];
            deg_set.insert(get_deg(source(e, g), g));
            deg_set.insert(get_deg(target(e, g), g));
        }

        if (cache)
        {
            // cache probabilities
            for (typeof(deg_set.begin()) s_iter = deg_set.begin();
                 s_iter != deg_set.end(); ++s_iter)
                for (typeof(deg_set.begin()) t_iter = deg_set.begin();
                     t_iter != deg_set.end(); ++t_iter)
                {
                    double p = _corr_prob(*s_iter, *t_iter);
                    if (isnan(p) || isinf(p) || p < 0)
                        p = 0;
                    _probs[make_pair(*s_iter, *t_iter)] = p;
                }
        }
    }

    double get_prob(const deg_t& s_deg, const deg_t& t_deg)
    {
        if (_probs.empty())
        {
            double p = _corr_prob(s_deg, t_deg);
            if (isnan(p) || isinf(p) || p < 0)
                p = 0;
            return p;
        }
        return _probs[make_pair(s_deg, t_deg)];
    }

    deg_t get_deg(vertex_t v, Graph& g)
    {
        return _blockdeg.get_block(v, g);
    }

    pair<size_t, bool> get_target_edge(size_t ei)
    {
        deg_t s_deg = get_deg(source(base_t::_edges[ei], _g), _g);
        deg_t t_deg = get_deg(target(base_t::_edges[ei], _g), _g);

        tr1::uniform_int<> sample(0, base_t::_edges.size() - 1);
        size_t epi = sample(base_t::_rng);
        pair<size_t, bool> ep = make_pair(epi, false);
        if (!is_directed::apply<Graph>::type::value)
        {
            // for undirected graphs we must select a random direction
            tr1::bernoulli_distribution coin(0.5);
            ep.second = coin(base_t::_rng);
        }

        deg_t ep_s_deg = get_deg(source(ep, base_t::_edges, _g), _g);
        deg_t ep_t_deg = get_deg(target(ep, base_t::_edges, _g), _g);

        double pi = (get_prob(s_deg, t_deg) *
                     get_prob(ep_s_deg, ep_t_deg));
        double pf = (get_prob(s_deg, ep_t_deg) *
                     get_prob(ep_s_deg, t_deg));

        if (pf >= pi || pi == 0)
            return ep;

        if (pf == 0)
            return make_pair(ei, false); // reject

        double a = exp(log(pf) - log(pi));

        tr1::variate_generator<rng_t&, tr1::uniform_real<> >
            rsample(base_t::_rng, tr1::uniform_real<>(0.0, 1.0));
        double r = rsample();
        if (r > a)
            return make_pair(ei, false); // reject
        else
            return ep;
    }

    void update_edge(size_t e, bool insert) {}

private:
    Graph& _g;
    EdgeIndexMap _edge_index;
    CorrProb _corr_prob;
    BlockDeg _blockdeg;
    tr1::unordered_map<pair<deg_t, deg_t>, double, hash<pair<deg_t, deg_t> > >
        _probs;
};

} // graph_tool namespace

#endif // GRAPH_REWIRING_HH
