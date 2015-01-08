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

#ifndef GRAPH_REWIRING_HH
#define GRAPH_REWIRING_HH

#include <unordered_set>
#include <tuple>

#include <boost/functional/hash.hpp>

#include <iostream>

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_util.hh"
#include "sampler.hh"

#include "random.hh"

#ifdef HAVE_SPARSEHASH
#include SPARSEHASH_INCLUDE(dense_hash_map)
#endif

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


template <class Nmap, class Graph>
void add_count(size_t s, size_t t, Nmap& nvmap, Graph& g)
{
    if (!is_directed::apply<Graph>::type::value && s > t)
        std::swap(s, t);
    auto& nmap = nvmap[s];
    nmap[t]++;
}

template <class Nmap, class Graph>
void remove_count(size_t s, size_t t, Nmap& nvmap, Graph& g)
{
    if (!is_directed::apply<Graph>::type::value && s > t)
        std::swap(s, t);
    auto& nmap = nvmap[s];
    auto iter = nmap.find(t);
    iter->second--;
    if (iter->second == 0)
        nmap.erase(iter);
}

template <class Nmap, class Graph>
size_t get_count(size_t s, size_t t, Nmap& nvmap, Graph& g)
{
    if (!is_directed::apply<Graph>::type::value && s > t)
        std::swap(s, t);
    auto& nmap = nvmap[s];
    auto iter = nmap.find(t);
    if (iter == nmap.end())
        return 0;
    return iter->second;
    // if (s != t)
    //     return iter->second;
    // else
    //     return iter->second / 2;
}

// this functor will swap the source of the edge e with the source of edge se
// and the target of edge e with the target of te
struct swap_edge
{
    template <class Nmap, class Graph>
    static bool
    parallel_check_target (const pair<size_t, bool>& e,
                           const pair<size_t, bool>& te,
                           vector<typename graph_traits<Graph>::edge_descriptor>& edges,
                           Nmap& nmap,
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
            s = source(e, edges, g),          // current source
            t = target(e, edges, g),          // current target
            nt = target(te, edges, g),        // new target
            te_s = source(te, edges, g);      // target edge source

        if (get_count(s,  nt, nmap, g) > 0)
            return true; // e would clash with an existing edge
        if (get_count(te_s, t, nmap, g) > 0)
            return true; // te would clash with an existing edge
        return false; // no parallel edges
    }

    template <class Graph>
    static void swap_target
        (const pair<size_t, bool>& e,
         const pair<size_t, bool>& te,
         vector<typename graph_traits<Graph>::edge_descriptor>& edges,
         Graph& g)
    {
        // swap the target of the edge 'e' with the target of edge 'te', as
        // such:
        //
        //  (s)    -e--> (t)          (s)    -e--> (nt)
        //  (te_s) -te-> (nt)   =>    (te_s) -te-> (t)

        if (e.first == te.first)
            return;

        // new edges which will replace the old ones
        typename graph_traits<Graph>::edge_descriptor ne, nte;
        typename graph_traits<Graph>::vertex_descriptor
            s_e = source(e, edges, g),
            t_e = target(e, edges, g),
            s_te = source(te, edges, g),
            t_te = target(te, edges, g);
        remove_edge(edges[e.first], g);
        remove_edge(edges[te.first], g);

        if (is_directed::apply<Graph>::type::value || !e.second)
            ne = add_edge(s_e, t_te, g).first;
        else // keep invertedness (only for undirected graphs)
            ne = add_edge(t_te, s_e, g).first;
        edges[e.first] = ne;
        if (is_directed::apply<Graph>::type::value || !te.second)
            nte = add_edge(s_te, t_e, g).first;
        else // keep invertedness (only for undirected graphs)
            nte = add_edge(t_e, s_te,  g).first;
        edges[te.first] = nte;
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
                      const Graph&) const
    {
        return get(_p, v);
    }

private:
    PropertyMap _p;
};

// select an appropriate "null" key for densehash
template <class Type>
struct get_null_key
{
    Type operator()() const
    {
        return numeric_limits<Type>::max();
    }
};

template <>
struct get_null_key<string>
{
    string operator()() const
    {
        return lexical_cast<string>(get_null_key<size_t>()());
    }
};

template <>
struct get_null_key<boost::python::object>
{
    boost::python::object operator()() const
    {
        return boost::python::object();
    }
};

template <class Type>
struct get_null_key<vector<Type>>
{
    vector<Type> operator()() const
    {
        vector<Type> v(1);
        v[0] = get_null_key<Type>()();
        return v;
    }
};

template <class Type1, class Type2>
struct get_null_key<pair<Type1, Type2>>
{
    pair<Type1, Type2> operator()() const
    {
        return make_pair(get_null_key<Type1>()(),
                         get_null_key<Type2>()());
    }
};


// main rewire loop
template <template <class Graph, class EdgeIndexMap, class CorrProb,
                    class BlockDeg>
          class RewireStrategy>
struct graph_rewire
{

    template <class Graph, class EdgeIndexMap, class CorrProb,
              class BlockDeg, class PinMap>
    void operator()(Graph& g, EdgeIndexMap edge_index, CorrProb corr_prob,
                    PinMap pin, bool self_loops, bool parallel_edges,
                    pair<size_t, bool> iter_sweep,
                    std::tuple<bool, bool, bool> cache_verbose,
                    size_t& pcount, rng_t& rng, BlockDeg bd)
        const
    {
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;
        bool persist = std::get<0>(cache_verbose);
        bool cache = std::get<1>(cache_verbose);
        bool verbose = std::get<2>(cache_verbose);

        vector<edge_t> edges;
        vector<size_t> edge_pos;
        typename graph_traits<Graph>::edge_iterator e, e_end;
        for (tie(e, e_end) = boost::edges(g); e != e_end; ++e)
        {
            if (pin[*e])
                continue;
            edges.push_back(*e);
            edge_pos.push_back(edge_pos.size());
        }

        typedef random_permutation_iterator<typename vector<size_t>::iterator,
                                            rng_t>
            random_edge_iter;

        RewireStrategy<Graph, EdgeIndexMap, CorrProb, BlockDeg>
            rewire(g, edge_index, edges, corr_prob, bd, cache, rng, parallel_edges);

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

            for (random_edge_iter ei = ei_begin; ei != ei_end; ++ei)
            {
                size_t e_pos = ei - ei_begin;
                if (verbose)
                    print_progress(i, niter, e_pos, no_sweep ? 1 : edges.size(),
                                   str);

                size_t e = *ei;

                bool success = false;
                do
                {
                    success = rewire(e, self_loops, parallel_edges);
                }
                while(persist && !success);

                if (!success)
                    ++pcount;

                if (no_sweep)
                    break;
            }
        }
        if (verbose)
            cout << endl;
    }

    template <class Graph, class EdgeIndexMap, class CorrProb, class PinMap>
    void operator()(Graph& g, EdgeIndexMap edge_index, CorrProb corr_prob,
                    PinMap pin, bool self_loops, bool parallel_edges,
                    pair<size_t, bool> iter_sweep,
                    std::tuple<bool, bool, bool> cache_verbose,
                    size_t& pcount, rng_t& rng)
        const
    {
        operator()(g, edge_index, corr_prob, pin, self_loops, parallel_edges,
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
                        bool, rng_t& rng, bool parallel_edges)
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
        std::uniform_int_distribution<size_t> sample(0, _vertices.size() - 1);
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

        remove_edge(_edges[ei], _g);
        edge_t ne = add_edge(s, t, _g).first;
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
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

    typedef typename EdgeIndexMap::value_type index_t;

    RewireStrategyBase(Graph& g, EdgeIndexMap edge_index, vector<edge_t>& edges,
                       rng_t& rng, bool parallel_edges)
        : _g(g), _edge_index(edge_index), _edges(edges), _rng(rng),
          _nmap(get(vertex_index, g), num_vertices(g))
    {
        if (!parallel_edges)
        {
            typename graph_traits<Graph>::vertex_iterator v, v_end;
            for (tie(v, v_end) = boost::vertices(g); v != v_end; ++v)
            {
    #ifdef HAVE_SPARSEHASH
                _nmap[*v].set_empty_key(get_null_key<size_t>()());
                _nmap[*v].set_deleted_key(get_null_key<size_t>()() - 1);
    #endif
            }

            for (size_t i = 0; i < edges.size(); ++i)
                add_count(source(edges[i], g), target(edges[i], g), _nmap, g);
        }
    }

    bool operator()(size_t ei, bool self_loops, bool parallel_edges)
    {
        RewireStrategy& self = *static_cast<RewireStrategy*>(this);

        // try randomly drawn pairs of edges and check if they satisfy all the
        // consistency checks

        pair<size_t, bool> e = make_pair(ei, false);

        // rewire target
        pair<size_t, bool> et = self.get_target_edge(e, parallel_edges);

        if (!self_loops) // reject self-loops if not allowed
        {
            if((source(e, _edges, _g) == target(et, _edges, _g)) ||
               (target(e, _edges, _g) == source(et, _edges, _g)))
                return false;
        }

        // reject parallel edges if not allowed
        if (!parallel_edges && (et.first != e.first))
        {
            if (swap_edge::parallel_check_target(e, et, _edges, _nmap, _g))
                return false;
        }

        if (e.first != et.first)
        {
            self.update_edge(e.first, false);
            self.update_edge(et.first, false);

            if (!parallel_edges)
            {
                remove_count(source(e, _edges, _g), target(e, _edges, _g), _nmap, _g);
                remove_count(source(et, _edges, _g), target(et, _edges, _g), _nmap, _g);
            }

            swap_edge::swap_target(e, et, _edges, _g);

            self.update_edge(e.first, true);
            self.update_edge(et.first, true);

            if (!parallel_edges)
            {
                add_count(source(e, _edges, _g), target(e, _edges, _g), _nmap, _g);
                add_count(source(et, _edges, _g), target(et, _edges, _g), _nmap, _g);
            }
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

#ifdef HAVE_SPARSEHASH
        typedef google::dense_hash_map<size_t, size_t> nmapv_t;
#else
        typedef unordered_map<size_t, size_t> nmapv_t;
#endif
    typedef typename property_map_type::apply<nmapv_t,
                                              typename property_map<Graph, vertex_index_t>::type>
        ::type::unchecked_t nmap_t;
    nmap_t _nmap;
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
                         bool, rng_t& rng, bool parallel_edges)
        : base_t(g, edge_index, edges, rng, parallel_edges), _g(g) {}

    pair<size_t,bool> get_target_edge(pair<size_t,bool>& e, bool parallel_edges)
    {
        std::uniform_int_distribution<> sample(0, base_t::_edges.size() - 1);
        pair<size_t, bool> et = make_pair(sample(base_t::_rng), false);
        if (!is_directed::apply<Graph>::type::value)
        {
            std::bernoulli_distribution coin(0.5);
            et.second = coin(base_t::_rng);
            e.second = coin(base_t::_rng);
        }
        return et;
    }

    void update_edge(size_t, bool) {}

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

    typedef typename BlockDeg::block_t deg_t;

    CorrelatedRewireStrategy(Graph& g, EdgeIndexMap edge_index,
                             vector<edge_t>& edges, CorrProb, BlockDeg blockdeg,
                             bool, rng_t& rng, bool parallel_edges)
        : base_t(g, edge_index, edges, rng, parallel_edges), _blockdeg(blockdeg), _g(g)
    {
#ifdef HAVE_SPARSEHASH
        _edges_by_target.set_empty_key(get_null_key<deg_t>()());
#endif

        for (size_t ei = 0; ei < base_t::_edges.size(); ++ei)
        {
            // For undirected graphs, there is no difference between source and
            // target, and each edge will appear _twice_ in the list below,
            // once for each different ordering of source and target.
            edge_t& e = base_t::_edges[ei];

            vertex_t t = target(e, _g);
            deg_t tdeg = get_deg(t, _g);;
            _edges_by_target[tdeg].push_back(make_pair(ei, false));

            if (!is_directed::apply<Graph>::type::value)
            {
                t = source(e, _g);
                deg_t tdeg = get_deg(t, _g);
                _edges_by_target[tdeg].push_back(make_pair(ei, true));
            }
        }
    }

    pair<size_t,bool> get_target_edge(pair<size_t, bool>& e, bool parallel_edge)
    {
        if (!is_directed::apply<Graph>::type::value)
        {
            std::bernoulli_distribution coin(0.5);
            e.second = coin(base_t::_rng);
        }

        vertex_t t = target(e, base_t::_edges, _g);
        deg_t tdeg = get_deg(t, _g);
        auto& elist = _edges_by_target[tdeg];
        std::uniform_int_distribution<> sample(0, elist.size() - 1);
        auto ep = elist[sample(base_t::_rng)];
        if (get_deg(target(ep, base_t::_edges, _g), _g) != tdeg)
            ep.second = not ep.second;
        return ep;
    }

    void update_edge(size_t, bool) {}

    deg_t get_deg(vertex_t v, const Graph& g)
    {
        return _blockdeg.get_block(v, g);
    }

private:
    BlockDeg _blockdeg;

#ifdef HAVE_SPARSEHASH
    typedef google::dense_hash_map<deg_t,
                                   vector<pair<size_t, bool>>,
                                   boost::hash<deg_t>>
        edges_by_end_deg_t;
#else
    typedef std::unordered_map<deg_t,
                               vector<pair<size_t, bool>>,
                               boost::hash<deg_t>>
        edges_by_end_deg_t;
#endif

    edges_by_end_deg_t _edges_by_target;

protected:
    const Graph& _g;
};


// general stochastic blockmodel
// this version is based on rejection sampling
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
                                BlockDeg blockdeg, bool cache, rng_t& rng,
                                bool parallel_edges)
        : base_t(g, edge_index, edges, rng, parallel_edges), _g(g), _corr_prob(corr_prob),
          _blockdeg(blockdeg)
    {
#ifdef HAVE_SPARSEHASH
        _probs.set_empty_key(get_null_key<pair<deg_t, deg_t>>()());
#endif
        if (cache)
        {
            // cache probabilities
            _corr_prob.get_probs(_probs);

            if (_probs.empty())
            {
                std::unordered_set<deg_t, boost::hash<deg_t> > deg_set;
                for (size_t ei = 0; ei < base_t::_edges.size(); ++ei)
                {
                    edge_t& e = base_t::_edges[ei];
                    deg_set.insert(get_deg(source(e, g), g));
                    deg_set.insert(get_deg(target(e, g), g));
                }

                for (auto s_iter = deg_set.begin(); s_iter != deg_set.end(); ++s_iter)
                    for (auto t_iter = deg_set.begin(); t_iter != deg_set.end(); ++t_iter)
                    {
                        double p = _corr_prob(*s_iter, *t_iter);
                        _probs[make_pair(*s_iter, *t_iter)] = p;
                    }
            }

            for (auto iter = _probs.begin(); iter != _probs.end(); ++iter)
            {
                double& p = iter->second;
                p = log(p);
            }
        }
    }

    double get_prob(const deg_t& s_deg, const deg_t& t_deg)
    {
        if (_probs.empty())
        {
            double p = _corr_prob(s_deg, t_deg);
            // avoid zero probability to not get stuck in rejection step
            if (std::isnan(p) || std::isinf(p) || p <= 0)
                p = numeric_limits<double>::min();
            return log(p);
        }
        auto k = make_pair(s_deg, t_deg);
        auto iter = _probs.find(k);
        if (iter == _probs.end())
            return log(numeric_limits<double>::min());
        return iter->second;
    }

    deg_t get_deg(vertex_t v, Graph& g)
    {
        return _blockdeg.get_block(v, g);
    }

    pair<size_t, bool> get_target_edge(pair<size_t, bool>& e,
                                       bool parallel_edges)
    {
        if (!is_directed::apply<Graph>::type::value)
        {
            std::bernoulli_distribution coin(0.5);
            e.second = coin(base_t::_rng);
        }

        deg_t s_deg = get_deg(source(e, base_t::_edges, _g), _g);
        deg_t t_deg = get_deg(target(e, base_t::_edges, _g), _g);

        std::uniform_int_distribution<> sample(0, base_t::_edges.size() - 1);
        size_t epi = sample(base_t::_rng);
        pair<size_t, bool> ep = make_pair(epi, false);
        if (!is_directed::apply<Graph>::type::value)
        {
            // for undirected graphs we must select a random direction
            std::bernoulli_distribution coin(0.5);
            ep.second = coin(base_t::_rng);
        }

        if (source(e, base_t::_edges, _g) == source(ep, base_t::_edges, _g) ||
            target(e, base_t::_edges, _g) == target(ep, base_t::_edges, _g))
            return ep; // rewiring is a no-op

        deg_t ep_s_deg = get_deg(source(ep, base_t::_edges, _g), _g);
        deg_t ep_t_deg = get_deg(target(ep, base_t::_edges, _g), _g);

        double pi = get_prob(s_deg, t_deg) + get_prob(ep_s_deg, ep_t_deg);
        double pf = get_prob(s_deg, ep_t_deg) + get_prob(ep_s_deg, t_deg);

        if (pf >= pi)
            return ep;

        double a = exp(pf - pi);

        std::uniform_real_distribution<> rsample(0.0, 1.0);
        double r = rsample(base_t::_rng);
        if (r > a)
            return e; // reject
        else
            return ep;
    }

    void update_edge(size_t, bool) {}

private:
    Graph& _g;
    EdgeIndexMap _edge_index;
    CorrProb _corr_prob;
    BlockDeg _blockdeg;


#ifdef HAVE_SPARSEHASH
    typedef google::dense_hash_map<pair<deg_t, deg_t>, double,
                                   boost::hash<pair<deg_t, deg_t>>> prob_map_t;
#else
    typedef std::unordered_map<pair<deg_t, deg_t>, double,
                               boost::hash<pair<deg_t, deg_t>>> prob_map_t;
#endif

    prob_map_t _probs;
};


// general "degree-corrected" stochastic blockmodel
// this version is based on the alias method
template <class Graph, class EdgeIndexMap, class CorrProb, class BlockDeg>
class AliasProbabilisticRewireStrategy:
    public RewireStrategyBase<Graph, EdgeIndexMap,
                              AliasProbabilisticRewireStrategy<Graph, EdgeIndexMap,
                                                               CorrProb, BlockDeg> >
{
public:
    typedef RewireStrategyBase<Graph, EdgeIndexMap,
                               AliasProbabilisticRewireStrategy<Graph, EdgeIndexMap,
                                                                CorrProb, BlockDeg> >
        base_t;

    typedef Graph graph_t;
    typedef EdgeIndexMap edge_index_t;

    typedef typename BlockDeg::block_t deg_t;

    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    typedef typename graph_traits<Graph>::edge_descriptor edge_t;
    typedef typename EdgeIndexMap::value_type index_t;

    AliasProbabilisticRewireStrategy(Graph& g, EdgeIndexMap edge_index,
                                     vector<edge_t>& edges, CorrProb corr_prob,
                                     BlockDeg blockdeg, bool, rng_t& rng,
                                     bool parallel_edges)
        : base_t(g, edge_index, edges, rng, parallel_edges), _g(g),
          _corr_prob(corr_prob), _blockdeg(blockdeg)
    {

#ifdef HAVE_SPARSEHASH
        _probs.set_empty_key(get_null_key<pair<deg_t, deg_t>>()());
        _sampler.set_empty_key(get_null_key<deg_t>()());
        _in_edges.set_empty_key(get_null_key<deg_t>()());
        _out_edges.set_empty_key(get_null_key<deg_t>()());
        _sprob.set_empty_key(get_null_key<deg_t>()());
#endif

        _in_pos.resize(base_t::_edges.size());
        if (!is_directed::apply<Graph>::type::value)
            _out_pos.resize(base_t::_edges.size());

        std::unordered_set<deg_t> deg_set;
        for (size_t ei = 0; ei < base_t::_edges.size(); ++ei)
        {
            edge_t& e = base_t::_edges[ei];
            deg_set.insert(get_deg(source(e, g), g));
            deg_set.insert(get_deg(target(e, g), g));

            vertex_t v = target(e, g);
            deg_t r = get_deg(v, g);
            _in_edges[r].push_back(ei);
            _in_pos[ei] = _in_edges[r].size() - 1;

            if (!is_directed::apply<Graph>::type::value)
            {
                v = source(e, g);
                deg_t s = get_deg(v, g);
                _out_edges[s].push_back(ei);
                _out_pos[ei] = _out_edges[s].size() - 1;
            }
        }

        _corr_prob.get_probs(_probs);

        if (_probs.empty())
        {
            vector<deg_t> items;

            for (auto s_iter = deg_set.begin(); s_iter != deg_set.end(); ++s_iter)
                items.push_back(*s_iter);

            for (auto s_iter = deg_set.begin(); s_iter != deg_set.end(); ++s_iter)
            {
                vector<double> probs;
                double sum = 0;
                for (auto t_iter = deg_set.begin(); t_iter != deg_set.end(); ++t_iter)
                {
                    double p = _corr_prob(*s_iter, *t_iter);
                    // avoid zero probability to not get stuck in rejection step
                    if (std::isnan(p) || std::isinf(p) || p <= 0)
                        continue;
                    probs.push_back(p);
                    _probs[make_pair(*s_iter, *t_iter)] = log(p);
                    sum += p;
                }

                _sampler[*s_iter] = new Sampler<deg_t, boost::mpl::false_>(items, probs);

                auto& ps = _sprob[*s_iter];

#ifdef HAVE_SPARSEHASH
                ps.set_empty_key(get_null_key<deg_t>()());
#endif
                for (size_t i = 0; i < items.size(); ++i)
                {
                    double er = 0;
                    if (!is_directed::apply<Graph>::type::value)
                        er = max(_in_edges[items[i]].size() + _out_edges[items[i]].size(),
                                 numeric_limits<double>::min());
                    else
                        er = max(_in_edges[items[i]].size(),
                                 numeric_limits<double>::min());
                    ps[items[i]] = exp(log(probs[i]) - log(sum) - log(er));
                }
            }
        }
        else
        {
            std::unordered_map<deg_t, vector<double>> sprobs;
            std::unordered_map<deg_t, vector<deg_t>> sitems;
            for (auto iter = _probs.begin(); iter != _probs.end(); ++iter)
            {
                deg_t s = iter->first.first;
                deg_t t = iter->first.second;
                double& p = iter->second;
                if (std::isnan(p) || std::isinf(p) || p <= 0)
                    p = numeric_limits<double>::min();
                sitems[s].push_back(t);
                sprobs[s].push_back(p);
                p = log(p);
            }


            for (auto iter = sitems.begin(); iter != sitems.end(); ++iter)
            {
                deg_t s = iter->first;
                _sampler[s] = new Sampler<deg_t, boost::mpl::false_>(iter->second, sprobs[s]);

                double sum = 0;
                for (size_t i = 0; i < iter->second.size(); ++i)
                {
                    deg_t t = iter->second[i];
                    sum += sprobs[s][i];
                }

                auto& pr = _sprob[s];

#ifdef HAVE_SPARSEHASH
                pr.set_empty_key(get_null_key<deg_t>()());
#endif

                for (size_t i = 0; i < iter->second.size(); ++i)
                {
                    deg_t t = iter->second[i];
                    pr[t] = exp(log(sprobs[s][i]) - log(sum));
                }
            }
        }
    }

    ~AliasProbabilisticRewireStrategy()
    {
        for (typeof(_sampler.begin()) iter = _sampler.begin();
             iter != _sampler.end(); ++iter)
            delete iter->second;
    }

    double get_prob(const deg_t& s_deg, const deg_t& t_deg)
    {
        static const double zero = log(numeric_limits<double>::min());
        auto k = make_pair(s_deg, t_deg);
        auto iter = _probs.find(k);
        if (iter == _probs.end())
            return zero;
        return iter->second;
    }

    double get_sprob(const deg_t& s_deg, const deg_t& t_deg)
    {
        auto& pr = _sprob[s_deg];
        auto iter = pr.find(t_deg);
        if (iter == pr.end())
            return numeric_limits<double>::min();
        return iter->second;
    }

    deg_t get_deg(vertex_t v, Graph& g)
    {
        return _blockdeg.get_block(v, g);
    }

    pair<size_t, bool> get_target_edge(pair<size_t, bool>& e,
                                       bool parallel_edges)
    {
        if (!is_directed::apply<Graph>::type::value)
        {
            std::bernoulli_distribution coin(0.5);
            e.second = coin(base_t::_rng);
        }

        vertex_t s = source(e, base_t::_edges, _g);
        vertex_t t = target(e, base_t::_edges, _g);

        deg_t s_deg = get_deg(s, _g);
        deg_t t_deg = get_deg(t, _g);

        auto iter = _sampler.find(s_deg);
        if (iter == _sampler.end())
            throw GraphException("Block label without defined connection probability!");

        deg_t nt = iter->second->sample(base_t::_rng);

        if (_in_edges[nt].empty() && _out_edges[nt].empty())
            return e; // reject

        pair<size_t, bool> ep;
        std::bernoulli_distribution coin(_in_edges[nt].size() /
                                         double(_in_edges[nt].size() +
                                                _out_edges[nt].size()));
        if (is_directed::apply<Graph>::type::value || coin(base_t::_rng))
        {
            vector<size_t>& ies = _in_edges[nt];

            std::uniform_int_distribution<> sample(0, ies.size() - 1);
            size_t epi = ies[sample(base_t::_rng)];

            ep = make_pair(epi, false);
        }
        else
        {
            vector<size_t>& oes = _out_edges[nt];

            std::uniform_int_distribution<> sample(0, oes.size() - 1);
            size_t epi = oes[sample(base_t::_rng)];

            ep = make_pair(epi, true);
        }

        if (source(e, base_t::_edges, _g) == source(ep, base_t::_edges, _g) ||
            target(e, base_t::_edges, _g) == target(ep, base_t::_edges, _g))
            return ep; // rewiring is a no-op

        vertex_t ep_s = source(ep, base_t::_edges, _g);
        vertex_t ep_t = target(ep, base_t::_edges, _g);

        deg_t ep_s_deg = get_deg(ep_s, _g);
        deg_t ep_t_deg = get_deg(ep_t, _g);

        //assert(ep_t_deg == nt);

        double pi = get_prob(s_deg, t_deg) + get_prob(ep_s_deg, ep_t_deg);
        double pf = get_prob(s_deg, ep_t_deg) + get_prob(ep_s_deg, t_deg);

        double a = exp(pf - pi);

        if (is_directed::apply<Graph>::type::value)
        {
            double e_ep_t = _in_edges[ep_t_deg].size();
            double e_t = _in_edges[t_deg].size();
            if (t == ep_t)
            {
                e_ep_t -= in_degreeS()(t, _g);
                e_t -= in_degreeS()(ep_t, _g);
            }

            a /= (get_sprob(s_deg, ep_t_deg) / e_ep_t +
                  get_sprob(ep_s_deg, t_deg) / e_t);        // forwards
            a *= (get_sprob(s_deg, t_deg) / e_t +
                  get_sprob(ep_s_deg, ep_t_deg) / e_ep_t);  // backwards
        }
        else
        {
            double e_ep_t = _in_edges[ep_t_deg].size() + _out_edges[ep_t_deg].size();
            double e_t = _in_edges[t_deg].size()  + _out_edges[t_deg].size();
            if (t == ep_t)
            {
                e_ep_t -= out_degree(t, _g);
                e_t -= out_degree(ep_t, _g);
            }

            double e_ep_s = _in_edges[ep_s_deg].size() + _out_edges[ep_s_deg].size();
            double e_s = _in_edges[s_deg].size() + _out_edges[s_deg].size();
            if (s == ep_s)
            {
                e_ep_s -= out_degree(s, _g);
                e_s -= out_degree(ep_s, _g);
            }

            a /= (get_sprob(s_deg, ep_t_deg) / e_ep_t +
                  get_sprob(ep_t_deg, s_deg) / e_s +
                  get_sprob(ep_s_deg, t_deg) / e_t +
                  get_sprob(t_deg, ep_s_deg) / e_ep_s);      // forwards
            a *= (get_sprob(s_deg, t_deg) / e_t +
                  get_sprob(t_deg, s_deg) / e_s +
                  get_sprob(ep_s_deg, ep_t_deg) / e_ep_t +
                  get_sprob(ep_t_deg, ep_s_deg) / e_ep_s);   // backwards
        }

        if (a > 1)
            return ep;

        std::uniform_real_distribution<> rsample(0.0, 1.0);
        double r = rsample(base_t::_rng);
        if (r > a)
            return e; // reject
        return ep;
    }

    void update_edge(size_t ei, bool insert)
    {
        if (insert)
        {
            vertex_t v = target(base_t::_edges[ei], _g);
            deg_t d = get_deg(v, _g);
            auto& in_vec = _in_edges[d];
            in_vec.push_back(ei);
            _in_pos[ei] = in_vec.size() - 1;

            if (!is_directed::apply<Graph>::type::value)
            {
                v = source(base_t::_edges[ei], _g);
                deg_t d = get_deg(v, _g);
                auto& out_vec = _out_edges[d];
                out_vec.push_back(ei);
                _out_pos[ei] = out_vec.size() - 1;
            }
        }
        else
        {
            vertex_t v = target(base_t::_edges[ei], _g);
            deg_t d = get_deg(v, _g);
            size_t j = _in_pos[ei];
            auto& in_vec = _in_edges[d];
            _in_pos[in_vec.back()] = j;
            in_vec[j] = in_vec.back();
            in_vec.pop_back();

            if (!is_directed::apply<Graph>::type::value)
            {
                v = source(base_t::_edges[ei], _g);
                deg_t d = get_deg(v, _g);
                size_t j = _out_pos[ei];
                auto& out_vec = _out_edges[d];
                _out_pos[out_vec.back()] = j;
                out_vec[j] = out_vec.back();
                out_vec.pop_back();
            }
        }
    }

private:
    Graph& _g;
    EdgeIndexMap _edge_index;
    CorrProb _corr_prob;
    BlockDeg _blockdeg;

#ifdef HAVE_SPARSEHASH
    typedef google::dense_hash_map<deg_t, Sampler<deg_t, boost::mpl::false_>*> sampler_map_t;
#else
    typedef std::unordered_map<deg_t, Sampler<deg_t, boost::mpl::false_>*> sampler_map_t;
#endif

    sampler_map_t _sampler;

#ifdef HAVE_SPARSEHASH
    typedef google::dense_hash_map<deg_t, google::dense_hash_map<deg_t, double>> sprob_map_t;
#else
    typedef std::unordered_map<deg_t, std::unordered_map<deg_t, double>> sprob_map_t;
#endif

    sprob_map_t _sprob;


#ifdef HAVE_SPARSEHASH
    typedef google::dense_hash_map<pair<deg_t, deg_t>, double, boost::hash<pair<deg_t, deg_t>>> prob_map_t;
#else
    typedef std::unordered_map<pair<deg_t, deg_t>, double, boost::hash<pair<deg_t, deg_t>>> prob_map_t;
#endif

    prob_map_t _probs;

#ifdef HAVE_SPARSEHASH
    typedef google::dense_hash_map<deg_t, vector<size_t>> edge_map_t;
#else
    typedef std::unordered_map<deg_t, vector<size_t>> edge_map_t;
#endif

    edge_map_t _in_edges;
    edge_map_t _out_edges;
    vector<size_t> _in_pos;
    vector<size_t> _out_pos;
};


// general "traditional" stochastic blockmodel
// this version is based on the alias method, and does not keep the degrees fixed
template <class Graph, class EdgeIndexMap, class CorrProb, class BlockDeg>
class TradBlockRewireStrategy
{
public:
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    typedef typename graph_traits<Graph>::edge_descriptor edge_t;
    typedef typename EdgeIndexMap::value_type index_t;
    typedef typename BlockDeg::block_t deg_t;

    TradBlockRewireStrategy(Graph& g, EdgeIndexMap edge_index,
                            vector<edge_t>& edges, CorrProb corr_prob,
                            BlockDeg blockdeg, bool, rng_t& rng,
                            bool parallel_edges)

        : _g(g), _edge_index(edge_index), _edges(edges), _corr_prob(corr_prob),
          _blockdeg(blockdeg), _rng(rng), _sampler(nullptr)
    {

#ifdef HAVE_SPARSEHASH
        _vertices.set_empty_key(get_null_key<deg_t>()());
#endif

        typename graph_traits<Graph>::vertex_iterator v, v_end;
        for (tie(v, v_end) = vertices(_g); v != v_end; ++v)
        {
            deg_t d = _blockdeg.get_block(*v, g);
            _vertices[d].push_back(*v);
        }

        std::unordered_map<pair<deg_t, deg_t>, double, boost::hash<pair<deg_t, deg_t> > >
            probs;
        _corr_prob.get_probs(probs);

        vector<double> dprobs;
        if (probs.empty())
        {
            for (auto s_iter = _vertices.begin(); s_iter != _vertices.end(); ++s_iter)
            {
                for (auto t_iter = _vertices.begin(); t_iter != _vertices.end(); ++t_iter)
                {
                    double p = _corr_prob(s_iter->first, t_iter->first);
                    if (std::isnan(p) || std::isinf(p) || p <= 0)
                        continue;

                    _items.push_back(make_pair(s_iter->first, t_iter->first));
                    dprobs.push_back(p);
                }
            }
        }
        else
        {
            for (auto iter = probs.begin(); iter != probs.end(); ++iter)
            {
                deg_t s = iter->first.first;
                deg_t t = iter->first.second;
                double p = iter->second;
                // avoid zero probability to not get stuck in rejection step
                if (std::isnan(p) || std::isinf(p) || p <= 0)
                    continue;
                _items.push_back(make_pair(s, t));
                dprobs.push_back(p);
            }
        }

        if (_items.empty())
            throw GraphException("No connection probabilities larger than zero!");

        _sampler = new Sampler<pair<deg_t, deg_t> >(_items, dprobs);
    }

    ~TradBlockRewireStrategy()
    {
        if (_sampler != nullptr)
            delete _sampler;
    }

    bool operator()(size_t ei, bool self_loops, bool parallel_edges)
    {
        typename graph_traits<Graph>::vertex_descriptor s, t;

        while (true)
        {
            const pair<deg_t, deg_t>& deg = _sampler->sample(_rng);

            vector<vertex_t>& svs = _vertices[deg.first];
            vector<vertex_t>& tvs = _vertices[deg.second];

            if (svs.empty() || tvs.empty())
                continue;

            std::uniform_int_distribution<size_t> s_sample(0, svs.size() - 1);
            std::uniform_int_distribution<size_t> t_sample(0, tvs.size() - 1);

            s = svs[s_sample(_rng)];
            t = tvs[t_sample(_rng)];

            break;
        }

        // reject self-loops if not allowed
        if (!self_loops && s == t)
            return false;

        // reject parallel edges if not allowed
        if (!parallel_edges && is_adjacent(s, t, _g))
            return false;

        remove_edge(_edges[ei], _g);
        edge_t ne = add_edge(s, t, _g).first;
        _edges[ei] = ne;

        return true;
    }

private:
    Graph& _g;
    EdgeIndexMap _edge_index;
    vector<edge_t>& _edges;
    CorrProb _corr_prob;
    BlockDeg _blockdeg;
    rng_t& _rng;

#ifdef HAVE_SPARSEHASH
    google::dense_hash_map<deg_t, vector<vertex_t>> _vertices;
#else
    std::unordered_map<deg_t, vector<vertex_t>> _vertices;
#endif

    vector<pair<deg_t, deg_t> > _items;
    Sampler<pair<deg_t, deg_t> >* _sampler;

};

} // graph_tool namespace

#endif // GRAPH_REWIRING_HH
