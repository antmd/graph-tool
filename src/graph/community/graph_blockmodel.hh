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

#ifndef GRAPH_BLOCKMODEL_HH
#define GRAPH_BLOCKMODEL_HH

#include <cmath>
#include <iostream>
#include <queue>
#include <boost/math/special_functions/zeta.hpp>
#include <boost/functional/hash.hpp>

#include "config.h"
#include <unordered_set>
#include <unordered_map>
#include <tuple>

#ifdef HAVE_SPARSEHASH
#include SPARSEHASH_INCLUDE(dense_hash_set)
#include SPARSEHASH_INCLUDE(dense_hash_map)
#endif

#include "../generation/sampler.hh"
#include "../generation/dynamic_sampler.hh"

#ifdef USING_OPENMP
#include <omp.h>
#endif


double spence(double);

namespace graph_tool
{

#ifdef HAVE_SPARSEHASH
using google::dense_hash_set;
using google::dense_hash_map;
#endif

using namespace boost;

template <class Key>
class IdentityArrayPropertyMap
    : public boost::put_get_helper<Key, IdentityArrayPropertyMap<Key>>
{
public:
    typedef std::array<Key, 1> value_type;
    typedef value_type reference;
    typedef Key key_type;
    typedef boost::readable_property_map_tag category;

    inline __attribute__((always_inline))
    const value_type operator[](const key_type& c) const { return {c}; }
};

// ====================
// Entropy calculation
// ====================

// Repeated computation of x*log(x) and log(x) actually adds up to a lot of
// time. A significant speedup can be made by caching pre-computed values. This
// is doable since the values of mrse are bounded in [0, 2E], where E is the
// total number of edges in the network. Here we simply grow the cache as
// needed.

extern vector<double> __safelog_cache;
extern vector<double> __xlogx_cache;
extern vector<double> __lgamma_cache;

template <class Type>
__attribute__((always_inline))
inline double safelog(Type x)
{
    if (x == 0)
        return 0;
    return log(x);
}

__attribute__((always_inline))
inline double safelog(size_t x)
{
    if (x >= __safelog_cache.size())
    {
        if (x == 0)
            return 0;
        return log(x);
    }
    return __safelog_cache[x];
}

__attribute__((always_inline))
inline double xlogx(size_t x)
{
    if (x >= __xlogx_cache.size())
        return x * safelog(x);
    return __xlogx_cache[x];
}

__attribute__((always_inline))
inline double lgamma_fast(size_t x)
{
    if (x >= __lgamma_cache.size())
        return lgamma(x);
    return __lgamma_cache[x];
}

// polylogarithm and degree-distribution description length (xi)

template <class Type>
Type polylog(int n, Type z, Type epsilon=1e-6)
{
    if (n == 2)
        return spence(1 - z);

    int k = 1;
    Type S = 0;
    Type delta = epsilon + 1;
    Type zk = z;
    while (delta > epsilon)
    {
        Type dS = zk / pow(k, n);
        k++;
        zk *= z;
        S += dS;
        delta = dS;
    }
    return S;
}

template <class NType, class EType>
void get_mu_l(NType N, EType E, double& mu, double& l,
              double epsilon=1e-8)
{
    mu = sqrt(polylog<double>(2, 1.) / double(E));
    l = 1. - exp(-double(N) * mu);

    double delta = epsilon + 1;
    while (delta > epsilon)
    {
        double nmu = sqrt(polylog<double>(2, l) / double(E));
        double nl = 1. - exp(-double(N) * mu);

        delta = abs(nmu - mu) + abs(nl - l);
        mu = nmu;
        l = nl;
    }

    l = -log(l);
}

template <class NType, class EType>
double get_xi(NType N, EType E, double epsilon=1e-8)
{
    if (E == 0 || N == 0)
        return 0;

    double mu = 0, l = 0;
    get_mu_l(N, E, mu, l, epsilon);
    double S = double(N) * l + 2 * double(E) * mu;
    return S;
}

template <class NType, class EType>
double get_xi_fast(NType N, EType E)
{
    if (E == 0 || N == 0)
        return 0;
    static const double z2 = boost::math::zeta(2.);
    return 2 * sqrt(z2 * E);
}

// Sparse entropy terms
// ====================

// "edge" term of the entropy
template <class Graph>
__attribute__((always_inline))
inline double eterm(size_t r, size_t s, size_t mrs, const Graph&)
{
    if (!is_directed::apply<Graph>::type::value && r == s)
        mrs *= 2;

    double val = xlogx(mrs);

    if (is_directed::apply<Graph>::type::value || r != s)
        return -val;
    else
        return -val / 2;
}

// "vertex" term of the entropy
template <class Graph>
inline double vterm(size_t mrp, size_t mrm, size_t wr, bool deg_corr,
                    Graph&)
{
    double one = 0.5;

    if (is_directed::apply<Graph>::type::value)
        one = 1;

    if (deg_corr)
        return one * (xlogx(mrm) + xlogx(mrp));
    else
        return one * (mrm * safelog(wr) + mrp * safelog(wr));
}

struct entropy
{
    template <class Graph, class Eprop, class Vprop>
    void operator()(Eprop& mrs, Vprop& mrp, Vprop& mrm, Vprop& wr,
                    bool deg_corr, Graph& g, double& S) const
    {
        S = 0;
        for (auto e : edges_range(g))
            S += eterm(source(e, g), target(e, g), mrs[e], g);
        for (auto v : vertices_range(g))
            S += vterm(mrp[v], mrm[v], wr[v], deg_corr, g);
    }
};

// Parallel edges
// ==============

template <class Vertex, class List, class Graph, class GetNode>
double get_parallel_neighbours_entropy(Vertex v, List& us, Graph&,
                                       const GetNode& get_node)
{
    double S = 0;
    for (auto& uc : us)
    {
        auto& u = uc.first;
        auto& m = uc.second;
        if (m > 1)
        {
            if (get_node(u) == size_t(v) && !is_directed::apply<Graph>::type::value)
            {
                assert(m % 2 == 0);
                S += lgamma_fast(m/2 + 1);
            }
            else
            {
                S += lgamma_fast(m + 1);
            }
        }
    }
    return S;
}

struct entropy_parallel_edges
{
    template <class Graph, class Weight>
    void operator()(Graph& g, Weight weight, double& S) const
    {
        S = 0;
        auto get_node = [](size_t i) {return i;};
        for (auto v : vertices_range(g))
        {
            unordered_map<decltype(v), int> us;
            for (auto e : out_edges_range(v, g))
            {
                auto u = target(e, g);
                if (u < v && !is_directed::apply<Graph>::type::value)
                    continue;
                us[u] += weight[e];
            }

            S += get_parallel_neighbours_entropy(v, us, g, get_node);
        }
    }
};


// Dense entropy
// =============

// Warning: lgamma(x) is not thread-safe! However, since in the context of this
// program the outputs should _always_ be positive, this can be overlooked.

__attribute__((always_inline))
inline double lbinom(double N, double k)
{
    if (N == 0 || k == 0 || k > N)
        return 0;
    return (lgamma(N + 1) - lgamma(k + 1)) - lgamma(N - k + 1);
}

__attribute__((always_inline))
inline double lbinom_fast(int N, int k)
{
    if (N == 0 || k == 0 || k > N)
        return 0;
    return lgamma_fast(N + 1) - lgamma_fast(N - k + 1) - lgamma_fast(k + 1);
}

__attribute__((always_inline))
inline double lbinom_careful(double N, double k)
{
    if (N == 0 || k == 0 || k >= N)
        return 0;
    double lgN = lgamma(N + 1);
    double lgk = lgamma(k + 1);
    if (lgN - lgk > 1e8)
    {
        // We have N >> k. Use Stirling's approximation: ln N! ~ N ln N - N
        // and reorder
        return - N * log1p(-k / N) - k * log1p(-k / N) - k - lgk + k * log(N);
    }
    else
    {
        return lgN - lgamma(N - k + 1) - lgk;
    }
}


// "edge" term of the entropy
template <class Graph>
__attribute__((always_inline))
inline double eterm_dense(size_t r, size_t s, int ers, double wr_r,
                          double wr_s, bool multigraph, const Graph&)
{
    // we should not use integers here, since they may overflow
    double nrns;

    if (ers == 0)
        return 0.;

    if (r != s || is_directed::apply<Graph>::type::value)
    {
        nrns = wr_r * wr_s;
    }
    else
    {
        if (multigraph)
            nrns = (wr_r * (wr_r + 1)) / 2;
        else
            nrns = (wr_r * (wr_r - 1)) / 2;
    }

    double S;
    if (multigraph)
        S = lbinom(nrns + ers - 1, ers); // do not use lbinom_fast!
    else
        S = lbinom(nrns, ers);
    return S;
}

struct entropy_dense
{
    template <class Graph, class Eprop, class Vprop>
    void operator()(Eprop mrs, Vprop& wr, bool multigraph, Graph& g, double& S) const
    {
        S = 0;
        for (auto e : edges_range(g))
        {
            auto r = source(e, g);
            auto s = target(e, g);
            S += eterm_dense(r, s, mrs[e], wr[r], wr[s], multigraph, g);
        }
    }
};

// ===============
// Partition stats
// ===============

class partition_stats_t
{
public:

#ifdef HAVE_SPARSEHASH
    typedef dense_hash_map<pair<size_t,size_t>, int, std::hash<pair<size_t,size_t>>> map_t;
#else
    typedef unordered_map<pair<size_t,size_t>, int> map_t;
#endif

    partition_stats_t() : _enabled(false) {}

    template <class Graph, class Vprop, class Eprop>
    partition_stats_t(Graph& g, Vprop b, Eprop eweight, size_t N, size_t B,
                      bool edges_dl)
        : _enabled(true), _N(N), _E(0), _B(B), _hist(B), _total(B), _ep(B),
          _em(B), _edges_dl(edges_dl)
    {

#ifdef HAVE_SPARSEHASH
        for (size_t r = 0; r < B; ++r)
        {
            _hist[r].set_empty_key(make_pair(numeric_limits<size_t>::max(),
                                             numeric_limits<size_t>::max()));
            _hist[r].set_deleted_key(make_pair(numeric_limits<size_t>::max() - 1,
                                               numeric_limits<size_t>::max() - 1));
        }
#endif

        for (auto v : vertices_range(g))
        {
            auto r = b[v];
            size_t kin = in_degreeS()(v, g, eweight);
            size_t kout = out_degreeS()(v, g, eweight);
            _hist[r][make_pair(kin, kout)]++;
            _total[r]++;
            _em[r] += kin;
            _ep[r] += kout;
            _E += kout;
        }

        if (!is_directed::apply<Graph>::type::value)
            _E /= 2;

        _actual_B = 0;
        for (size_t r = 0; r < B; ++r)
            if (_total[r] > 0)
                _actual_B++;
    }

    double get_partition_dl()
    {
        double S = 0;
        S += lbinom(_actual_B + _N - 1, _N);
        S += lgamma(_N + 1);
        for (auto nr : _total)
            S -= lgamma(nr + 1);
        return S;
    }

    double get_deg_dl(bool ent, bool dl_alt, bool xi_fast)
    {
        double S = 0;
        for (size_t r = 0; r < _B; ++r)
        {
            if (ent)
            {
                for (auto& k_c : _hist[r])
                {
                    double p = k_c.second / double(_total[r]);
                    S -= p * log(p) * _total[r];
                }
            }
            else
            {
                double S1 = 0;

                if (xi_fast)
                {
                    S1 += get_xi_fast(_total[r], _ep[r]);
                    S1 += get_xi_fast(_total[r], _em[r]);
                }
                else
                {
                    S1 += get_xi(_total[r], _ep[r]);
                    S1 += get_xi(_total[r], _em[r]);
                }

                S1 += lgamma(_total[r] + 1);
                for (auto& k_c : _hist[r])
                    S1 -= lgamma(k_c.second + 1);

                if (dl_alt)
                {
                    double S2 = 0;
                    S2 += lbinom(_total[r] + _ep[r] - 1, _ep[r]);
                    S2 += lbinom(_total[r] + _em[r] - 1, _em[r]);
                    S += min(S1, S2);
                }
                else
                {
                    S += S1;
                }
            }
        }
        return S;
    }

    template <class Graph, class OStats>
    double get_delta_dl(size_t, size_t r, size_t nr, OStats&, Graph&)
    {
        if (r == nr)
            return 0;

        double S_b = 0, S_a = 0;
        S_b += -lgamma_fast(_total[r] + 1) - lgamma_fast(_total[nr] + 1);
        S_a += -lgamma_fast(_total[r]    ) - lgamma_fast(_total[nr] + 2);

        int dB = 0;
        if (_total[r] == 1)
            dB--;
        if (_total[nr] == 0)
            dB++;

        if (dB != 0)
        {
            S_b += lbinom(_actual_B + _N - 1, _N);
            S_a += lbinom(_actual_B + dB + _N - 1, _N);

            if (_edges_dl)
            {
                auto get_x = [](size_t B) -> size_t
                {
                    if (is_directed::apply<Graph>::type::value)
                        return B * B;
                    else
                        return (B * (B + 1)) / 2;
                };

                S_b += lbinom(get_x(_actual_B) + _E - 1, _E);
                S_a += lbinom(get_x(_actual_B + dB) + _E - 1, _E);
            }
        }

        return S_a - S_b;
    }

    template <class Graph, class EWeight, class OStats>
    double get_delta_deg_dl(size_t v, size_t r, size_t nr, EWeight& eweight,
                            OStats&, Graph& g)
    {
        if (r == nr)
            return 0;

        double S_b = 0, S_a = 0;

        int kin = in_degreeS()(v, g, eweight);
        int kout = out_degreeS()(v, g, eweight);

        auto get_Se = [&](size_t s, int delta, int kin, int kout) -> double
            {
                double S = 0;
                S += get_xi_fast(_total[s] + delta, _em[s] + kin);
                S += get_xi_fast(_total[s] + delta, _ep[s] + kout);
                return S;
            };

        S_b += get_Se(r,  0,    0,     0) + get_Se(nr, 0,   0,    0);
        S_a += get_Se(r, -1, -kin, -kout) + get_Se(nr, 1, kin, kout);

        auto get_Sr = [&](size_t s, int delta) -> double
            {
                return lgamma_fast(_total[s] + delta + 1);
            };

        S_b += get_Sr(r,  0) + get_Sr(nr, 0);
        S_a += get_Sr(r, -1) + get_Sr(nr, 1);


        auto get_Sk = [&](size_t s, pair<size_t, size_t>& deg, int delta) -> double
            {
                size_t nd = 0;
                auto iter = _hist[s].find(deg);
                if (iter != _hist[s].end())
                    nd = iter->second;

                return -lgamma_fast(nd + delta + 1);
            };

        auto deg = make_pair(size_t(kin), size_t(kout));
        S_b += get_Sk(r, deg,  0) + get_Sk(nr, deg, 0);
        S_a += get_Sk(r, deg, -1) + get_Sk(nr, deg, 1);

        return S_a - S_b;
    }

    template <class Graph, class OStats, class EWeight>
    void move_vertex(size_t v, size_t r, size_t nr, bool deg_corr, OStats&,
                     Graph& g, EWeight& eweight, size_t kin = 0,
                     size_t kout = 0)
    {
        if (r == nr)
            return;

        _total[r]--;
        _total[nr]++;

        if (_total[r] == 0)
            _actual_B--;
        if (_total[nr] == 1)
            _actual_B++;

        if (deg_corr)
        {
            if (kin + kout == 0)
            {
                kin = in_degreeS()(v, g, eweight);
                kout = out_degreeS()(v, g, eweight);
            }
            auto deg = make_pair(kin, kout);
            auto iter = _hist[r].find(deg);
            iter->second--;
            if (iter->second == 0)
                _hist[r].erase(iter);
            _hist[nr][deg]++;
            _em[r] -= deg.first;
            _ep[r] -= deg.second;
            _em[nr] += deg.first;
            _ep[nr] += deg.second;
        }
    }

    bool is_enabled() { return _enabled; }
    void set_enabled(bool enabled) { _enabled = enabled; }

private:
    bool _enabled;
    size_t _N;
    size_t _E;
    size_t _B;
    size_t _actual_B;
    vector<map_t> _hist;
    vector<int> _total;
    vector<int> _ep;
    vector<int> _em;
    bool _edges_dl;
};

// ===============================
// Block moves
// ===============================

// this structure speeds up the access to the edges between given blocks,
// since we're using an adjacency list to store the block structure (the emat_t
// is simply a corresponding adjacency matrix)
struct get_emat_t
{
    template <class Graph>
    struct apply
    {
        typedef multi_array<pair<typename graph_traits<Graph>::edge_descriptor, bool>, 2> type;
    };
};


struct create_emat
{
    template <class Graph>
    void operator()(Graph& g, boost::any& oemap) const
    {
        typedef typename get_emat_t::apply<Graph>::type emat_t;
        size_t B = num_vertices(g);
        emat_t emat(boost::extents[B][B]);

        for (size_t i = 0; i < B; ++i)
            for (size_t j = 0; j < B; ++j)
                emat[i][j].second = false;
        for (auto e : edges_range(g))
        {
            if (source(e, g) >= B || target(e, g) >= B)
                throw GraphException("incorrect number of blocks when creating emat!");
            emat[source(e, g)][target(e, g)] = make_pair(e, true);
            if (!is_directed::apply<Graph>::type::value)
                emat[target(e, g)][source(e, g)] = make_pair(e, true);
        }

        oemap = emat;
    }
};

template <class Graph>
inline __attribute__((always_inline))
pair<typename graph_traits<Graph>::edge_descriptor, bool>
get_me(typename graph_traits<Graph>::vertex_descriptor r,
       typename graph_traits<Graph>::vertex_descriptor s,
       const typename get_emat_t::apply<Graph>::type& emat, const Graph&)
{
    return emat[r][s];
}

template <class Graph>
inline __attribute__((always_inline))
void
put_me(typename graph_traits<Graph>::vertex_descriptor r,
       typename graph_traits<Graph>::vertex_descriptor s,
       const typename graph_traits<Graph>::edge_descriptor& e,
       typename get_emat_t::apply<Graph>::type& emat, const Graph&)
{
    emat[r][s] = make_pair(e, true);
    if (!is_directed::apply<Graph>::type::value && r != s)
        emat[s][r] = make_pair(e, true);
}

template <class Graph>
inline __attribute__((always_inline))
void
remove_me(typename graph_traits<Graph>::vertex_descriptor r,
          typename graph_traits<Graph>::vertex_descriptor s,
          const typename graph_traits<Graph>::edge_descriptor&,
          typename get_emat_t::apply<Graph>::type& emat, Graph&,
          bool delete_edge=true)
{
    if (!delete_edge)
    {
        emat[r][s].second = false;
        if (!is_directed::apply<Graph>::type::value && r != s)
            emat[s][r].second = false;
    }
}


// this structure speeds up the access to the edges between given blocks, since
// we're using an adjacency list to store the block structure (this is like
// emat_t above, but takes less space and is slower)
struct get_ehash_t
{
    template <class Graph>
    struct apply
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;
#ifdef HAVE_SPARSEHASH
        typedef dense_hash_map<vertex_t, edge_t, std::hash<vertex_t>> map_t;
#else
        typedef unordered_map<vertex_t, edge_t> map_t;
#endif
        typedef vector<map_t> type;
    };
};


template<class Graph>
inline __attribute__((always_inline))
pair<typename graph_traits<Graph>::edge_descriptor, bool>
get_me(typename graph_traits<Graph>::vertex_descriptor r,
       typename graph_traits<Graph>::vertex_descriptor s,
       const typename get_ehash_t::apply<Graph>::type& ehash, const Graph&)
{
    assert(r < ehash.size());
    const auto& map = ehash[r];
    auto iter = map.find(s);
    if (iter == map.end())
        return (make_pair(typename graph_traits<Graph>::edge_descriptor(), false));
    return make_pair(iter->second, true);
}

template<class Graph>
inline __attribute__((always_inline))
void
put_me(typename graph_traits<Graph>::vertex_descriptor r,
       typename graph_traits<Graph>::vertex_descriptor s,
       const typename graph_traits<Graph>::edge_descriptor& e,
       typename get_ehash_t::apply<Graph>::type& ehash,
       const Graph&)
{
    assert(r < ehash.size());
    ehash[r][s] = e;
    if (!is_directed::apply<Graph>::type::value)
        ehash[s][r] = e;
}

template<class Graph>
inline __attribute__((always_inline))
void
remove_me(typename graph_traits<Graph>::vertex_descriptor r,
          typename graph_traits<Graph>::vertex_descriptor s,
          const typename graph_traits<Graph>::edge_descriptor& e,
          typename get_ehash_t::apply<Graph>::type& ehash, Graph& bg,
          bool delete_edge=true)
{
    assert(r < ehash.size());
    ehash[r].erase(s);
    if (!is_directed::apply<Graph>::type::value)
        ehash[s].erase(r);
    if (delete_edge)
        remove_edge(e, bg);
}

struct create_ehash
{
    template <class Graph>
    void operator()(Graph& g, boost::any& oemap) const
    {
        typedef typename get_ehash_t::apply<Graph>::type emat_t;
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

        emat_t emat(num_vertices(g));

#ifdef HAVE_SPARSEHASH
        for (auto v : vertices_range(g))
        {
            emat[v].set_empty_key(numeric_limits<vertex_t>::max());
            emat[v].set_deleted_key(numeric_limits<vertex_t>::max() - 1);
        }
#endif

        for (auto e : edges_range(g))
            put_me(source(e, g), target(e, g), e, emat, g);

#ifdef HAVE_SPARSEHASH
        for (auto v : vertices_range(g))
            emat[v].resize(0);
#endif
        oemap = emat;
    }
};

template <class Vertex, class Eprop, class Emat, class BGraph>
__attribute__((always_inline))
inline size_t get_mrs(Vertex r, Vertex s, const Eprop& mrs, Emat& emat,
                      BGraph& bg)
{
    const pair<typename graph_traits<BGraph>::edge_descriptor, bool> me =
        get_me(r, s, emat, bg);
    if (me.second)
        return mrs[me.first];
    else
        return 0;
}

struct standard_neighbours_policy
{
    template <class Graph, class Vertex>
    IterRange<typename out_edge_iteratorS<Graph>::type>
    get_out_edges(Vertex v, Graph& g) const
    {
        return out_edges_range(v, g);
    }

    template <class Graph, class Vertex>
    IterRange<typename in_edge_iteratorS<Graph>::type>
    get_in_edges(Vertex v, Graph& g) const
    {
        return in_edges_range(v, g);
    }

    template <class Graph, class Vertex, class Weight>
    int get_out_degree(Vertex& v, Graph& g, Weight& eweight) const
    {
        return out_degreeS()(v, g, eweight);
    }

    template <class Graph, class Vertex, class Weight>
    int get_in_degree(Vertex& v, Graph& g, Weight& eweight) const
    {
        return in_degreeS()(v, g, eweight);
    }
};

// remove a vertex from its current block
template <class Graph, class BGraph, class Eprop, class Vprop, class EWprop,
          class VWprop, class EMat, class OStats,
          class NPolicy = standard_neighbours_policy>
void remove_vertex(size_t v, Eprop& mrs, Vprop& mrp, Vprop& mrm, Vprop& wr,
                   Vprop& b, const EWprop& eweight, const VWprop& vweight,
                   Graph& g, BGraph& bg, EMat& emat, OStats& overlap_stats,
                   const NPolicy& npolicy = NPolicy())
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    vertex_t r = b[v];

    int self_weight = 0;
    for (auto e : npolicy.get_out_edges(v, g))
    {
        vertex_t u = target(e, g);
        vertex_t s = b[u];

        // if (!emat[r][s].second)
        //     throw GraphException("no edge? " + lexical_cast<string>(r) +
        //                          " " + lexical_cast<string>(s));

        auto me = get_me(r, s, emat, bg).first;

        size_t ew = eweight[e];
        if (u == v && !is_directed::apply<Graph>::type::value)
        {
            self_weight += ew;
        }
        else
        {
            mrs[me] -= ew;

            assert(mrs[me] >= 0);

            mrp[r] -= ew;
            mrm[s] -= ew;

            if (mrs[me] == 0)
                remove_me(r, s, me, emat, bg);
        }
    }

    if (self_weight > 0)
    {
        assert(self_weight % 2 == 0);
        auto me = get_me(r, r, emat, bg).first;
        mrs[me] -= self_weight / 2;
        mrp[r] -= self_weight / 2;
        mrm[r] -= self_weight / 2;
        assert(mrs[me] >= 0);
        if (mrs[me] == 0)
            remove_me(r, r, me, emat, bg);
    }

    for (auto e : npolicy.get_in_edges(v, g))
    {
        vertex_t u = source(e, g);
        if (u == v)
            continue;
        vertex_t s = b[u];

        // if (!emat[s][r].second)
        //     throw GraphException("no edge? " + lexical_cast<string>(s) +
        //                          " " + lexical_cast<string>(r));

        typename graph_traits<BGraph>::edge_descriptor me =
            get_me(s, r, emat, bg).first;

        size_t ew = eweight[e];
        mrs[me] -= ew;

        mrp[s] -= ew;
        mrm[r] -= ew;

        if (mrs[me] == 0)
            remove_me(s, r, me, emat, bg);
    }

    if (!overlap_stats.is_enabled())
    {
        wr[r] -= vweight[v];
    }
    else
    {
        overlap_stats.remove_half_edge(v, r, b, g);
        wr[r] = overlap_stats.get_block_size(r);
    }
}

// add a vertex to block rr
template <class Graph, class BGraph, class Eprop, class Vprop, class EWprop,
          class VWprop, class EMat, class OStats,
          class NPolicy = standard_neighbours_policy>
void add_vertex(size_t v, size_t r, Eprop& mrs, Vprop& mrp, Vprop& mrm,
                Vprop& wr, Vprop& b, const EWprop& eweight,
                const VWprop& vweight, Graph& g, BGraph& bg, EMat& emat,
                OStats& overlap_stats, const NPolicy& npolicy = NPolicy())
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

    int self_weight = 0;
    for (auto e : npolicy.get_out_edges(v, g))
    {
        vertex_t u = target(e, g);
        vertex_t s;

        if (u != v)
            s = b[u];
        else
            s = r;

        typename graph_traits<BGraph>::edge_descriptor me;

        auto mep = get_me(r, s, emat, bg);

        if (!mep.second)
        {
            mep = add_edge(r, s, bg);
            put_me(r, s, mep.first, emat, bg);
            mrs[mep.first] = 0;
        }
        me = mep.first;

        size_t ew = eweight[e];

        if (u == v && !is_directed::apply<Graph>::type::value)
        {
            self_weight += ew;
        }
        else
        {
            mrs[me] += ew;
            mrp[r] += ew;
            mrm[s] += ew;
        }
    }

    if (self_weight > 0)
    {
        assert(self_weight % 2 == 0);
        auto me = get_me(r, r, emat, bg).first;
        mrs[me] += self_weight / 2;
        mrp[r] += self_weight / 2;
        mrm[r] += self_weight / 2;
        assert(mrs[me] >= 0);
    }

    for (auto e : npolicy.get_in_edges(v, g))
    {
        vertex_t u = source(e, g);
        if (u == v)
            continue;

        vertex_t s = b[u];

        typename graph_traits<BGraph>::edge_descriptor me;
        pair<typename graph_traits<BGraph>::edge_descriptor, bool> mep =
                get_me(s, r, emat, bg);


        if (!mep.second)
        {
            mep = add_edge(s, r, bg);
            put_me(s, r, mep.first, emat, bg);
            mrs[mep.first] = 0;
        }
        me = mep.first;

        size_t ew = eweight[e];

        mrs[me] += ew;

        mrp[s] += ew;
        mrm[r] += ew;
    }

    if (!overlap_stats.is_enabled())
    {
        wr[r] += vweight[v];
    }
    else
    {
        overlap_stats.add_half_edge(v, r, b, g);
        wr[r] = overlap_stats.get_block_size(r);
    }

    b[v] = r;
}


// move a vertex from its current block to block nr
template <class Graph, class BGraph, class Eprop, class Vprop, class EWprop,
          class VWprop, class EMat, class OStats, class PStats, class Vec,
          class NPolicy = standard_neighbours_policy>
void move_vertex(const Vec& vs, size_t nr, Eprop& mrs, Vprop& mrp, Vprop& mrm,
                 Vprop& wr, Vprop& b, bool deg_corr, const EWprop& eweight,
                 const VWprop& vweight, Graph& g, BGraph& bg, EMat& emat,
                 OStats& overlap_stats,  PStats& partition_stats,
                 const NPolicy& npolicy = NPolicy())
{
    if (b[vs[0]] == int(nr))
        return;

    size_t kin = 0, kout = 0;
    for (auto v : vs)
    {
        kin += in_degreeS()(v, g, eweight);
        kout += out_degreeS()(v, g, eweight);
    }

    if (partition_stats.is_enabled())
        partition_stats.move_vertex(vs[0], b[vs[0]], nr, deg_corr, overlap_stats,
                                    g, eweight, kin, kout);

    for (auto v : vs)
    {
        remove_vertex(v, mrs, mrp, mrm, wr, b, eweight, vweight, g, bg, emat,
                      overlap_stats, npolicy);
        add_vertex(v, nr, mrs, mrp, mrm, wr, b, eweight, vweight, g, bg, emat,
                   overlap_stats, npolicy);
    }
}

template <class Graph, class BGraph, class Eprop, class Vprop, class EWprop,
          class VWprop, class EMat, class OStats, class PStats,
          class NPolicy = standard_neighbours_policy>
void move_vertex(size_t v, size_t nr, Eprop& mrs, Vprop& mrp, Vprop& mrm,
                 Vprop& wr, Vprop& b, bool deg_corr, const EWprop& eweight,
                 const VWprop& vweight, Graph& g, BGraph& bg, EMat& emat,
                 OStats& overlap_stats, PStats& partition_stats,
                 const NPolicy& npolicy = NPolicy())
{
    std::array<size_t, 1> vs = {v};
    move_vertex(vs, nr, mrs, mrp, mrm, wr, b, deg_corr, eweight, vweight, g, bg,
                emat, overlap_stats, partition_stats, npolicy);
}


template <class Type1, class Type2, class Graph>
__attribute__((always_inline))
inline pair<Type1,Type2> make_ordered_pair(const Type1& v1, const Type2& v2, const Graph&)
{
    if (!is_directed::apply<Graph>::type::value)
    {
        if (v1 < v2)
            return make_pair(v1, v2);
        else
            return make_pair(v2, v1);
    }
    return make_pair(v1, v2);
}

template <class Graph>
class EntrySet
{
public:
    EntrySet() {}
    EntrySet(size_t B)
    {
        _null = numeric_limits<size_t>::max();
        _r_field_t.resize(B, _null);
        _nr_field_t.resize(B, _null);

        if (is_directed::apply<Graph>::type::value)
        {
            _r_field_s.resize(B, _null);
            _nr_field_s.resize(B, _null);
        }
        _entries.reserve(B);
        _delta.reserve(B);
    }

    __attribute__((always_inline))
    void set_move(size_t r, size_t nr)
    {
        _rnr = make_pair(r, nr);
    }

    __attribute__((always_inline))
    void insert_delta(size_t r, size_t s, int delta, bool source)
    {
        if (s == _rnr.first || s == _rnr.second)
        {
            if ((!is_directed::apply<Graph>::type::value && s < r) || source)
                std::swap(r, s);
            if (source)
                source = false;
        }

        if (source && (s == r))
            source = false;

        auto& r_field = (source) ? _r_field_s : _r_field_t;
        auto& nr_field = (source) ? _nr_field_s : _nr_field_t;

        vector<size_t>& field = (_rnr.first == r) ? r_field : nr_field;
        if (field[s] == _null)
        {
            field[s] = _entries.size();
            if ((!is_directed::apply<Graph>::type::value && s < r) || source)
                _entries.emplace_back(s, r);
            else
                _entries.emplace_back(r, s);
            _delta.push_back(delta);
        }
        else
        {
            _delta[field[s]] += delta;
        }
    }

    __attribute__((always_inline))
    int get_delta(size_t t, size_t s)
    {
        if (is_directed::apply<Graph>::type::value)
        {
            if (t == _rnr.first || t == _rnr.second)
                return get_delta_target(t, s);
            if (s == _rnr.first || s == _rnr.second)
                return get_delta_source(t, s);
        }
        else
        {
            if (t == _rnr.first || t == _rnr.second)
                return get_delta_target(t, s);
            if (s == _rnr.first || s == _rnr.second)
                return get_delta_target(s, t);
        }
        return 0;
    }

    __attribute__((always_inline))
    int get_delta_target(size_t r, size_t s)
    {
        vector<size_t>& field = (_rnr.first == r) ? _r_field_t : _nr_field_t;
        if (field[s] == _null)
            return 0;
        else
            return _delta[field[s]];
    }

    __attribute__((always_inline))
    int get_delta_source(size_t s, size_t r)
    {
        vector<size_t>& field = (_rnr.first == r) ? _r_field_s : _nr_field_s;
        if (field[s] == _null)
            return 0;
        else
            return _delta[field[s]];
    }

    __attribute__((always_inline))
    void clear()
    {
        for (size_t i = 0; i < _entries.size(); ++i)
        {
            size_t r = _entries[i].first;
            size_t s = _entries[i].second;
            _r_field_t[r] = _nr_field_t[r] = _null;
            _r_field_t[s] = _nr_field_t[s] = _null;
            if (is_directed::apply<Graph>::type::value)
            {
                _r_field_s[r] = _nr_field_s[r] = _null;
                _r_field_s[s] = _nr_field_s[s] = _null;
            }
        }
        _entries.clear();
        _delta.clear();
    }

    vector<pair<size_t, size_t> >& get_entries() { return _entries; }
    vector<int>& get_delta() { return _delta; }

private:
    pair<size_t, size_t> _rnr;
    size_t _null;
    vector<size_t> _r_field_t;
    vector<size_t> _nr_field_t;
    vector<size_t> _r_field_s;
    vector<size_t> _nr_field_s;
    vector<pair<size_t, size_t> > _entries;
    vector<int> _delta;
};

struct is_loop_nop
{
    bool operator()(size_t) const { return false; }
};


// obtain the necessary entries in the e_rs matrix which need to be modified
// after the move
template <class Graph, class BGraph, class Vertex, class Vprop, class Eprop,
          class MEntries,  class NPolicy = standard_neighbours_policy,
          class IL = is_loop_nop>
void move_entries(Vertex v, Vertex nr, Vprop& b, Eprop& eweights, Graph& g,
                  BGraph&, MEntries& m_entries,
                  const NPolicy& npolicy = NPolicy(), IL is_loop = IL())
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    vertex_t r = b[v];

    m_entries.set_move(r, nr);

    int self_weight = 0;
    for (auto e : npolicy.get_out_edges(v, g))
    {
        vertex_t u = target(e, g);
        vertex_t s = b[u];
        int ew = eweights[e];
        //assert(ew > 0);

        m_entries.insert_delta(r, s, -ew, false);

        if (u == v || is_loop(v))
        {
            s = nr;
            if (!is_directed::apply<Graph>::type::value)
                self_weight += ew;
        }

        m_entries.insert_delta(nr, s, +ew, false);
    }

    if (self_weight > 0 && self_weight % 2 == 0)
    {
        m_entries.insert_delta(r,   r,  self_weight / 2, false);
        m_entries.insert_delta(nr, nr, -self_weight / 2, false);
    }

    for (auto e : npolicy.get_in_edges(v, g))
    {
        vertex_t u = source(e, g);
        if (u == v)
            continue;
        vertex_t s = b[u];
        int ew = eweights[e];

        m_entries.insert_delta(r,  s, -ew, true);
        m_entries.insert_delta(nr, s, +ew, true);
    }
}

// obtain the entropy difference given a set of entries in the e_rs matrix
template <class MEntries, class Eprop, class BGraph, class EMat>
double entries_dS(MEntries& m_entries, Eprop& mrs, BGraph& bg, EMat& emat)
{
    auto& entries = m_entries.get_entries();
    auto& delta = m_entries.get_delta();

    double dS = 0;
    for (size_t i = 0; i < entries.size(); ++i)
    {
        auto er = entries[i].first;
        auto es = entries[i].second;
        int d = delta[i];

        int ers = get_mrs(er, es, mrs, emat, bg);
        assert(ers + d >= 0);
        dS += eterm(er, es, ers + d, bg) - eterm(er, es, ers, bg);
    }
    return dS;
}

// compute the entropy difference of a virtual move of vertex from block r to nr
template <class Graph, class BGraph, class Eprop, class Vprop, class EWprop,
          class VWprop, class EMat, class MEntries, class OStats, class Vec,
          class NPolicy = standard_neighbours_policy, class IL = is_loop_nop>
double virtual_move_sparse(const Vec& vs, size_t nr, Eprop& mrs, Vprop& mrp,
                           Vprop& mrm, Vprop& wr, Vprop& b, bool deg_corr,
                           const EWprop& eweight, const VWprop& vweight,
                           Graph& g, BGraph& bg, EMat& emat,
                           MEntries& m_entries, OStats& overlap_stats,
                           bool parallel_edges,
                           const NPolicy& npolicy = NPolicy(),
                           IL is_loop = IL())

{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    vertex_t r = b[vs[0]];

    if (r == nr)
        return 0.;

    m_entries.clear();
    int kin = 0, kout = 0;
    int self_weight = 0;
    for (auto v : vs)
    {
        move_entries(v, nr, b, eweight, g, bg, m_entries, npolicy, is_loop);
        kout += npolicy.get_out_degree(v, g, eweight);
        if (is_directed::apply<Graph>::type::value)
        {
            kin += npolicy.get_in_degree(v, g, eweight);
        }
        else
        {
            if (is_loop(v))
                self_weight++;
        }
    }

    if (self_weight > 0 && self_weight % 2 == 0)
    {
        m_entries.insert_delta(r,   r,  self_weight / 2, false);
        m_entries.insert_delta(nr, nr, -self_weight / 2, false);
    }

    double dS = entries_dS(m_entries, mrs, bg, emat);

    int dwr = 0, dwnr = 0;
    if (!overlap_stats.is_enabled())
    {
        for (auto v : vs)
            dwr += vweight[v];
        dwnr = dwr;
    }
    else
    {
        dwr = wr[r] - overlap_stats.virtual_remove_size(vs[0], r, kin, kout);
        dwnr = overlap_stats.virtual_add_size(vs[0], nr) - wr[nr];
        if (deg_corr)
            dS += overlap_stats.virtual_move_dS(vs[0], r, nr, g, kin, kout);
        if (parallel_edges)
            dS += overlap_stats.virtual_move_parallel_dS(vs[0], r, nr, b, g);
    }

    if (!is_directed::apply<Graph>::type::value)
        kin = kout;

    dS += vterm(mrp[r]  - kout, mrm[r]  - kin, wr[r]  - dwr , deg_corr, bg);
    dS += vterm(mrp[nr] + kout, mrm[nr] + kin, wr[nr] + dwnr, deg_corr, bg);
    dS -= vterm(mrp[r]        , mrm[r]       , wr[r]        , deg_corr, bg);
    dS -= vterm(mrp[nr]       , mrm[nr]      , wr[nr]       , deg_corr, bg);

    return dS;
}

template <class Graph, class BGraph, class Eprop, class Vprop, class EWprop,
          class VWprop, class EMat, class MEntries, class OStats,
          class NPolicy = standard_neighbours_policy, class IL = is_loop_nop>
double virtual_move_sparse(size_t v, size_t nr, Eprop& mrs, Vprop& mrp,
                           Vprop& mrm, Vprop& wr, Vprop& b, bool deg_corr,
                           const EWprop& eweight, const VWprop& vweight,
                           Graph& g, BGraph& bg, EMat& emat,
                           MEntries& m_entries, OStats& overlap_stats,
                           bool parallel_edges,
                           const NPolicy& npolicy = NPolicy(),
                           IL is_loop = IL())

{
    std::array<size_t, 1> vs = {v};
    return virtual_move_sparse(vs, nr, mrs, mrp, mrm, wr, b, deg_corr, eweight,
                               vweight, g, bg, emat, m_entries, overlap_stats,
                               parallel_edges, npolicy, is_loop);
}


// compute the entropy difference of a virtual move of vertex r to block nr
template <class Graph, class BGraph, class Eprop, class Vprop, class EWprop,
          class VWprop, class EMat, class MEntries, class OStats, class Vec,
          class NPolicy = standard_neighbours_policy, class IL = is_loop_nop>
double virtual_move_dense(const Vec& vs, size_t nr, Eprop& mrs, Vprop&,
                          Vprop&, Vprop& wr, Vprop& b, bool deg_corr,
                          const EWprop& eweight, const VWprop& vweight,
                          Graph& g, BGraph& bg, EMat& emat, MEntries& m_entries,
                          OStats& overlap_stats, bool multigraph,
                          const NPolicy& npolicy = NPolicy(), IL is_loop = IL())
{
    if (deg_corr)
        throw GraphException("Dense entropy for degree corrected model not implemented!");

    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    vertex_t r = b[vs[0]];

    if (r == nr)
        return 0;

    // m_entries is not used in the computation below, but it is expected afterwards
    m_entries.clear();
    int kin = 0, kout = 0;
    int self_weight = 0;
    for (auto v : vs)
    {
        move_entries(v, nr, b, eweight, g, bg, m_entries, npolicy, is_loop);
        kout += npolicy.get_out_degree(v, g, eweight);
        if (is_directed::apply<Graph>::type::value)
        {
            kin += npolicy.get_in_degree(v, g, eweight);
        }
        else
        {
            if (is_loop(v))
                self_weight++;
        }
    }

    if (self_weight > 0 && self_weight % 2 == 0)
    {
        m_entries.insert_delta(r,   r,  self_weight / 2, false);
        m_entries.insert_delta(nr, nr, -self_weight / 2, false);
    }

    vector<int> deltap(num_vertices(bg), 0);
    int deltal = 0;
    for (auto v : vs)
    {
        for (auto e : npolicy.get_out_edges(v, g))
        {
            vertex_t u = target(e, g);
            vertex_t s = b[u];
            if (u == v)
            {
                deltal += eweight[e];
            }
            else
            {
                deltap[s] += eweight[e];
            }
        }
    }
    if (!is_directed::apply<Graph>::type::value)
        deltal /= 2;

    vector<int> deltam(num_vertices(bg), 0);
    for (auto v : vs)
    {
        for (auto e : npolicy.get_in_edges(v, g))
        {
            vertex_t u = source(e, g);
            if (u == v)
                continue;
            vertex_t s = b[u];
            deltam[s] += eweight[e];
        }
    }

    double dS = 0;
    int dwr = 0, dwnr = 0;
    if (!overlap_stats.is_enabled())
    {
        for (auto v : vs)
            dwr += vweight[v];
        dwnr = dwr;
    }
    else
    {
        dwr = wr[r] - overlap_stats.virtual_remove_size(vs[0], r, kin, kout);
        dwnr = overlap_stats.virtual_add_size(vs[0], nr) - wr[nr];
        if (deg_corr)
            dS += overlap_stats.virtual_move_dS(vs[0], r, nr, g, kin, kout);
        if (multigraph)
            dS += overlap_stats.virtual_move_parallel_dS(vs[0], r, nr, b, g);
    }

    double Si = 0, Sf = 0;
    for (vertex_t s = 0; s < num_vertices(bg); ++s)
    {
        int ers = get_mrs(r, s, mrs, emat, bg);
        int enrs = get_mrs(nr, s, mrs, emat, bg);

        if (!is_directed::apply<Graph>::type::value)
        {
            if (s != nr && s != r)
            {
                Si += eterm_dense(r,  s, ers,              wr[r],         wr[s], multigraph, bg);
                Sf += eterm_dense(r,  s, ers - deltap[s],  wr[r] - dwr,   wr[s], multigraph, bg);
                Si += eterm_dense(nr, s, enrs,             wr[nr],        wr[s], multigraph, bg);
                Sf += eterm_dense(nr, s, enrs + deltap[s], wr[nr] + dwnr, wr[s], multigraph, bg);
            }

            if (s == r)
            {
                Si += eterm_dense(r, r, ers,                      wr[r],       wr[r],       multigraph, bg);
                Sf += eterm_dense(r, r, ers - deltap[r] - deltal, wr[r] - dwr, wr[r] - dwr, multigraph, bg);
            }

            if (s == nr)
            {
                Si += eterm_dense(nr, nr, enrs,                       wr[nr],        wr[nr],        multigraph, bg);
                Sf += eterm_dense(nr, nr, enrs + deltap[nr] + deltal, wr[nr] + dwnr, wr[nr] + dwnr, multigraph, bg);

                Si += eterm_dense(r, nr, ers,                          wr[r],       wr[nr],        multigraph, bg);
                Sf += eterm_dense(r, nr, ers - deltap[nr] + deltap[r], wr[r] - dwr, wr[nr] + dwnr, multigraph, bg);
            }
        }
        else
        {
            int esr = get_mrs(s, r, mrs, emat, bg);
            int esnr = get_mrs(s, nr, mrs, emat, bg);

            if (s != nr && s != r)
            {
                Si += eterm_dense(r, s, ers            , wr[r]      , wr[s]      , multigraph, bg);
                Sf += eterm_dense(r, s, ers - deltap[s], wr[r] - dwr, wr[s]      , multigraph, bg);
                Si += eterm_dense(s, r, esr            , wr[s]      , wr[r]      , multigraph, bg);
                Sf += eterm_dense(s, r, esr - deltam[s], wr[s]      , wr[r] - dwr, multigraph, bg);

                Si += eterm_dense(nr, s, enrs            , wr[nr]       , wr[s]        , multigraph, bg);
                Sf += eterm_dense(nr, s, enrs + deltap[s], wr[nr] + dwnr, wr[s]        , multigraph, bg);
                Si += eterm_dense(s, nr, esnr            , wr[s]        , wr[nr]       , multigraph, bg);
                Sf += eterm_dense(s, nr, esnr + deltam[s], wr[s]        , wr[nr] + dwnr, multigraph, bg);
            }

            if(s == r)
            {
                Si += eterm_dense(r, r, ers                                  , wr[r]      , wr[r]      , multigraph, bg);
                Sf += eterm_dense(r, r, ers - deltap[r]  - deltam[r] - deltal, wr[r] - dwr, wr[r] - dwr, multigraph, bg);

                Si += eterm_dense(r, nr, esnr                         , wr[r]      , wr[nr]       , multigraph, bg);
                Sf += eterm_dense(r, nr, esnr - deltap[nr] + deltam[r], wr[r] - dwr, wr[nr] + dwnr, multigraph, bg);
            }

            if(s == nr)
            {
                Si += eterm_dense(nr, nr, esnr                                   , wr[nr]       , wr[nr]       , multigraph, bg);
                Sf += eterm_dense(nr, nr, esnr + deltap[nr] + deltam[nr] + deltal, wr[nr] + dwnr, wr[nr] + dwnr, multigraph, bg);

                Si += eterm_dense(nr, r, esr                         , wr[nr]       , wr[r]      , multigraph, bg);
                Sf += eterm_dense(nr, r, esr + deltap[r] - deltam[nr], wr[nr] + dwnr, wr[r] - dwr, multigraph, bg);
            }
        }
    }

    return Sf - Si + dS;
}

// compute the entropy difference of a virtual move of vertex r to block nr
template <class Graph, class BGraph, class Eprop, class Vprop, class EWprop,
          class VWprop, class EMat, class MEntries, class OStats,
          class NPolicy = standard_neighbours_policy, class IL = is_loop_nop>
double virtual_move_dense(size_t v, size_t nr, Eprop& mrs, Vprop& mrp,
                          Vprop& mrm, Vprop& wr, Vprop& b, bool deg_corr,
                          const EWprop& eweight, const VWprop& vweight,
                          Graph& g, BGraph& bg, EMat& emat, MEntries& m_entries,
                          OStats& overlap_stats, bool multigraph,
                          const NPolicy& npolicy = NPolicy(), IL is_loop = IL())
{
    std::array<size_t, 1> vs = {v};
    return virtual_move_dense(vs, nr, mrs, mrp, mrm, wr, b, deg_corr, eweight,
                              vweight, g, bg, emat, m_entries, overlap_stats,
                              multigraph, npolicy, is_loop);
}

template <class Graph, class BGraph, class Eprop, class Vprop, class EWprop,
          class VWprop, class EMat, class MEntries, class OStats, class Vec,
          class NPolicy = standard_neighbours_policy, class IL = is_loop_nop>
double virtual_move(const Vec& vs, size_t nr, bool dense, Eprop& mrs,
                    Vprop& mrp, Vprop& mrm, Vprop& wr, Vprop& b, bool deg_corr,
                    const EWprop& eweight, const VWprop& vweight, Graph& g,
                    BGraph& bg, EMat& emat, MEntries& m_entries,
                    OStats& overlap_stats, bool parallel_edges,
                    const NPolicy& npolicy = NPolicy(), IL is_loop = IL())
{
    double S = 0;
    if (dense)
    {
        S = virtual_move_dense(vs, nr, mrs, mrp, mrm, wr, b, deg_corr, eweight,
                               vweight, g, bg, emat, m_entries, overlap_stats,
                               parallel_edges, npolicy, is_loop);
    }
    else
    {
        S = virtual_move_sparse(vs, nr, mrs, mrp, mrm, wr, b, deg_corr, eweight,
                                vweight, g, bg, emat, m_entries, overlap_stats,
                                parallel_edges, npolicy, is_loop);
    }

    return S;
}

// ====================================
// Construct and manage half-edge lists
// ====================================

//the following guarantees a stable (source, target) ordering even for
//undirected graphs
template <class Edge, class Graph>
inline typename graph_traits<Graph>::vertex_descriptor
get_source(const Edge& e, const Graph &g)
{
    if (is_directed::apply<Graph>::type::value)
        return source(e, g);
    return min(source(e, g), target(e, g));
}

template <class Edge, class Graph>
inline typename graph_traits<Graph>::vertex_descriptor
get_target(const Edge& e, const Graph &g)
{
    if (is_directed::apply<Graph>::type::value)
        return target(e, g);
    return max(source(e, g), target(e, g));
}

struct egroups_manage
{

    template <class Graph, class Weighted>
    struct get_sampler
    {
        typedef typename mpl::if_<Weighted,
                                  DynamicSampler<std::tuple<typename graph_traits<Graph>::edge_descriptor, bool>>,
                                  vector<std::tuple<typename graph_traits<Graph>::edge_descriptor, bool>>>::type type;
    };

    template <class Eprop, class Vprop, class VEprop, class Graph, class VertexIndex>
    static void build(Vprop b, boost::any& oegroups, VEprop esrcpos,
                      VEprop etgtpos, Eprop eweight, Graph& g,
                      VertexIndex vertex_index, size_t B, bool weighted,
                      bool empty)
    {
        if (weighted)
        {
            typedef typename get_sampler<Graph, mpl::true_>::type sampler_t;
            typedef typename property_map_type::apply<sampler_t,
                                                      VertexIndex>::type vemap_t;
            vemap_t egroups_checked(vertex_index);
            oegroups = egroups_checked;
            if (empty)
                return;
            build_dispatch(b, egroups_checked.get_unchecked(B),
                           esrcpos, etgtpos, eweight, g, vertex_index, B);
        }
        else
        {
            typedef typename get_sampler<Graph, mpl::false_>::type sampler_t;
            typedef typename property_map_type::apply<sampler_t,
                                                      VertexIndex>::type vemap_t;
            vemap_t egroups_checked(vertex_index);
            oegroups = egroups_checked;
            if (empty)
                return;
            build_dispatch(b, egroups_checked.get_unchecked(B),
                           esrcpos, etgtpos, eweight, g, vertex_index, B);
        }
    }

    template <class Eprop, class Vprop, class VEprop, class Graph,
              class VertexIndex, class Egroups>
    static void build_dispatch(Vprop b, Egroups egroups, VEprop esrcpos,
                               VEprop etgtpos, Eprop eweight, Graph& g,
                               VertexIndex, size_t)
    {
        for (auto e : edges_range(g))
        {
            size_t r = b[get_source(e, g)];
            assert (r < B);
            auto& r_elist = egroups[r];
            esrcpos[e] = insert_edge(std::make_tuple(e, true), r_elist,
                                     eweight[e]);

            size_t s = b[get_target(e, g)];
            assert (s < B);
            auto& s_elist = egroups[s];
            etgtpos[e] = insert_edge(std::make_tuple(e, false), s_elist,
                                     eweight[e]);
        }
    }

    template <class Edge, class EV>
    static size_t insert_edge(const Edge& e, EV& elist, size_t)
    {
        elist.push_back(e);
        return elist.size() - 1;
    }

    template <class Edge>
    static size_t insert_edge(const Edge& e, DynamicSampler<Edge>& elist,
                              size_t weight)
    {
        return elist.insert(e, weight);
    }


    template <class Edge, class Epos>
    static void remove_edge(size_t pos, Epos& esrcpos, Epos& etgtpos,
                            vector<Edge>& elist)
    {
        typedef typename property_traits<Epos>::value_type val_t;
        if (get<1>(elist.back()))
            esrcpos[get<0>(elist.back())] = pos;
        else
            etgtpos[get<0>(elist.back())] = pos;
        if (get<1>(elist[pos]))
            esrcpos[get<0>(elist[pos])] = numeric_limits<val_t>::max();
        else
            etgtpos[get<0>(elist[pos])] = numeric_limits<val_t>::max();
        elist[pos] = elist.back();
        elist.pop_back();
    }

    template <class Edge, class Epos>
    static void remove_edge(size_t pos, Epos& esrcpos, Epos& etgtpos,
                            DynamicSampler<Edge>& elist)
    {
        typedef typename property_traits<Epos>::value_type val_t;
        if (get<1>(elist[pos]))
            esrcpos[get<0>(elist[pos])] = numeric_limits<val_t>::max();
        else
            etgtpos[get<0>(elist[pos])] = numeric_limits<val_t>::max();
        elist.remove(pos);
    }

    template <class Vertex, class Graph, class EVprop, class Eprop, class VEprop>
    static void remove_egroups(Vertex v, Vertex r, Eprop&,
                               EVprop& egroups, VEprop& esrcpos, VEprop& etgtpos,
                               Graph& g)
    {
        typedef Vertex vertex_t;

        auto& elist = egroups[r];

        //update the half-edge lists
        for (auto e : all_edges_range(v, g))
        {
            vertex_t src = get_source(e, g);
            vertex_t tgt = get_target(e, g);

            bool is_src = (src == v);

            // self-loops will appear twice; we must disambiguate
            if (src == tgt)
            {
                size_t pos = esrcpos[e];
                is_src = (pos < elist.size() && get<0>(elist[pos]) == e);
            }

            size_t pos = (is_src) ? esrcpos[e] : etgtpos[e];
            assert(pos < elist.size());
            assert(get<0>(elist[pos]) == e);
            remove_edge(pos, esrcpos, etgtpos, elist);
        }
    }

    template <class Vertex, class Graph, class EVprop, class Eprop, class VEprop>
    static void add_egroups(Vertex v, Vertex s, Eprop& eweight, EVprop& egroups,
                            VEprop& esrcpos, VEprop& etgtpos, Graph& g)
    {
        typedef Vertex vertex_t;

        auto& elist = egroups[s];

        //update the half-edge lists
        for (auto e : all_edges_range(v, g))
        {
            vertex_t src = get_source(e, g);
            vertex_t tgt = get_target(e, g);

            bool is_src = (src == v);

            // self-loops will appear twice; we must disambiguate
            if (src == tgt)
            {
                size_t pos = esrcpos[e];
                is_src = !(pos < elist.size() && get<0>(elist[pos]) == e);
            }

            typedef typename tuple_element<0, typename property_traits<EVprop>::value_type::value_type>::type e_type;
            size_t pos = insert_edge(std::make_tuple(e_type(e), is_src),
                                     elist, size_t(eweight[e]));
            if (is_src)
                esrcpos[e] = pos;
            else
                etgtpos[e] = pos;
        }
    }

    template <class Vertex, class Graph, class EVprop, class Eprop, class VEprop>
    static void update_egroups(Vertex v, Vertex r, Vertex s, Eprop& eweight,
                               EVprop& egroups, VEprop& esrcpos, VEprop& etgtpos,
                               Graph& g)
    {
        if (r == s)
            return;
        remove_egroups(v, r, eweight, egroups, esrcpos, etgtpos, g);
        add_egroups(v, s, eweight, egroups, esrcpos, etgtpos, g);
    }

    template <class Edge, class RNG>
    static typename std::tuple_element<0, Edge>::type
    sample_edge(DynamicSampler<Edge>& elist, RNG& rng)
    {
        return get<0>(elist.sample(rng));
    }

    template <class Edge, class RNG>
    static typename std::tuple_element<0, Edge>::type
    sample_edge(vector<Edge>& elist, RNG& rng)
    {
        std::uniform_int_distribution<size_t> urand(0, elist.size() - 1);
        size_t ur = urand(rng);
        return get<0>(elist[ur]);
    }
};

//============
// Main loops
//============

//computes the move proposal probability
template <class Vertex, class Graph, class Vprop, class Eprop, class MEprop,
          class Emat, class BGraph, class MEntries, class OStats>
inline double
get_move_prob(Vertex v, Vertex r, Vertex s, double c, Vprop& b, MEprop& mrs,
              Vprop& mrp, Vprop& mrm, Emat& emat, Eprop& eweight, Graph& g,
              BGraph& bg, MEntries& m_entries, bool reverse,
              OStats& overlap_stats, size_t kin = 0, size_t kout = 0)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    size_t B = num_vertices(bg);
    double p = 0;
    size_t w = 0;

    if (reverse && (kin + kout) == 0)
    {
        kout = out_degreeS()(v, g, eweight);
        kin = kout;
        if (is_directed::apply<Graph>::type::value)
            kin = in_degreeS()(v, g, eweight);
    }

    vector<Vertex> vertices;  // FIXME: Suboptimal performance
    if (!overlap_stats.is_enabled())
    {
        vertices = {v};
    }
    else
    {
        size_t vi = overlap_stats.get_node(v);
        auto& ns = overlap_stats.get_half_edges(vi);
        vertices.insert(vertices.end(), ns.begin(), ns.end());
    }

    for (auto vv : vertices)
    {
        for (auto e : all_edges_range(vv, g))
        {
            vertex_t u = target(e, g);
            if (is_directed::apply<Graph>::type::value && u == vv)
                u = source(e, g);
            vertex_t t = b[u];
            if (u == v)
                t = r;
            size_t ew = eweight[e];
            w += ew;

            int mts = get_mrs(t, s, mrs, emat, bg);
            int mtp = mrp[t];
            int mst = mts;
            int mtm = mtp;

            if (is_directed::apply<Graph>::type::value)
            {
                mst = get_mrs(s, t, mrs, emat, bg);
                mtm = mrm[t];
            }

            if (reverse)
            {
                int dts = m_entries.get_delta(t, s);
                int dst = dts;
                if (is_directed::apply<Graph>::type::value)
                    dst = m_entries.get_delta(s, t);

                mts += dts;
                mst += dst;

                if (t == s)
                {
                    mtp -= kout;
                    mtm -= kin;
                }

                if (t == r)
                {
                    mtp += kout;
                    mtm += kin;
                }
            }

            if (is_directed::apply<Graph>::type::value)
            {
                p += ew * ((mts + mst + c) / (mtp + mtm + c * B));
            }
            else
            {
                if (t == s)
                    mts *= 2;
                p += ew * (mts + c) / (mtp + c * B);
            }
        }
    }
    return p / w;
}


// merge vertex u into v
template <class Graph, class Eprop, class Vprop>
void merge_vertices(size_t u, size_t v, Eprop& eweight_u, Vprop& vweight, Graph& g)
{
    auto eweight = eweight_u.get_checked();

    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    typedef typename graph_traits<Graph>::edge_descriptor edge_t;

    unordered_map<vertex_t, vector<edge_t>> ns_u, ns_v;
    for(auto e : out_edges_range(u, g))
        ns_u[target(e, g)].push_back(e);
    for(auto e : out_edges_range(v, g))
        ns_v[target(e, g)].push_back(e);

    for(auto& kv : ns_u)
    {
        vertex_t t = kv.first;
        auto& es = kv.second;

        size_t w = 0;
        for (auto& e : es)
            w += eweight_u[e];

        if (t == u)
        {
            t = v;
            if (!is_directed::apply<Graph>::type::value)
            {
                assert(w % 2 == 0);
                w /= 2;
            }
        }

        auto iter = ns_v.find(t);
        if (iter != ns_v.end())
        {
            auto& e = iter->second.front();
            eweight_u[e] += w;
        }
        else
        {
            auto e = add_edge(v, t, g).first;
            ns_v[t].push_back(e);
            eweight[e] = w;
        }
    }

    if (is_directed::apply<Graph>::type::value)
    {
        ns_u.clear();
        ns_v.clear();

        for(auto e : in_edges_range(v, g))
            ns_v[source(e, g)].push_back(e);
        for(auto e : in_edges_range(u, g))
            ns_u[source(e, g)].push_back(e);

        for(auto& kv : ns_u)
        {
            size_t s = kv.first;
            auto& es = kv.second;

            if (s == u)
                continue;

            size_t w = 0;
            for (auto& e : es)
                w += eweight_u[e];

            auto iter = ns_v.find(s);
            if (iter != ns_v.end())
            {
                auto& e = iter->second.front();
                eweight_u[e] += w;
            }
            else
            {
                auto e = add_edge(s, v, g).first;
                ns_v[s].push_back(e);
                eweight[e] = w;
            }
        }
    }

    vweight[v] += vweight[u];
    vweight[u] = 0;
    clear_vertex(u, g);
}

template <class Graph, class BGraph, class EMprop, class Eprop, class Vprop,
          class EMat, class EVprop, class VEprop, class SamplerMap,
          class PartitionStats, class OStats, class BMap, class BRMap>
class BlockState
{
public:
    Graph& g;
    Eprop eweight;
    Vprop vweight;
    Vprop b;
    BGraph& bg;
    EMat& emat;
    EMprop mrs;
    Vprop mrp;
    Vprop mrm;
    Vprop wr;
    EVprop egroups;
    VEprop esrcpos;
    VEprop etgtpos;
    SamplerMap& neighbour_sampler;
    SamplerMap& cavity_neighbour_sampler;
    PartitionStats& partition_stats;
    OStats& overlap_stats;
    BMap& block_map;
    BRMap block_rmap;
    vector<size_t>& free_blocks;
    bool master;
    bool slave;
    bool single;
};

template <class Graph, class BGraph, class EMprop, class Eprop, class Vprop,
          class EMat, class EVprop, class VEprop, class SamplerMap,
          class PartitionStats, class OStats, class BMap, class BRMap>
BlockState<Graph, BGraph, EMprop, Eprop, Vprop, EMat, EVprop, VEprop,
           SamplerMap, PartitionStats, OStats, BMap, BRMap>
make_block_state(Graph& g, Eprop eweight, Vprop vweight, Vprop b, BGraph& bg,
                 EMat& emat, EMprop mrs, Vprop mrp, Vprop mrm, Vprop wr,
                 EVprop egroups, VEprop esrcpos, VEprop etgtpos,
                 SamplerMap& neighbour_sampler,
                 SamplerMap& cavity_neighbour_sampler,
                 PartitionStats& partition_stats, OStats& ostats,
                 BMap& block_map, BRMap block_rmap, vector<size_t>& free_blocks,
                 bool master, bool slave, bool single)
{
    BlockState<Graph, BGraph, EMprop, Eprop, Vprop, EMat, EVprop, VEprop,
               SamplerMap, PartitionStats, OStats, BMap, BRMap>
        state = {g, eweight, vweight, b, bg, emat, mrs, mrp, mrm, wr, egroups,
                 esrcpos, etgtpos, neighbour_sampler, cavity_neighbour_sampler,
                 partition_stats, ostats, block_map, block_rmap, free_blocks,
                 master, slave, single};
    return state;
}


template <class BGraph, class Eprop, class EMat, class ES>
double virtual_move_covariate(Eprop& mrs, BGraph& bg, EMat& emat,
                              ES& m_entries)
{
    typedef typename graph_traits<BGraph>::vertex_descriptor vertex_t;

    auto& entries = m_entries.get_entries();
    auto& delta = m_entries.get_delta();

    double dS = 0;
    for (size_t i = 0; i < entries.size(); ++i)
    {
        vertex_t er = entries[i].first;
        vertex_t es = entries[i].second;
        int d = delta[i];

        int ers = get_mrs(er, es, mrs, emat, bg);
        assert(ers + d >= 0);
        dS -= -lgamma_fast(ers + 1);
        dS += -lgamma_fast(ers + d + 1);
    }
    return dS;
}

template <class VLprop, class VVprop, class Vprop, class BlockState,
          class MEntries, class NPolicy = standard_neighbours_policy>
double virtual_move(size_t v, size_t s, Vprop& b, VLprop& cv, VVprop& vmap,
                    vector<BlockState>& states, vector<MEntries>& m_entries,
                    bool dense, bool deg_corr, double multigraph,
                    const NPolicy& npolicy = NPolicy())
{
    if (s == size_t(b[v]))
        return 0;

    double dS = 0;
    auto& ls = cv[v];
    auto& vs = vmap[v];
    for (size_t j = 0; j < ls.size(); ++j)
    {
        size_t l = ls[j] + 1;
        size_t u = vs[j];

        auto& state = states[l];

        size_t r_u = state.b[u];

        assert(r_u < num_vertices(state.bg));

        size_t s_u = (l == 0) ? s : get_block_map(state, state.block_map, s);

        if ((!state.slave && l > 0) || state.master || state.single)
        {
                dS += virtual_move(u, s_u, dense, state.mrs, state.mrp, state.mrm,
                                   state.wr, state.b, deg_corr, state.eweight,
                                   state.vweight, state.g, state.bg, state.emat,
                                   m_entries[l], state.overlap_stats, multigraph,
                                   npolicy);
        }

        if (state.partition_stats.is_enabled())
        {
            if (l == 0)
                dS += state.partition_stats.get_delta_dl(u, r_u, s_u,
                                                         state.overlap_stats,
                                                         state.g);
            if (deg_corr && ((!state.slave && l > 0) || state.master || state.single))
                dS += state.partition_stats.get_delta_deg_dl(u, r_u, s_u,
                                                             state.eweight,
                                                             state.overlap_stats,
                                                             state.g);
        }

        if (state.master || state.slave)
        {
            if (state.slave)
            {
                m_entries[l].clear();
                move_entries(u, s_u, state.b, state.eweight, state.g,
                             state.bg, m_entries[l], npolicy);
            }

            double ddS = virtual_move_covariate(state.mrs, state.bg, state.emat,
                                                m_entries[l]);
            if (state.master)
                dS -= ddS;
            else
                dS += ddS;
        }

        if (l > 0 && state.wr[s_u] == 0)
            remove_block_map(state, state.block_map, s);
    }
    return dS;
}

// The following function will only be called for overlapping moves, where
// several half-edges are moved at once
template <class VLprop, class VVprop, class Vprop, class BlockState,
          class MEntries, class Vec, class NPolicy = standard_neighbours_policy>
double virtual_move(const Vec& vs, size_t s, Vprop& b, VLprop& cv, VVprop& vmap,
                    vector<BlockState>& states, vector<MEntries>& m_entries,
                    bool dense, bool deg_corr, double multigraph,
                    const NPolicy& npolicy = NPolicy())
{
    if (s == size_t(b[vs[0]]))
        return 0;

    unordered_map<int, vector<size_t>> us;
    for (size_t i = 0; i < vs.size(); ++i)
    {
        size_t v = vs[i];
        auto& ls = cv[v];
        auto& ws = vmap[v];
        for (size_t j = 0; j < ls.size(); ++j)
        {
            int l = ls[j];
            size_t u = ws[j];
            us[l].push_back(u);
        }
    }

    double dS = 0;
    for (auto& iter : us)
    {
        int l = iter.first + 1;
        auto& state = states[l];
        size_t u = iter.second[0];
        size_t r_u = state.b[u];
        size_t s_u = (l == 0) ? s : get_block_map(state, state.block_map, s);

        if ((!state.slave && l > 0) || state.master || state.single)
        {

            auto is_loop = [&] (size_t w) -> bool
                {
                    int t = state.overlap_stats.get_out_neighbour(w);
                    if (t == -1)
                    t = state.overlap_stats.get_in_neighbour(w);
                    return ((state.overlap_stats.get_node(w) ==
                             state.overlap_stats.get_node(t)) &&
                            (size_t(state.b[t]) == r_u) &&
                            (size_t(state.b[w]) == r_u));
                };

            dS += virtual_move(iter.second, s_u, dense, state.mrs, state.mrp,
                               state.mrm, state.wr, state.b, deg_corr,
                               state.eweight, state.vweight, state.g, state.bg,
                               state.emat, m_entries[l], state.overlap_stats,
                               multigraph && (dense || !state.overlap_stats.is_enabled()),
                               npolicy, is_loop);

            if (multigraph && !dense && state.overlap_stats.is_enabled())
            {
                std::unordered_map<int, size_t> bundles;
                auto& mi = state.overlap_stats.get_mi();
                for (auto w : iter.second)
                {
                    if (mi[w] != -1)
                        bundles[mi[w]] = w;
                }

                for (auto& kv : bundles)
                {
                    dS += state.overlap_stats.virtual_move_parallel_dS(kv.second,
                                                                       r_u, s_u,
                                                                       state.b,
                                                                       state.g,
                                                                       true);
                }
            }
        }

        size_t kin = 0, kout = 0;
        for (auto w : iter.second)
        {
            kin += in_degreeS()(w, state.g, state.eweight);
            kout += out_degreeS()(w, state.g, state.eweight);
        }

        if (state.partition_stats.is_enabled())
        {
            assert(r_u < num_vertices(state.bg));

            if (l == 0)
                dS += state.partition_stats.get_delta_dl(u, r_u, s_u,
                                                         state.overlap_stats,
                                                         state.g, kin, kout);
            if (deg_corr && ((!state.slave && l > 0) || state.master || state.single))
                dS += state.partition_stats.get_delta_deg_dl(u, r_u, s_u,
                                                             state.eweight,
                                                             state.overlap_stats,
                                                             state.g, kin, kout);
        }

        if (state.master || state.slave)
        {
            if (state.slave)
            {
                m_entries[l].clear();
                for (auto w : iter.second)
                    move_entries(w, s_u, state.b, state.eweight, state.g,
                                 state.bg, m_entries[l], npolicy);
            }

            double ddS = virtual_move_covariate(state.mrs, state.bg, state.emat,
                                                m_entries[l]);
            if (state.master)
                dS -= ddS;
            else
                dS += ddS;
        }

        if (l > 0 && state.wr[s_u] == 0)
            remove_block_map(state, state.block_map, s);
    }
    return dS;
}

template <class BlockState, class BlockMap>
size_t get_block_map(const BlockState&, const BlockMap&, size_t r)
{
    return r;
}

template <class BlockState, class BlockMap>
void remove_block_map(const BlockState&, const BlockMap&, size_t)
{
}


template <class BlockState, class VLprop, class VVprop, class Vprop,
          class NPolicy = standard_neighbours_policy>
void move_vertex(size_t v, size_t s, Vprop& b, VLprop& cv, VVprop& vmap,
                 bool deg_corr, vector<BlockState>& states, bool update_egroups,
                 const NPolicy& npolicy = NPolicy())
{
    size_t r = b[v];

    if (r == s)
        return;

    auto& ls = cv[v];
    auto& vs = vmap[v];
    for (size_t j = 0; j < ls.size(); ++j)
    {
        int l = ls[j] + 1;
        size_t u = vs[j];

        auto& state = states[l];
        size_t r_u = state.b[u];

        assert(r_u < num_vertices(state.bg));

        size_t s_u = (l == 0) ? s : get_block_map(state, state.block_map, s);

        move_vertex(u, s_u, state.mrs, state.mrp, state.mrm, state.wr, state.b,
                    deg_corr, state.eweight, state.vweight, state.g, state.bg,
                    state.emat, state.overlap_stats, state.partition_stats,
                    npolicy);

        if (l == 0 && update_egroups)
        {
            egroups_manage::update_egroups(u, r_u, s_u,
                                           state.eweight,
                                           state.egroups,
                                           state.esrcpos,
                                           state.etgtpos,
                                           state.g);
        }

        if (l > 0 && state.wr[r_u] == 0)
            remove_block_map(state, state.block_map, r);
    }
}


template <class BlockState, class VLprop, class VVprop, class Vprop,
          class Vec, class NPolicy = standard_neighbours_policy>
void move_vertex(const Vec& vs, size_t s, Vprop& b, VLprop& cv, VVprop& vmap,
                 bool deg_corr, vector<BlockState>& states, bool update_egroups,
                 const NPolicy& npolicy = NPolicy())
{
    size_t r = b[vs[0]];
    if (s == r)
        return;

    unordered_map<int, vector<size_t>> us;
    for (size_t i = 0; i < vs.size(); ++i)
    {
        auto& ls = cv[vs[i]];
        auto& ws = vmap[vs[i]];
        for (size_t j = 0; j < ls.size(); ++j)
        {
            int l = ls[j];
            size_t u = ws[j];
            us[l].push_back(u);
        }
    }

    for (auto iter : us)
    {
        int l = iter.first + 1;
        size_t u = iter.second[0];

        auto& state = states[l];
        size_t r_u = state.b[u];

        assert(r_u < num_vertices(state.bg));

        size_t s_u = (l == 0) ? s : get_block_map(state, state.block_map, s);

        move_vertex(iter.second, s_u, state.mrs, state.mrp, state.mrm, state.wr,
                    state.b, deg_corr, state.eweight, state.vweight, state.g,
                    state.bg, state.emat, state.overlap_stats,
                    state.partition_stats, npolicy);

        if (l == 0 && update_egroups)
        {
            for (auto w : iter.second)
                egroups_manage::update_egroups(w, r_u, s_u,
                                               state.eweight,
                                               state.egroups,
                                               state.esrcpos,
                                               state.etgtpos,
                                               state.g);
        }

        if (l > 0 && state.wr[r_u] == 0)
            remove_block_map(state, state.block_map, r);
    }
}


template <class Type>
void insert_vec(vector<Type>& v, size_t i, const Type& x)
{
    v.insert(v.begin() + i, x);
}

template <class Type>
void insert_vec(const std::array<Type, 1>&, size_t, const Type&)
{
}

template <class Type>
void clear_vec(vector<Type>& v)
{
    v.clear();
}

template <class Type>
void clear_vec(const std::array<Type, 1>&)
{
}


//A single Monte Carlo Markov chain sweep
template <class Graph, class Vprop, class VVprop, class VLprop,
          class Eprop, class RNG, class BlockState, class MEntries>
void move_sweep(vector<BlockState>& states, vector<MEntries>& m_entries_r,
                Vprop wr, Vprop b, VLprop cv, VVprop vmap, Vprop clabel,
                vector<int>& vlist, bool deg_corr, bool dense, bool multigraph,
                double beta, Eprop eweight, Vprop vweight, Graph& g,
                bool sequential, bool parallel, bool random_move, double c,
                size_t nmerges, Vprop merge_map, size_t niter, size_t B,
                bool verbose, RNG& rng, double& S, size_t& nmoves)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

    if (vlist.size() < 100)
        parallel = false;

    nmoves = 0;
    S = 0;

    vector<rng_t*> rngs;
    size_t num_threads = 1;
    if (parallel)
    {
#ifdef USING_OPENMP
        num_threads = omp_get_max_threads();
#endif
        for (size_t i = 0; i < num_threads; ++i)
        {
            std::array<int, std::mt19937::state_size> seed_data;
            std::generate_n(seed_data.data(), seed_data.size(), std::ref(rng));
            std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
            rngs.push_back(new rng_t(seq));
        }
    }
    else
    {
        rngs.push_back(&rng);
    }

    // used only if merging
    std::unordered_set<vertex_t> past_moves;
    vector<pair<vertex_t, double> > best_move;
    if (nmerges > 0 || parallel)
        best_move.resize(num_vertices(g), make_pair(vertex_t(0), numeric_limits<double>::max()));

    std::uniform_int_distribution<size_t> s_rand(0, B - 1);

    vector<MEntries> m_entries = m_entries_r;

    for (size_t iter = 0; iter < niter; ++iter)
    {
        if (nmerges == 0 && !parallel)
            std::shuffle(vlist.begin(), vlist.end(), rng);

        int i = 0, N = vlist.size();
        #pragma omp parallel for default(shared) private(i) \
            firstprivate(past_moves, m_entries) \
            schedule(runtime) if (parallel)
        for (i = 0; i < N; ++i)
        {
            size_t tid = 0;
            if (parallel)
            {
    #ifdef USING_OPENMP
                tid = omp_get_thread_num();
    #endif
            }

            typedef std::uniform_real_distribution<> rdist_t;
            auto rand_real = std::bind(rdist_t(), std::ref(*rngs[tid]));

            vertex_t v;
            if (sequential)
            {
                v = vertex(vlist[i], g);
            }
            else
            {
                std::uniform_int_distribution<size_t> v_rand(0, N - 1);
                v = vertex(vlist[v_rand(*rngs[tid])], g);
            }

            vertex_t r = b[v];

            // blocks can't become empty (if not merging)
            if (nmerges == 0 && wr[r] == vweight[v] && std::isinf(beta))
                continue;

            if (vweight[v] == 0)
                continue;

            if (nmerges > 0)
                past_moves.clear();

            // attempt random block
            vertex_t s = s_rand(*rngs[tid]);

            if (!random_move && total_degreeS()(v, g) > 0)
            {
                auto& state = states[0];

                if (nmerges == 0)
                {
                    vertex_t u = state.neighbour_sampler[v].sample(*rngs[tid]);

                    vertex_t t = state.b[u];

                    double p_rand = 0;
                    if (c > 0)
                    {
                        if (is_directed::apply<Graph>::type::value)
                            p_rand = c * B / double(state.mrp[t] + state.mrm[t] + c * B);
                        else
                            p_rand = c * B / double(state.mrp[t] + c * B);
                    }

                    if (c == 0 || rand_real() >= p_rand)
                    {
                        const auto& e = egroups_manage::sample_edge(state.egroups[t], *rngs[tid]);
                        s = state.b[target(e, state.g)];
                        if (s == t)
                            s = state.b[source(e, state.g)];
                    }
                }
                else
                {
                    // neighbour sampler points to the *block graph*
                    s = state.neighbour_sampler[r].sample(*rngs[tid]);
                    if (s == r)
                        s = state.cavity_neighbour_sampler[s].sample(*rngs[tid]);
                    else
                        s = state.neighbour_sampler[s].sample(*rngs[tid]);
                }
            }

            if (s == r)
                continue;

            if (wr[s] == 0 && std::isinf(beta)) // don't populate empty blocks
                continue;

            if (clabel[s] != clabel[r])
                continue;

            if (nmerges > 0)
            {
                if (wr[s] == 0)
                    continue;
                if (past_moves.find(s) != past_moves.end())
                    continue;
                past_moves.insert(s);
            }

            double dS = virtual_move(v, s, b, cv, vmap, states, m_entries,
                                     dense, deg_corr, multigraph);

            if (nmerges > 0)
            {
                if (dS < best_move[v].second)
                {
                    best_move[v].first = s;
                    best_move[v].second = dS;
                }
            }
            else
            {
                bool accept = false;
                if (std::isinf(beta))
                {
                    accept = dS < 0;
                }
                else
                {
                    double pf = 0, pb = 0;
                    if (random_move)
                    {
                        pf = pb = 1;
                    }
                    else
                    {
                        auto& state = states[0];
                        double p = get_move_prob(v, r, s, c, state.b,
                                                 state.mrs, state.mrp,
                                                 state.mrm, state.emat,
                                                 state.eweight, state.g,
                                                 state.bg, m_entries[0],
                                                 false,
                                                 state.overlap_stats);
                        pf += p;

                        p = get_move_prob(v, s, r, c, state.b,
                                          state.mrs, state.mrp, state.mrm,
                                          state.emat, state.eweight,
                                          state.g, state.bg,
                                          m_entries[0], true,
                                          state.overlap_stats);
                        pb += p;
                    }

                    double a = -beta * dS + log(pb) - log(pf);

                    if (a > 0)
                    {
                        accept = true;
                    }
                    else
                    {
                        double sample = rand_real();
                        accept = sample < exp(a);
                    }
                }

                if (accept)
                {
                    if (!parallel)
                    {
                        assert(b[v] == int(r));
                        move_vertex(v, s, b, cv, vmap, deg_corr, states,
                                    not random_move);

                        S += dS;
                        ++nmoves;
                        assert(b[v] == int(s));
                        if (verbose)
                            cout << v << ": " << r << " -> " << s << " " << S << " "
                                 << vlist.size() << " " << wr[r] << " " << wr[s] << endl;
                    }
                    else
                    {
                        best_move[v].first = s;
                        best_move[v].second = dS;
                    }
                }
            }
        }

        if (parallel && (nmerges == 0))
        {
            for (vertex_t v : vlist)
            {
                vertex_t r = b[v];
                vertex_t s = best_move[v].first;
                double dS = best_move[v].second;
                if (s != r && dS != numeric_limits<double>::max())
                {
                    dS = virtual_move(v, s, b, cv, vmap, states,
                                      m_entries, dense, deg_corr,
                                      multigraph);

                    if (dS > 0 && std::isinf(beta))
                        continue;

                    move_vertex(v, s, b, cv, vmap, deg_corr, states,
                                not random_move);

                    S += dS;
                    ++nmoves;
                }
            }
        }
    }

    if (parallel && (nmerges == 0))
    {
        for (auto r : rngs)
            delete r;
        rngs.clear();
    }

    if (nmerges > 0)
    {
        auto merge_cmp_less =
            [] (const std::tuple<vertex_t, vertex_t, double>& a,
                const std::tuple<vertex_t, vertex_t, double>& b) -> bool
            {
                return get<2>(a) < get<2>(b);
            };

        // top is the merge with _largest_ dS
        priority_queue<std::tuple<vertex_t, vertex_t, double>,
                       vector<std::tuple<vertex_t, vertex_t, double> >,
                       decltype(merge_cmp_less)> move_heap(merge_cmp_less);

        for (size_t i = 0; i < num_vertices(g); ++i)
        {
            vertex_t v = vertex(i, g);
            vertex_t r = b[v];
            vertex_t s = best_move[v].first;
            double dS = best_move[v].second;

            if (r != s && dS < numeric_limits<double>::max() &&
                (move_heap.size() < nmerges * 2 || dS < get<2>(move_heap.top())))
            {
                move_heap.emplace(v, s, dS);
                if (move_heap.size() > nmerges * 2)
                    move_heap.pop();
            }

            // if (r != s && dS < numeric_limits<double>::max())
            // {
            //     assert(vweight[v] > 0);
            //     move_heap.push(std::make_tuple(v, s, dS));
            // }
        }

        // back is the merge with _smallest_ dS
        vector<pair<vertex_t, vertex_t>> best_moves;
        best_moves.reserve(move_heap.size());

        while (!move_heap.empty())
        {
            auto& v = move_heap.top();
            best_moves.emplace_back(get<0>(v), get<1>(v));
            move_heap.pop();
        }

        vector<bool> touched(B, false);

        while (!best_moves.empty() && nmoves < nmerges)
        {
            vertex_t v = best_moves.back().first;
            vertex_t s = best_moves.back().second;
            best_moves.pop_back();

            vertex_t r = b[v];

            if (touched[r] || touched[s] || vweight[v] == 0 || wr[s] == 0)
                continue;

            double dS = virtual_move(v, s, b, cv, vmap, states, m_entries,
                                     dense, deg_corr, multigraph);
            move_vertex(v, s, b, cv, vmap, deg_corr, states, false);

            size_t i = 0, j = 0;
            auto l_v = cv[v][0];
            auto l_s = cv[s][0];
            while (i <  cv[v].size() && j < cv[s].size())
            {
                for (; i < cv[v].size(); ++i)
                {
                    l_v = cv[v][i];
                    if (l_v >= l_s)
                        break;
                    insert_vec(cv[s], j, l_v);
                    insert_vec(vmap[s], j, vmap[v][i]);
                    ++j;
                }

                for (; j < cv[s].size(); ++j)
                {
                    l_s = cv[s][j];
                    if (l_s >= l_v)
                        break;
                }

                if (l_s != l_v)
                    continue;

                if (i < cv[v].size() && j < cv[s].size())
                {
                    assert(l_s == l_v);
                    size_t u = vmap[v][i];
                    size_t w = vmap[s][j];
                    int l = l_v + 1;
                    auto& state = states[l];
                    assert(u < num_vertices(state.g));
                    assert(w < num_vertices(state.g));
                    merge_vertices(u, w, state.eweight, state.vweight, state.g);
                    i++;
                    j++;
                }
            }
            merge_vertices(v, s, eweight, vweight, g);

            clear_vec(cv[v]);
            clear_vec(vmap[v]);

            merge_map[v] = s;

            touched[r] = touched[s] = true;
            ++nmoves;

            S += dS;
        }

        // collapse merge tree across multiple calls
        for (auto v : vertices_range(g))
        {
            vertex_t u = merge_map[v];
            while (merge_map[u] != int(u))
                u = merge_map[u];
            merge_map[v] = u;
            b[v] = b[u];
        }
    }
}


template <class Vertex, class Graph, class Eprop, class SMap>
void build_neighbour_sampler(Vertex v, SMap& sampler, Eprop& eweight,
                             bool self_loops, Graph& g)
{
    vector<Vertex> neighbours;
    vector<double> probs;
    neighbours.reserve(total_degreeS()(v, g));
    probs.reserve(total_degreeS()(v, g));

    for (auto e : all_edges_range(v, g))
    {
        Vertex u = target(e, g);
        if (is_directed::apply<Graph>::type::value && u == v)
            u = source(e, g);
        if (!self_loops && u == v)
            continue;
        neighbours.push_back(u);
        probs.push_back(eweight[e]);  // Self-loops will be added twice, and
                                      // hence will be sampled with probability
                                      // 2 * eweight[e]
    }

    if (probs.empty())
    {
        neighbours.push_back(v);
        probs.push_back(1.);
    }

    sampler[v] = Sampler<Vertex, mpl::false_>(neighbours, probs);
};

struct init_neighbour_sampler
{
    template <class Graph, class Eprop>
    void operator()(Graph& g, Eprop eweight, bool self_loops,
                    bool empty, boost::any& asampler) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

        typedef typename property_map<Graph, vertex_index_t>::type vindex_map_t;
        typedef typename property_map_type::apply<Sampler<vertex_t, mpl::false_>,
                                                  vindex_map_t>::type::unchecked_t
            sampler_map_t;

        sampler_map_t sampler(get(vertex_index_t(), g), num_vertices(g));
        asampler = sampler;

        if (!empty)
        {
            for (auto v : vertices_range(g))
                build_neighbour_sampler(v, sampler, eweight, self_loops, g);
        }
    }
};


// Sampling marginal probabilities on the edges
template <class Graph, class Vprop, class MEprop>
void collect_edge_marginals(size_t B, Vprop b, MEprop p, Graph& g, Graph&)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

    typename graph_traits<Graph>::edge_iterator e, e_end;
    for(tie(e, e_end) = edges(g); e != e_end; ++e)
    {
        vertex_t u = min(source(*e, g), target(*e, g));
        vertex_t v = max(source(*e, g), target(*e, g));

        vertex_t r = b[u];
        vertex_t s = b[v];

        typename property_traits<MEprop>::value_type& pv = p[*e];
        if (pv.size() < B * B)
            pv.resize(B * B);
        size_t j = r + B * s;
        pv[j]++;
    }
}

struct bethe_entropy
{
    template <class Graph, class MEprop, class MVprop>
    void operator()(Graph& g, size_t B, MEprop p, MVprop pv, double& H,
                    double& sH, double& Hmf, double& sHmf) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

        typename graph_traits<Graph>::vertex_iterator v, v_end;
        for(tie(v, v_end) = vertices(g); v != v_end; ++v)
        {
            pv[*v].resize(B);
            for (size_t i = 0; i < B; ++i)
                pv[*v][i] = 0;
        }

        H = Hmf = sH = sHmf =  0;

        typename graph_traits<Graph>::edge_iterator e, e_end;
        for(tie(e, e_end) = edges(g); e != e_end; ++e)
        {
            vertex_t u = min(source(*e, g), target(*e, g));
            vertex_t v = max(source(*e, g), target(*e, g));

            double sum = 0;
            for (size_t r = 0; r < B; ++r)
                for (size_t s = 0; s < B; ++s)
                {
                    size_t i = r + B * s;
                    pv[u][r] += p[*e][i];
                    pv[v][s] += p[*e][i];
                    sum += p[*e][i];
                }

            for (size_t i = 0; i < B * B; ++i)
            {
                if (p[*e][i] == 0)
                    continue;
                double pi = double(p[*e][i]) / sum;
                H -= pi * log(pi);
                sH += pow((log(pi) + 1) * sqrt(pi / sum), 2);
            }
        }

        for(tie(v, v_end) = vertices(g); v != v_end; ++v)
        {
            double sum = 0;
            for (size_t i = 0; i < B; ++i)
                sum += pv[*v][i];
            for (size_t i = 0; i < B; ++i)
            {
                if (pv[*v][i] == 0)
                    continue;
                pv[*v][i] /= sum;
                double pi = pv[*v][i];
                double kt = (1 - double(in_degreeS()(*v, g)) - double(out_degree(*v, g)));
                if (kt != 0)
                {
                    H -= kt * (pi * log(pi));
                    sH += pow(kt * (log(pi) + 1) * sqrt(pi / sum), 2);
                }

                Hmf -= pi * log(pi);
                sHmf += pow((log(pi) + 1) * sqrt(pi / sum), 2);
            }
        }
    }
};

template <class Graph, class Vprop, class VVprop>
void collect_vertex_marginals(Vprop b, VVprop p, Graph& g)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

    typename graph_traits<Graph>::vertex_iterator v, v_end;
    for(tie(v, v_end) = vertices(g); v != v_end; ++v)
    {
        vertex_t r = b[*v];
        if (p[*v].size() <= r)
            p[*v].resize(r + 1);
        p[*v][r]++;
    }
}


} // graph_tool namespace

#endif //GRAPH_BLOCKMODEL_HH
