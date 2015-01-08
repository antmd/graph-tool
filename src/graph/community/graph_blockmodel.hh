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
    assert(x < __safelog_cache.size());
    return __safelog_cache[x];
}

__attribute__((always_inline))
inline double xlogx(size_t x)
{
    if (x >= __xlogx_cache.size())
        cout << x << " " << __xlogx_cache.size();
    assert(x < __xlogx_cache.size());
    return __xlogx_cache[x];
}

__attribute__((always_inline))
inline double lgamma_fast(size_t x)
{
    if (x >= __lgamma_cache.size())
        return lgamma(x);
    return __lgamma_cache[x];
}

template <class Type>
Type poisson_entropy(Type l, Type epsilon=1e-8)
{
    if (l == 0)
        return 0;

    Type S = l * (1 - log(l));
    Type delta = 1 + epsilon;
    Type ll = log(l);
    int k = 2;

    Type x = 0;
    while (delta > epsilon)
    {
        Type old = x;
        x += exp(k * ll + log(lgamma(k + 1)) - lgamma(k + 1));
        k++;
        delta = abs(x - old);
    }

    S += exp(-l + log(x));
    return S;
}

template <class Type>
inline Type lpoisson(Type l, int k)
{
    if (l == 0)
        return (k == 0) ? 0 : -numeric_limits<Type>::infinity();
    return -l + k * log(l) - lgamma(k + 1);
}

template <class Type>
inline Type poisson(Type l, int k)
{
    return exp(lpoisson(l, k));
}

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
    get_mu_l(N, E, mu, l);
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

//
// Sparse entropy
//

// "edge" term of the entropy
template <class Graph>
__attribute__((always_inline))
inline double eterm(size_t r, size_t s, size_t mrs, const Graph& g)
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
template <class Graph, class Vertex>
inline double vterm(Vertex v, size_t mrp, size_t mrm, size_t wr, bool deg_corr,
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
            S += vterm(v, mrp[v], mrm[v], wr[v], deg_corr, g);
    }
};


template <class Vertex, class List, class Graph, class GetNode>
double get_parallel_neighbours_entropy(Vertex v, List& us, Graph& g,
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


//
// Dense entropy
//

__attribute__((always_inline))
inline double lbinom(double N, double k)
{
    if (N == 0 || k == 0 || k > N)
        return 0;
    return lgamma(N + 1) - lgamma(N - k + 1) - lgamma(k + 1);
}

__attribute__((always_inline))
inline double lbinom_fast(int N, int k)
{
    if (N == 0 || k == 0 || k > N)
        return 0;
    return lgamma_fast(N + 1) - lgamma_fast(N - k + 1) - lgamma_fast(k + 1);
}


// "edge" term of the entropy
template <class Graph>
__attribute__((always_inline))
inline double eterm_dense(size_t r, size_t s, int ers, double wr_r,
                          double wr_s, bool multigraph, const Graph& g)
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
        S = lbinom(nrns + ers - 1, ers);
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

// ===============================
// Partition stats
// ===============================

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
    partition_stats_t(Graph& g, Vprop b, Eprop eweight, size_t N, size_t B)
        : _enabled(true), _N(N), _B(B), _hist(B), _total(B), _ep(B), _em(B)
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
            _hist[r][make_pair(in_degreeS()(v, g), out_degree(v, g))]++;
            _total[r]++;
            _em[r] += in_degreeS()(v, g);
            _ep[r] += out_degreeS()(v, g);
        }
    }

    double get_partition_dl()
    {
        double S = 0;
        S += lbinom(_B + _N - 1, _N);
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
    double get_delta_dl(size_t v, size_t r, size_t nr, bool deg_corr, OStats&,
                        Graph& g)
    {
        if (r == nr)
            return 0;

        double S_b = 0, S_a = 0;
        S_b += -lgamma_fast(_total[r] + 1) - lgamma_fast(_total[nr] + 1);
        S_a += -lgamma_fast(_total[r]    ) - lgamma_fast(_total[nr] + 2);

        if (deg_corr)
        {
            int kin = in_degreeS()(v, g);
            int kout = out_degreeS()(v, g);

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
        }

        return S_a - S_b;
    }

    template <class Graph, class OStats>
    void move_vertex(size_t v, size_t r, size_t nr,
                     bool deg_corr, OStats& overlap_stats, Graph& g)
    {
        if (r == nr)
            return;

        _total[r]--;
        _total[nr]++;

        if (deg_corr)
        {
            auto deg = make_pair(in_degreeS()(v, g), out_degree(v, g));
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

private:
    bool _enabled;
    size_t _N;
    size_t _B;
    vector<map_t> _hist;
    vector<int> _total;
    vector<int> _ep;
    vector<int> _em;

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
       const typename get_emat_t::apply<Graph>::type& emat, const Graph& g)
{
    return emat[r][s];
}

template <class Graph>
inline __attribute__((always_inline))
void
put_me(typename graph_traits<Graph>::vertex_descriptor r,
       typename graph_traits<Graph>::vertex_descriptor s,
       const typename graph_traits<Graph>::edge_descriptor& e,
       typename get_emat_t::apply<Graph>::type& emat, const Graph& g)
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
          const typename graph_traits<Graph>::edge_descriptor& e,
          typename get_emat_t::apply<Graph>::type& emat, Graph& g,
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
        typedef dense_hash_map<vertex_t, edge_t> map_t;
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
       const typename get_ehash_t::apply<Graph>::type& ehash, const Graph& bg)
{
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
       const Graph& bg)
{
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
        typedef typename get_ehash_t::apply<Graph>::map_t map_t;
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

        emat_t emat(num_vertices(g));

#ifdef HAVE_SPARSEHASH
        for (auto v : vertices_range(g))
        {
            emat[v].set_empty_key(numeric_limits<vertex_t>::max());
            emat[v].set_deleted_key(num_vertices(g));
        }
#endif

        for (auto e : edges_range(g))
            put_me(source(e, g), target(e, g), e, emat, g);

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
                   Graph& g, BGraph& bg, EMat& emat,
                   OStats& overlap_stats,
                   const NPolicy& npolicy = NPolicy())
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    vertex_t r = b[v];

    int self_count = 0;
    for (auto e : npolicy.get_out_edges(v, g))
    {
        vertex_t u = target(e, g);
        if (u == v && !is_directed::apply<Graph>::type::value)
        {
            ++self_count;
            if (self_count % 2 == 0)
                continue;
        }

        vertex_t s = b[u];

        // if (!emat[r][s].second)
        //     throw GraphException("no edge? " + lexical_cast<string>(r) +
        //                          " " + lexical_cast<string>(s));

        typename graph_traits<BGraph>::edge_descriptor me =
            get_me(r, s, emat, bg).first;

        size_t ew = eweight[e];
        mrs[me] -= ew;

        mrp[r] -= ew;
        mrm[s] -= ew;

        if (mrs[me] == 0)
            remove_me(r, s, me, emat, bg);
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
        overlap_stats.remove_half_edge(v, r, g);
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

    int self_count = 0;
    for (auto e : npolicy.get_out_edges(v, g))
    {
        vertex_t u = target(e, g);
        if (u == v && !is_directed::apply<Graph>::type::value )
        {
            ++self_count;
            if (self_count % 2 == 0)
                continue;
        }
        vertex_t s;

        if (u != v)
            s = b[u];
        else
            s = r;

        typename graph_traits<BGraph>::edge_descriptor me;

        pair<typename graph_traits<BGraph>::edge_descriptor, bool> mep =
                get_me(r, s, emat, bg);

        if (!mep.second)
        {
            mep = add_edge(r, s, bg);
            put_me(r, s, mep.first, emat, bg);
            mrs[mep.first] = 0;
        }
        me = mep.first;

        size_t ew = eweight[e];

        mrs[me] += ew;

        mrp[r] += ew;
        mrm[s] += ew;
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
        overlap_stats.add_half_edge(v, r, g);
        wr[r] = overlap_stats.get_block_size(r);
    }

    b[v] = r;
}


// move a vertex from its current block to block nr
template <class Graph, class BGraph, class Eprop, class Vprop, class EWprop,
          class VWprop, class EMat, class OStats, class PStats,
          class NPolicy = standard_neighbours_policy>
void move_vertex(size_t v, size_t nr, Eprop& mrs, Vprop& mrp, Vprop& mrm,
                 Vprop& wr, Vprop& b, bool deg_corr, const EWprop& eweight,
                 const VWprop& vweight, Graph& g, BGraph& bg, EMat& emat,
                 OStats& overlap_stats,  PStats& partition_stats,
                 const NPolicy& npolicy = NPolicy())
{
    if (b[v] == int(nr))
        return;

    if (partition_stats.is_enabled())
        partition_stats.move_vertex(v, b[v], nr, deg_corr, overlap_stats, g);

    remove_vertex(v, mrs, mrp, mrm, wr, b, eweight, vweight, g, bg, emat,
                  overlap_stats, npolicy);
    add_vertex(v, nr, mrs, mrp, mrm, wr, b, eweight, vweight, g, bg, emat,
               overlap_stats, npolicy);
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
    }

    void SetMove(size_t r, size_t nr)
    {
        _rnr = make_pair(r, nr);
    }

    void InsertDeltaTarget(size_t r, size_t s, int delta)
    {
        if (!is_directed::apply<Graph>::type::value &&
            (s == _rnr.first || s == _rnr.second) && s < r)
        {
            InsertDeltaTarget(s, r, delta);
            return;
        }

        if (is_directed::apply<Graph>::type::value &&
            (s == _rnr.first || s == _rnr.second) && s < r)
        {
            InsertDeltaSource(r, s, delta);
            return;
        }

        vector<size_t>& field = (_rnr.first == r) ? _r_field_t : _nr_field_t;
        if (field[s] == _null)
        {
            field[s] = _entries.size();
            _entries.push_back(make_pair(r, s));
            _delta.push_back(delta);
        }
        else
        {
            _delta[field[s]] += delta;
        }
    }

    void InsertDeltaSource(size_t s, size_t r, int delta)
    {
        if (s == r)
        {
            InsertDeltaTarget(r, s, delta);
            return;
        }

        if ((s == _rnr.first || s == _rnr.second) && s < r)
        {
            InsertDeltaTarget(s, r, delta);
            return;
        }

        vector<size_t>& field = (_rnr.first == r) ? _r_field_s : _nr_field_s;
        if (field[s] == _null)
        {
            field[s] = _entries.size();
            _entries.push_back(make_pair(s, r));
            _delta.push_back(delta);
        }
        else
        {
            _delta[field[s]] += delta;
        }
    }

    int GetDelta(size_t t, size_t s)
    {
        if (is_directed::apply<Graph>::type::value)
        {
            if (t == _rnr.first || t == _rnr.second)
                return GetDeltaTarget(t, s);
            if (s == _rnr.first || s == _rnr.second)
                return GetDeltaSource(t, s);
        }
        else
        {
            if (t == _rnr.first || t == _rnr.second)
                return GetDeltaTarget(t, s);
            if (s == _rnr.first || s == _rnr.second)
                return GetDeltaTarget(s, t);
        }
        return 0;
    }

    int GetDeltaTarget(size_t r, size_t s)
    {
        vector<size_t>& field = (_rnr.first == r) ? _r_field_t : _nr_field_t;
        if (field[s] == _null)
        {
            return 0;
        }
        else
        {
            return _delta[field[s]];
        }
    }

    int GetDeltaSource(size_t s, size_t r)
    {
        vector<size_t>& field = (_rnr.first == r) ? _r_field_s : _nr_field_s;
        if (field[s] == _null)
        {
            return 0;
        }
        else
        {
            return _delta[field[s]];
        }
    }

    void Clear()
    {
        size_t r, s;
        for (size_t i = 0; i < _entries.size(); ++i)
        {
            tie(r, s) = _entries[i];
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

    vector<pair<size_t, size_t> >& GetEntries() { return _entries; }
    vector<int>& GetDelta() { return _delta; }

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

// obtain the necessary entries in the e_rs matrix which need to be modified
// after the move
template <class Graph, class BGraph, class Vertex, class Vprop, class Eprop,
          class NPolicy = standard_neighbours_policy>
void move_entries(Vertex v, Vertex nr, Vprop& b, Eprop& eweights, Graph& g,
                  BGraph& bg, EntrySet<Graph>& m_entries,
                  const NPolicy& npolicy = NPolicy())
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    vertex_t r = b[v];

    m_entries.SetMove(r, nr);

    int self_count = 0;
    for (auto e : npolicy.get_out_edges(v, g))
    {
        vertex_t u = target(e, g);
        if (u == v && !is_directed::apply<Graph>::type::value)
        {
            ++self_count;
            if (self_count % 2 == 0)
                continue;
        }
        vertex_t s = b[u];
        int ew = eweights[e];
        //assert(ew > 0);

        m_entries.InsertDeltaTarget(r, s, -ew);

        if (u == v)
            s = nr;

        m_entries.InsertDeltaTarget(nr, s, +ew);
    }

    for (auto e : npolicy.get_in_edges(v, g))
    {
        vertex_t u = source(e, g);
        if (u == v)
            continue;
        vertex_t s = b[u];
        int ew = eweights[e];

        m_entries.InsertDeltaSource(s,  r, -ew);
        m_entries.InsertDeltaSource(s, nr, +ew);
    }
}

// obtain the entropy difference given a set of entries in the e_rs matrix
template <class Graph, class Eprop, class BGraph, class EMat>
double entries_dS(EntrySet<Graph>& m_entries, Eprop& mrs, BGraph& bg, EMat& emat)
{
    typedef typename graph_traits<BGraph>::vertex_descriptor vertex_t;
    vector<pair<size_t, size_t> >& entries = m_entries.GetEntries();
    vector<int>& delta = m_entries.GetDelta();

    double dS = 0;
    for (size_t i = 0; i < entries.size(); ++i)
    {
        vertex_t er = entries[i].first;
        vertex_t es = entries[i].second;
        int d = delta[i];

        int ers = get_mrs(er, es, mrs, emat, bg);
        assert(ers + d >= 0);
        dS += eterm(er, es, ers + d, bg) - eterm(er, es, ers, bg);
    }
    return dS;
}

// compute the entropy difference of a virtual move of vertex r to block nr
template <class Graph, class BGraph, class Eprop, class Vprop, class EWprop,
          class VWprop, class EMat, class OStats,
          class NPolicy = standard_neighbours_policy>
double virtual_move_sparse(size_t v, size_t nr, Eprop& mrs, Vprop& mrp,
                           Vprop& mrm, Vprop& wr, Vprop& b, bool deg_corr,
                           const EWprop& eweight, const VWprop& vweight,
                           Graph& g, BGraph& bg, EMat& emat,
                           EntrySet<Graph>& m_entries, OStats& overlap_stats,
                           bool parallel_edges,
                           const NPolicy& npolicy = NPolicy())

{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    vertex_t r = b[v];

    if (r == nr)
        return 0.;

    m_entries.Clear();
    move_entries(v, nr, b, eweight, g, bg, m_entries, npolicy);
    double dS = entries_dS(m_entries, mrs, bg, emat);
    int kout = npolicy.get_out_degree(v, g, eweight);
    int kin = kout;
    if (is_directed::apply<Graph>::type::value)
        kin = npolicy.get_in_degree(v, g, eweight);
    //assert(mrm[r]  - kin >= 0);
    //assert(mrp[r]  - kout >= 0);

    int dwr, dwnr;
    set<int> bv, nbv;
    if (!overlap_stats.is_enabled())
    {
        dwr = dwnr = vweight[v];
    }
    else
    {
        dwr = wr[r] - overlap_stats.virtual_remove_size(v, r, g);
        dwnr = overlap_stats.virtual_add_size(v, nr, g) - wr[nr];

        if (deg_corr)
            dS += overlap_stats.virtual_move_dS(v, r, nr, g);

        if (parallel_edges)
            dS += overlap_stats.virtual_move_parallel_dS(v, r, nr, b, g);
    }

    dS += vterm(r,  mrp[r]  - kout, mrm[r]  - kin, wr[r]  - dwr , deg_corr, bg);
    dS += vterm(nr, mrp[nr] + kout, mrm[nr] + kin, wr[nr] + dwnr, deg_corr, bg);
    dS -= vterm(r,  mrp[r]        , mrm[r]       , wr[r]        , deg_corr, bg);
    dS -= vterm(nr, mrp[nr]       , mrm[nr]      , wr[nr]       , deg_corr, bg);

    return dS;
}


// compute the entropy difference of a virtual move of vertex r to block nr
template <class Graph, class BGraph, class Eprop, class Vprop, class EWprop,
          class VWprop, class EMat, class OStats,
          class NPolicy = standard_neighbours_policy>
double virtual_move_dense(size_t v, size_t nr, Eprop& mrs, Vprop& mrp,
                          Vprop& mrm, Vprop& wr, Vprop& b, bool deg_corr,
                          const EWprop& eweight, const VWprop& vweight,
                          Graph& g, BGraph& bg, EMat& emat, bool multigraph,
                          OStats& overlap_stats, const NPolicy& npolicy = NPolicy())
{
    if (deg_corr)
        throw GraphException("Dense entropy for degree corrected model not implemented!");

    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    vertex_t r = b[v];

    if (r == nr)
        return 0;

    vector<int> deltap(num_vertices(bg), 0);
    int deltal = 0;
    int self_count = 0;
    for (auto e : npolicy.get_out_edges(v, g))
    {
        vertex_t u = target(e, g);
        vertex_t s = b[u];
        if (u == v)
        {
            ++self_count;
            if (self_count % 2 == 0)
                continue;
            deltal += eweight[e];
        }
        else
        {
            deltap[s] += eweight[e];
        }
    }

    vector<int> deltam(num_vertices(bg), 0);
    for (auto e : npolicy.get_in_edges(v, g))
    {
        vertex_t u = source(e, g);
        if (u == v)
            continue;
        vertex_t s = b[u];
        deltam[s] += eweight[e];
    }

    double dS = 0;
    int dwr, dwnr;
    if (!overlap_stats.is_enabled())
    {
        dwr = dwnr = vweight[v];
    }
    else
    {
        dwr = wr[r] - overlap_stats.virtual_remove_size(v, r, g);
        dwnr = overlap_stats.virtual_add_size(v, nr, g) - wr[nr];

        if (deg_corr)
            dS += overlap_stats.virtual_move_dS(v, r, nr, g);
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

template <class Graph, class BGraph, class Eprop, class Vprop, class EWprop,
          class VWprop, class EMat, class OStats, class PStats,
          class NPolicy = standard_neighbours_policy>
double virtual_move(size_t v, size_t nr, bool dense, Eprop& mrs, Vprop& mrp,
                    Vprop& mrm, Vprop& wr, Vprop& b, bool deg_corr,
                    const EWprop& eweight, const VWprop& vweight, Graph& g,
                    BGraph& bg, EMat& emat, EntrySet<Graph>& m_entries,
                    OStats& overlap_stats, bool parallel_edges,
                    PStats& partition_stats, const NPolicy& npolicy = NPolicy())
{
    double S = 0;
    if (dense)
    {
        m_entries.Clear();
        move_entries(v, nr, b, eweight, g, bg, m_entries, npolicy); // expected!
        S = virtual_move_dense(v, nr, mrs, mrp, mrm, wr, b, deg_corr,
                               eweight, vweight, g, bg, emat,
                               parallel_edges, overlap_stats, npolicy);
    }
    else
    {
        S = virtual_move_sparse(v, nr, mrs, mrp, mrm, wr, b, deg_corr,
                                eweight, vweight, g, bg, emat, m_entries,
                                overlap_stats, parallel_edges,
                                npolicy);
    }

    if (partition_stats.is_enabled())
        S += partition_stats.get_delta_dl(v, b[v], nr, deg_corr,  overlap_stats, g);

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
                      VertexIndex vertex_index, size_t B, bool weighted, bool empty)
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
                           esrcpos, etgtpos, eweight, g, vertex_index, B,
                           mpl::true_());
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
                           esrcpos, etgtpos, eweight, g, vertex_index, B,
                           mpl::true_());
        }
    }

    template <class Eprop, class Vprop, class VEprop, class Graph, class VertexIndex, class Egroups>
    static void build_dispatch(Vprop b, Egroups egroups, VEprop esrcpos,
                               VEprop etgtpos, Eprop eweight, Graph& g,
                               VertexIndex vertex_index, size_t B, mpl::true_)
    {
        for (auto e : edges_range(g))
        {
            size_t r = b[get_source(e, g)];
            assert (r < B);
            auto& r_elist = egroups[r];
            esrcpos[e] = insert_edge(std::make_tuple(e, true), r_elist, eweight[e]);

            size_t s = b[get_target(e, g)];
            assert (s < B);
            auto& s_elist = egroups[s];
            etgtpos[e] = insert_edge(std::make_tuple(e, false), s_elist, eweight[e]);
        }
    }

    template <class Edge, class EV>
    static size_t insert_edge(const Edge& e, EV& elist, size_t weight)
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
                            DynamicSampler<Edge>& elist)
    {
        elist.remove(pos);
    }

    template <class Edge, class Epos>
    static void remove_edge(size_t pos, Epos& esrcpos, Epos& etgtpos,
                            vector<Edge>& elist)
    {
        if (get<1>(elist.back()))
            esrcpos[get<0>(elist.back())] = pos;
        else
            etgtpos[get<0>(elist.back())] = pos;
        elist[pos] = elist.back();
        elist.pop_back();
    }

    template <class Vertex, class Graph, class EVprop, class Eprop, class VEprop>
    static void remove_egroups(Vertex v, Vertex r, Eprop& eweight,
                               EVprop& egroups, VEprop& esrcpos, VEprop& etgtpos,
                               Graph& g)
    {
        typedef Vertex vertex_t;
        int self_count = 0;

        //update the half-edge lists
        for (auto e : all_edges_range(v, g))
        {
            vertex_t src = get_source(e, g);
            vertex_t tgt = get_target(e, g);

            bool is_src = (src == v);

            // self-loops will appear twice
            if (src == tgt)
            {
                is_src = self_count % 2 == 0;
                ++self_count;
            }

            auto& elist = egroups[r];
            size_t pos = (is_src) ? esrcpos[e] : etgtpos[e];
            remove_edge(pos, esrcpos, etgtpos, elist);
        }
    }

    template <class Vertex, class Graph, class EVprop, class Eprop, class VEprop>
    static void add_egroups(Vertex v, Vertex s, Eprop& eweight, EVprop& egroups,
                            VEprop& esrcpos, VEprop& etgtpos, Graph& g)
    {
        typedef Vertex vertex_t;
        int self_count = 0;

        //update the half-edge lists
        for (auto e : all_edges_range(v, g))
        {
            vertex_t src = get_source(e, g);
            vertex_t tgt = get_target(e, g);

            bool is_src = (src == v);

            // self-loops will appear twice
            if (src == tgt)
            {
                is_src = self_count % 2 == 0;
                ++self_count;
            }

            auto& elist = egroups[s];
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
template <class Vertex, class Graph, class Vprop, class Eprop, class Emat,
          class BGraph, class OStats>
inline double
get_move_prob(Vertex v, Vertex r, Vertex s, double c, Vprop& b, Eprop& mrs,
              Vprop& mrp, Vprop& mrm, Emat& emat, Eprop& eweight, Graph& g,
              BGraph& bg, EntrySet<Graph>& m_entries, bool reverse,
              OStats& overlap_stats)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    size_t B = num_vertices(bg);
    double p = 0;
    size_t w = 0;

    int kout = 0, kin = 0;
    if (reverse)
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
                int dts = m_entries.GetDelta(t, s);
                int dst = dts;
                if (is_directed::apply<Graph>::type::value)
                    dst = m_entries.GetDelta(s, t);

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
void merge_vertices(size_t u, size_t v, Eprop& eweight_u, Vprop& vweight,
                    Vprop& merge_map, Graph& g)
{
    auto eweight = eweight_u.get_checked();

    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    int self_count = 0;
    for(auto e : out_edges_range(u, g))
    {
        vertex_t t = target(e, g);
        if (t == u && !is_directed::apply<Graph>::type::value)
        {
            ++self_count;
            if (self_count % 2 == 0)
                continue;
        }

        if (t == u)
            t = v;

        auto ne = edge(v, t, g);
        if (!ne.second)
        {
            ne = add_edge(v, t, g);
            eweight[ne.first] = 0;
        }
        eweight[ne.first] += eweight[e];
    }

    for (auto e : in_edges_range(u, g))
    {
        vertex_t t = source(e, g);
        if (t == u)
            continue;

        auto ne = edge(t, v, g);
        if (!ne.second)
        {
            ne = add_edge(t, v, g);
            eweight[ne.first] = 0;
        }
        eweight[ne.first] += eweight[e];
    }

    vweight[v] += vweight[u];
    vweight[u] = 0;
    merge_map[u] = v;
    clear_vertex(u, g);
}

struct merge_cmp_less
{
    template<class Vertex>
    double operator()(const std::tuple<Vertex, Vertex, double>& a,
                      const std::tuple<Vertex, Vertex, double>& b)
    {
        return get<2>(a) < get<2>(b);
    }
};

//A single Monte Carlo Markov chain sweep
template <class Graph, class BGraph, class EMprop, class Eprop, class Vprop,
          class EMat, class EVprop, class VEprop, class SamplerMap, class OStats,
          class RNG>
void move_sweep(EMprop mrs, Vprop mrp, Vprop mrm, Vprop wr, Vprop b,
                Vprop clabel, vector<int>& vlist, bool deg_corr, bool dense,
                bool multigraph, double beta, Eprop eweight, Vprop vweight,
                EVprop egroups, VEprop esrcpos, VEprop etgtpos, Graph& g,
                BGraph& bg, EMat& emat, SamplerMap neighbour_sampler,
                SamplerMap cavity_neighbour_sampler, bool sequential,
                bool parallel, bool random_move, double c, size_t nmerges,
                size_t ntries, Vprop merge_map,
                partition_stats_t& partition_stats, bool verbose, RNG& rng,
                double& S, size_t& nmoves, OStats overlap_stats)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

    size_t B = num_vertices(bg);

    if (vlist.size() < 100)
        parallel = false;

    // used only if merging
    std::unordered_set<vertex_t> past_moves;
    vector<pair<vertex_t, double> > best_move;
    if (nmerges > 0 || parallel)
    {
        best_move.resize(num_vertices(g), make_pair(vertex_t(0), numeric_limits<double>::max()));
    }
    else
    {
        std::shuffle(vlist.begin(), vlist.end(), rng);
    }

    std::uniform_int_distribution<size_t> s_rand(0, B - 1);

    nmoves = 0;
    S = 0;

    EntrySet<Graph> m_entries(B);

    vector<rng_t*> rngs;
    if (parallel)
    {
        size_t num_threads = 1;
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

    int i = 0, N = vlist.size();
    #pragma omp parallel for default(shared) private(i) \
        firstprivate(m_entries, past_moves) schedule(runtime) if (parallel)
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

        size_t j = 0;
        while (j < ntries)
        {
            ++j;

            // attempt random block
            vertex_t s = s_rand(*rngs[tid]);

            if (!random_move && total_degreeS()(v, g) > 0)
            {
                if (nmerges == 0)
                {
                    vertex_t u = neighbour_sampler[v].sample(*rngs[tid]);

                    vertex_t t = b[u];

                    double p_rand = 0;
                    if (c > 0)
                    {
                        if (is_directed::apply<Graph>::type::value)
                            p_rand = c * B / double(mrp[t] + mrm[t] + c * B);
                        else
                            p_rand = c * B / double(mrp[t] + c * B);
                    }

                    if (c == 0 || rand_real() >= p_rand)
                    {
                        const auto& e = egroups_manage::sample_edge(egroups[t], *rngs[tid]);
                        s = b[target(e, g)];
                        if (s == t)
                            s = b[source(e, g)];
                    }
                }
                else
                {
                    // neighbour sampler points to the *block graph*
                    s = neighbour_sampler[r].sample(*rngs[tid]);
                    if (s == r)
                        s = cavity_neighbour_sampler[s].sample(*rngs[tid]);
                    else
                        s = neighbour_sampler[s].sample(*rngs[tid]);
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

            double dS = virtual_move(v, s, dense, mrs, mrp, mrm, wr, b,
                                     deg_corr, eweight, vweight, g, bg, emat,
                                     m_entries, overlap_stats, multigraph,
                                     partition_stats);

            if (nmerges > 0)
            {
                if (dS < best_move[v].second)
                {
                    best_move[v].first = s;
                    best_move[v].second = dS;
                    j = 0; // reset try counter
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
                    double pf = random_move ? 1 :
                        get_move_prob(v, r, s, c, b, mrs, mrp, mrm, emat,
                                      eweight, g, bg, m_entries, false,
                                      overlap_stats);

                    double pb = random_move ? 1 :
                        get_move_prob(v, s, r, c, b, mrs, mrp, mrm, emat,
                                      eweight, g, bg, m_entries, true,
                                      overlap_stats);

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
                        move_vertex(v, s, mrs, mrp, mrm, wr, b, deg_corr,
                                    eweight, vweight, g, bg, emat,
                                    overlap_stats, partition_stats);
                        if (!random_move)
                            egroups_manage::update_egroups(v, r, s, eweight, egroups,
                                                           esrcpos, etgtpos, g);
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
    }

    if (parallel && (nmerges == 0))
    {
        for (auto r : rngs)
            delete r;
        rngs.clear();

        for (vertex_t v : vlist)
        {
            vertex_t r = b[v];
            vertex_t s = best_move[v].first;
            double dS = best_move[v].second;
            if (s != r && dS != numeric_limits<double>::max())
            {
                double dS = virtual_move(v, s, dense, mrs, mrp, mrm, wr, b,
                                         deg_corr, eweight, vweight, g, bg,
                                         emat, m_entries, overlap_stats,
                                         multigraph, partition_stats);
                if (dS > 0 && std::isinf(beta))
                    continue;
                move_vertex(v, s, mrs, mrp, mrm, wr, b, deg_corr, eweight,
                            vweight, g, bg, emat, overlap_stats,
                            partition_stats);
                if (!random_move)
                    egroups_manage::update_egroups(v, r, s, eweight, egroups,
                                                   esrcpos, etgtpos, g);
                S += dS;
                ++nmoves;
            }
        }
    }

    if (nmerges > 0)
    {
        // top is the merge with _largest_ dS
        priority_queue<std::tuple<vertex_t, vertex_t, double>,
                       vector<std::tuple<vertex_t, vertex_t, double> >,
                       merge_cmp_less> move_heap;

        for (size_t i = 0; i < num_vertices(g); ++i)
        {
            vertex_t v = vertex(i, g);
            vertex_t r = b[v];
            vertex_t s = best_move[v].first;
            double dS = best_move[v].second;

            if (r != s && dS < numeric_limits<double>::max() &&
                (move_heap.size() < nmerges || dS < get<2>(move_heap.top())))
            {
                move_heap.push(std::make_tuple(v, s, dS));
                if (move_heap.size() > nmerges)
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
            std::tuple<vertex_t, vertex_t, double> v = move_heap.top();
            best_moves.push_back(make_pair(get<0>(v), get<1>(v)));
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

            double dS = virtual_move(v, s, dense, mrs, mrp, mrm, wr, b,
                                     deg_corr, eweight, vweight, g, bg, emat,
                                     m_entries, overlap_stats, multigraph,
                                     partition_stats);

            move_vertex(v, s, mrs, mrp, mrm, wr, b, deg_corr, eweight, vweight,
                        g, bg, emat, overlap_stats, partition_stats);

            merge_vertices(v, s, eweight, vweight, merge_map, g);

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
