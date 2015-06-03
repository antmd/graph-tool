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

#ifndef GRAPH_BLOCKMODEL_OVERLAP_HH
#define GRAPH_BLOCKMODEL_OVERLAP_HH

#include "config.h"
#include <unordered_set>
#include <unordered_map>
#include <tuple>

#include "graph_blockmodel.hh"

#ifdef HAVE_SPARSEHASH
#include SPARSEHASH_INCLUDE(dense_hash_set)
#include SPARSEHASH_INCLUDE(dense_hash_map)
#endif

namespace graph_tool
{
#ifdef HAVE_SPARSEHASH
using google::dense_hash_set;
using google::dense_hash_map;
#endif

using namespace boost;

//===============
// Overlap stats
//===============

class overlap_stats_t
{
public:
    typedef pair<size_t, size_t> deg_t;

    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type::unchecked_t
        vmap_t;
    typedef property_map_type::apply<int64_t,
                                     GraphInterface::vertex_index_map_t>::type::unchecked_t
        vimap_t;
    typedef property_map_type::apply<vector<int64_t>,
                                     GraphInterface::vertex_index_map_t>::type::unchecked_t
        vvmap_t;

    overlap_stats_t(): _enabled(false) {}

    template <class Graph>
    overlap_stats_t(Graph& g, vmap_t b, vvmap_t half_edges, vimap_t node_index,
                    size_t B)
        : _half_edges(half_edges), _node_index(node_index),
          _out_neighbours(num_vertices(g), -1),
          _in_neighbours(num_vertices(g), -1)
    {
        _enabled = true;
        _block_nodes.resize(B);

#ifdef HAVE_SPARSEHASH
        for (auto& bnodes : _block_nodes)
        {
            bnodes.set_empty_key(numeric_limits<size_t>::max());
            bnodes.set_deleted_key(numeric_limits<size_t>::max() - 1);
        }
#endif

        size_t N = 0;
        for (auto v : vertices_range(g))
        {
            size_t vi = node_index[v];
            N = std::max(N, vi + 1);
            size_t kin = in_degreeS()(v, g);
            size_t kout = out_degreeS()(v, g);

            size_t r = b[v];
            auto& bnodes = _block_nodes[r];
            auto& k = bnodes[vi];
            k.first += kin;
            k.second += kout;

            for (auto e : out_edges_range(v, g))
                _out_neighbours[v] = target(e, g);
            for (auto e : in_edges_range(v, g))
                _in_neighbours[v] = source(e, g);
        }

        // parallel edges
        _mi.resize(num_vertices(g), -1);

        for (size_t i = 0; i < N; ++i)
        {
            auto& he = half_edges[i];

            unordered_map<size_t, vector<size_t>> out_us;
            for (auto u : he)
            {
                auto w = _out_neighbours[u];
                if (w == -1)
                    continue;
                if (!is_directed::apply<Graph>::type::value && size_t(node_index[w]) < i)
                    continue;
                out_us[node_index[w]].push_back(u);
            }

            for (auto& uc : out_us)
            {
                if (uc.second.size() > 1)
                {
                    _parallel_bundles.resize(_parallel_bundles.size() + 1);
                    auto& h = _parallel_bundles.back();
#ifdef HAVE_SPARSEHASH
                    h.set_empty_key(make_pair(numeric_limits<size_t>::max(),
                                              numeric_limits<size_t>::max()));
                    h.set_deleted_key(make_pair(numeric_limits<size_t>::max() - 1,
                                                numeric_limits<size_t>::max() - 1));
#endif
                    for (auto u : uc.second)
                    {
                        auto w = _out_neighbours[u];
                        assert(w != -1);
                        _mi[u] = _mi[w] = _parallel_bundles.size() - 1;
                        size_t r = b[u];
                        size_t s = b[w];
                        if (!is_directed::apply<Graph>::type::value && r > s)
                            std::swap(r, s);
                        h[make_pair(r, s)]++;
                    }
                }
            }
        }
    }

    template <class Graph, class VProp>
    void add_half_edge(size_t v, size_t v_r, VProp& b, Graph&)
    {
        size_t u = _node_index[v];
        size_t kin = (_in_neighbours[v] == -1) ? 0 : 1;
        size_t kout = (_out_neighbours[v] == -1) ? 0 : 1;
        assert(kin + kout == 1);
        auto& k = _block_nodes[v_r][u];
        k.first += kin;
        k.second += kout;

        int m = _mi[v];
        if (m != -1)
        {
            size_t r, s;
            int u = _out_neighbours[v];
            if (u == -1)
            {
                u = _in_neighbours[v];
                r = b[u];
                s = v_r;
            }
            else
            {
                r = v_r;
                s = b[u];
            }
            auto& h = _parallel_bundles[m];
            if (!is_directed::apply<Graph>::type::value && r > s)
                std::swap(r, s);
            h[make_pair(r, s)]++;
        }
    }

    template <class Graph, class VProp>
    void remove_half_edge(size_t v, size_t v_r, VProp& b, Graph&)
    {
        size_t u = _node_index[v];
        size_t kin = (_in_neighbours[v] == -1) ? 0 : 1;
        size_t kout = (_out_neighbours[v] == -1) ? 0 : 1;
        assert(kin + kout == 1);
        auto& k = _block_nodes[v_r][u];
        k.first -= kin;
        k.second -= kout;

        if (k.first + k.second == 0)
            _block_nodes[v_r].erase(u);

        int m = _mi[v];
        if (m != -1)
        {
            size_t r, s;
            int u = _out_neighbours[v];
            if (u == -1)
            {
                u = _in_neighbours[v];
                r = b[u];
                s = v_r;
            }
            else
            {
                r = v_r;
                s = b[u];
            }
            auto& h = _parallel_bundles[m];
            if (!is_directed::apply<Graph>::type::value && r > s)
                std::swap(r, s);
            auto iter = h.find(make_pair(r, s));
            assert(iter->second > 0);
            iter->second--;
            if (iter->second == 0)
                h.erase(iter);
        }
    }

    size_t get_block_size(size_t r) const
    {
        return _block_nodes[r].size();
    }

    size_t virtual_remove_size(size_t v, size_t r, size_t in_deg = 0,
                               size_t out_deg = 0) const
    {
        size_t nr = _block_nodes[r].size();
        size_t u = _node_index[v];
        size_t kin = (in_deg + out_deg) > 0 ?
            in_deg : ((_in_neighbours[v] == -1) ? 0 : 1);
        size_t kout = (in_deg + out_deg) > 0 ?
            out_deg : ((_out_neighbours[v] == -1) ? 0 : 1);
        const auto iter = _block_nodes[r].find(u);
        const auto& deg = iter->second;
        if (deg.first == kin && deg.second == kout)
            nr--;
        return nr;
    }

    size_t virtual_add_size(size_t v, size_t r) const
    {
        size_t nr = _block_nodes[r].size();
        size_t u = _node_index[v];
        const auto& bnodes = _block_nodes[r];
        if (bnodes.find(u) == bnodes.end())
            nr++;
        return nr;
    }

    template <class Graph>
    double virtual_move_dS(size_t v, size_t r, size_t nr, Graph& g,
                           size_t in_deg = 0, size_t out_deg = 0) const
    {
        double dS = 0;

        size_t u = _node_index[v];
        size_t u_kin = ((in_deg + out_deg) > 0) ? in_deg : in_degreeS()(v, g);
        size_t u_kout = ((in_deg + out_deg) > 0) ? out_deg : out_degreeS()(v, g);

        auto deg =  _block_nodes[r].find(u)->second;
        auto ndeg = deg;
        ndeg.first -= u_kin;
        ndeg.second -= u_kout;

        dS -= lgamma_fast(ndeg.first + 1) + lgamma_fast(ndeg.second + 1);
        dS += lgamma_fast(deg.first + 1) + lgamma_fast(deg.second + 1);

        const auto iter = _block_nodes[nr].find(u);
        if (iter != _block_nodes[nr].end())
            deg = iter->second;
        else
            deg = make_pair(0, 0);
        ndeg = deg;
        ndeg.first += u_kin;
        ndeg.second += u_kout;

        dS -= lgamma_fast(ndeg.first + 1) + lgamma_fast(ndeg.second + 1);
        dS += lgamma_fast(deg.first + 1) + lgamma_fast(deg.second + 1);

        return dS;
    }

    template <class Graph, class VProp>
    double virtual_move_parallel_dS(size_t v, size_t v_r, size_t v_nr, VProp& b,
                                    Graph&, bool coherent=false) const
    {
        int m = _mi[v];
        if (m == -1)
            return 0;

        size_t r, s, nr, ns;
        int u = _out_neighbours[v];
        if (u == -1)
        {
            u = _in_neighbours[v];
            r = b[u];
            s = v_r;
            nr = r;
            ns = v_nr;
        }
        else
        {
            r = v_r;
            s = b[u];
            nr = v_nr;
            ns = s;
        }
        if (!is_directed::apply<Graph>::type::value && r > s)
            std::swap(r, s);
        if (!is_directed::apply<Graph>::type::value && nr > ns)
            std::swap(nr, ns);

        auto& h = _parallel_bundles[m];

        auto get_h = [&](const pair<size_t, size_t>& k) -> int
            {
                const auto iter = h.find(k);
                if (iter == h.end())
                    return 0;
                return iter->second;
            };

        int c  = get_h(make_pair(r,  s));
        int nc = get_h(make_pair(nr, ns));

        assert(c > 0);
        assert(nc >= 0);
        assert(v_r != v_nr);
        assert(make_pair(r, s) != make_pair(nr, ns));

        double S = 0;
        S -= lgamma_fast(c + 1) + lgamma_fast(nc + 1);
        if (!coherent)
            S += lgamma_fast(c) + lgamma_fast(nc + 2);
        else
            S += lgamma_fast(c + nc + 1);

        return S;
    }

    // sample another half-edge adjacent to the node w
    template <class RNG>
    size_t sample_half_edge(size_t w, RNG& rng) const
    {
        auto& half_edges = _half_edges[w];
        return uniform_sample(half_edges, rng);
    }

    size_t get_node(size_t v) const { return _node_index[v]; }
    const vector<int64_t>& get_half_edges(size_t v) const { return _half_edges[v]; }

    int64_t get_out_neighbour(size_t v) const { return _out_neighbours[v]; }
    int64_t get_in_neighbour(size_t v) const { return _in_neighbours[v]; }


    bool is_enabled() const { return _enabled; }


#ifdef HAVE_SPARSEHASH
    typedef dense_hash_map<pair<size_t, size_t>, int,
                           std::hash<pair<size_t, size_t>>> phist_t;
#else
    typedef unordered_map<pair<size_t, size_t>, int> phist_t;
#endif

    const vector<phist_t>& get_parallel_bundles() { return _parallel_bundles; }
    const vector<int>& get_mi() { return _mi; }

private:
    bool _enabled;
    vvmap_t _half_edges;     // half-edges to each node
    vimap_t _node_index;     // node to each half edges

#ifdef HAVE_SPARSEHASH
    typedef dense_hash_map<size_t, deg_t, std::hash<size_t>> node_map_t;
#else
    typedef unordered_map<size_t, deg_t> node_map_t;
#endif

    vector<node_map_t> _block_nodes; // nodes (and degrees) in each block

    vector<int64_t> _out_neighbours;
    vector<int64_t> _in_neighbours;


    vector<int> _mi;
    vector<phist_t> _parallel_bundles; // parallel edge bundles
};


//=============================
// Partition Description length
//=============================

extern double lbinom(double N, double k);
extern double xlogx(size_t x);

struct overlap_partition_stats_t
{
    typedef std::tuple<int, int> deg_t;
    typedef vector<deg_t> cdeg_t;

    typedef vector<int> bv_t;

#ifdef HAVE_SPARSEHASH
    typedef dense_hash_map<bv_t, size_t, std::hash<bv_t>> bhist_t;
    typedef dense_hash_map<cdeg_t, size_t, std::hash<cdeg_t>> cdeg_hist_t;

    typedef dense_hash_map<bv_t, cdeg_hist_t, std::hash<bv_t>> deg_hist_t;

    typedef dense_hash_map<bv_t, vector<size_t>, std::hash<bv_t>> ebhist_t;
#else
    typedef unordered_map<bv_t, size_t> bhist_t;
    typedef unordered_map<cdeg_t, size_t, std::hash<cdeg_t>> cdeg_hist_t;

    typedef unordered_map<bv_t, cdeg_hist_t> deg_hist_t;

    typedef unordered_map<bv_t, vector<size_t>> ebhist_t;
#endif

    typedef unordered_map<int, int> dmap_t;

    overlap_partition_stats_t() : _enabled(false) {}

    template <class Graph, class Vprop, class Eprop>
    overlap_partition_stats_t(Graph& g, Vprop b, overlap_stats_t& overlap_stats,
                              Eprop eweight, size_t N, size_t B, bool edges_dl)
        : _enabled(true), _N(N), _B(B), _D(0),
          _dhist(B + 1), _r_count(B), _bhist(N), _emhist(B), _ephist(B),
          _embhist(N), _epbhist(N), _deg_hist(N), _bvs(N), _nbvs(N), _degs(N),
          _ndegs(N), _deg_delta(B), _edges_dl(edges_dl)
    {
#ifdef HAVE_SPARSEHASH
        bv_t empty = {numeric_limits<int>::max()};
        bv_t deleted = {numeric_limits<int>::max() - 1};
        _bhist.set_empty_key(empty);
        _bhist.set_deleted_key(deleted);
        _deg_hist.set_empty_key(empty);
        _deg_hist.set_deleted_key(deleted);
        _embhist.set_empty_key(empty);
        _embhist.set_deleted_key(deleted);
        _epbhist.set_empty_key(empty);
        _epbhist.set_deleted_key(deleted);
#endif

        dmap_t in_hist, out_hist;
        for (size_t v = 0; v < N; ++v)
        {
            dmap_t in_hist, out_hist;
            set<size_t> rs;

            get_bv_deg(v, b, eweight, overlap_stats, g, rs, in_hist, out_hist);

            cdeg_t cdeg;
            for (auto r : rs)
            {
                deg_t deg = std::make_tuple(in_hist[r], out_hist[r]);
                cdeg.push_back(deg);
            }

            bv_t bv(rs.begin(), rs.end());

            _bvs[v] = bv;
            _bvs[v].reserve(bv.size() * 2);
            _nbvs[v] = bv;
            _nbvs[v].reserve(bv.size() * 2);
            _degs[v] = cdeg;
            _degs[v].reserve(cdeg.size() * 2);
            _ndegs[v] = cdeg;
            _ndegs[v].reserve(cdeg.size() * 2);

            auto & cdh = _deg_hist[bv];
#ifdef HAVE_SPARSEHASH
            if (cdh.empty())
            {
                cdeg_t empty = {deg_t(numeric_limits<int>::max(),
                                      numeric_limits<int>::max())};
                cdeg_t deleted = {deg_t(numeric_limits<int>::max() - 1,
                                        numeric_limits<int>::max() - 1)};
                cdh.set_empty_key(empty);
                cdh.set_deleted_key(deleted);
            }

#endif
            cdh[cdeg]++;

            size_t d = bv.size();
            _D = max(_D, d);
            _dhist[d]++;
            _bhist[bv]++;

            auto& bmh = _embhist[bv];
            auto& bph = _epbhist[bv];
            bmh.resize(bv.size());
            bph.resize(bv.size());

            for (size_t i = 0; i < bv.size(); ++i)
            {
                size_t r = bv[i];
                _emhist[r] += get<0>(cdeg[i]);
                _ephist[r] += get<1>(cdeg[i]);
                bmh[i] += get<0>(cdeg[i]);
                bph[i] += get<1>(cdeg[i]);
            }

        }

        for (auto& bv_c : _bhist)
        {
            assert(bv_c.second > 0);
            for (auto r : bv_c.first)
                _r_count[r] += 1;
        }

        _actual_B = 0;
        for (size_t r = 0; r < _B; ++r)
        {
            if (overlap_stats.get_block_size(r) > 0)
                _actual_B++;
        }
    }

    template <class Graph, class Vprop, class Eprop>
    void get_bv_deg(size_t v, Vprop& b, Eprop& eweight,
                    overlap_stats_t& overlap_stats, Graph& g, set<size_t>& rs,
                    dmap_t& in_hist, dmap_t& out_hist) const
    {
        if (!overlap_stats.is_enabled())
        {
            in_hist[b[v]] += in_degreeS()(v, g, eweight);
            out_hist[b[v]] += out_degreeS()(v, g, eweight);
        }
        else
        {
            auto& half_edges = overlap_stats.get_half_edges(v);
            for (size_t u : half_edges)
            {
                size_t kin = in_degreeS()(u, g);
                size_t kout = out_degreeS()(u, g);

                in_hist[b[u]] += kin;
                out_hist[b[u]] += kout;
            }
        }

        for (auto& rk : in_hist)
            rs.insert(rk.first);
    }


    double get_partition_dl() const
    {
        double S = 0;
        for (size_t d = 1; d < _dhist.size(); ++d)
        {
            size_t nd = _dhist[d];
            if (nd == 0)
                continue;
            double x = lbinom_fast(_actual_B, d);
            double ss = lbinom_careful((exp(x) + nd) - 1, nd); // not fast
            if (std::isinf(ss) || std::isnan(ss))
                ss = nd * x - lgamma_fast(nd + 1);
            assert(!std::isinf(ss));
            assert(!std::isnan(ss));
            S += ss;
        }

        S += lbinom_fast(_D + _N - 1, _N) + lgamma_fast(_N + 1);

        for (auto& bh : _bhist)
            S -= lgamma_fast(bh.second + 1);

        return S;
    }

    double get_deg_dl(bool ent, bool deg_alt, bool xi_fast) const
    {
        double S = 0;
        if (ent)
        {
            for (auto& ch : _deg_hist)
            {
                auto& bv = ch.first;
                auto& cdeg_hist = ch.second;

                size_t n_bv = _bhist.find(bv)->second;

                S += xlogx(n_bv);
                for (auto& dh : cdeg_hist)
                    S -= xlogx(dh.second);
            }
        }
        else
        {
            S = 0;

            for (auto& ch : _deg_hist)
            {
                auto& bv = ch.first;
                auto& cdeg_hist = ch.second;

                size_t n_bv = _bhist.find(bv)->second;

                if (n_bv == 0)
                    continue;

                const auto& bmh = _embhist.find(bv)->second;
                const auto& bph = _epbhist.find(bv)->second;

                double S1 = 0;
                for (size_t i = 0; i < bv.size(); ++i)
                {
                    if (xi_fast)
                    {
                        S1 += get_xi_fast(n_bv, bmh[i]);
                        S1 += get_xi_fast(n_bv, bph[i]);
                    }
                    else
                    {
                        S1 += get_xi(n_bv, bmh[i]);
                        S1 += get_xi(n_bv, bph[i]);
                    }
                }

                S1 += lgamma_fast(n_bv + 1);

                for (auto& dh : cdeg_hist)
                    S1 -= lgamma_fast(dh.second + 1);

                if (deg_alt)
                {
                    double S2 = 0;
                    for (size_t i = 0; i < bv.size(); ++i)
                    {
                        S2 += lbinom(n_bv + bmh[i] - 1, bmh[i]);
                        S2 += lbinom(n_bv + bph[i] - 1, bph[i]);
                    }
                    S += min(S1, S2);
                }
                else
                {
                    S += S1;
                }
            }

            for (size_t r = 0; r < _B; ++r)
            {
                if (_r_count[r] == 0)
                    continue;
                S += lbinom(_r_count[r] + _emhist[r] - 1,  _emhist[r]);
                S += lbinom(_r_count[r] + _ephist[r] - 1,  _ephist[r]);
            }
        }
        return S;
    }


    template <class Graph>
    bool get_n_bv(size_t v, size_t r, size_t nr, const bv_t& bv,
                  const cdeg_t& deg, bv_t& n_bv, cdeg_t& n_deg, Graph& g,
                  size_t in_deg = 0, size_t out_deg = 0)
    {
        size_t kin = (in_deg + out_deg == 0) ? in_degreeS()(v, g) : in_deg;
        size_t kout = (in_deg + out_deg == 0) ? out_degreeS()(v, g) : out_deg;

        auto& d_r = _deg_delta[r];
        d_r.first -= kin;
        d_r.second -= kout;

        auto& d_nr = _deg_delta[nr];
        d_nr.first += kin;
        d_nr.second += kout;

        n_deg.clear();
        n_bv.clear();
        bool is_same_bv = true;
        bool has_r = false, has_nr = false;
        for (size_t i = 0; i < bv.size(); ++i)
        {
            size_t s = bv[i];
            auto k_s = deg[i];

            auto& d_s = _deg_delta[s];
            get<0>(k_s) += d_s.first;
            get<1>(k_s) += d_s.second;

            d_s.first = d_s.second = 0;

            if (s == r)
                has_r = true;

            if (s == nr)
                has_nr = true;

            if ((get<0>(k_s) + get<1>(k_s)) > 0)
            {
                n_bv.push_back(s);
                n_deg.push_back(k_s);
            }
            else
            {
                is_same_bv = false;
            }
        }

        if (!has_r || !has_nr)
        {
            is_same_bv = false;
            std::array<size_t, 2> ss = {r, nr};
            for (auto s : ss)
            {
                auto& d_s = _deg_delta[s];
                if (d_s.first + d_s.second == 0)
                    continue;
                size_t kin = d_s.first;
                size_t kout = d_s.second;
                auto pos = std::lower_bound(n_bv.begin(), n_bv.end(), s);
                auto dpos = n_deg.begin();
                std::advance(dpos, pos - n_bv.begin());
                n_bv.insert(pos, s);
                n_deg.insert(dpos, make_pair(kin, kout));
                d_s.first = d_s.second = 0;
            }
        }
        return is_same_bv;
    }

    // get deg counts without increasing the container
    size_t get_deg_count(const bv_t& bv, const cdeg_t& deg) const
    {
        auto iter = _deg_hist.find(bv);
        if (iter == _deg_hist.end())
            return 0;
        auto& hist = iter->second;
        if (hist.empty())
            return 0;
        auto diter = hist.find(deg);
        if (diter == hist.end())
            return 0;
        return diter->second;
    }

    // get bv counts without increasing the container
    size_t get_bv_count(const bv_t& bv) const
    {
        auto iter = _bhist.find(bv);
        if (iter == _bhist.end())
            return 0;
        return iter->second;
    }

    template <class Graph>
    double get_delta_dl(size_t v, size_t r, size_t nr,
                        overlap_stats_t& overlap_stats, Graph& g,
                        size_t in_deg = 0, size_t out_deg = 0)
    {
        if (r == nr)
            return 0;

        size_t u = overlap_stats.get_node(v);
        auto& bv = _bvs[u];
        auto& n_bv = _nbvs[u];
        size_t d = bv.size();
        cdeg_t& deg = _degs[u];
        cdeg_t& n_deg = _ndegs[u];

        bool is_same_bv = get_n_bv(v, r, nr, bv, deg, n_bv, n_deg, g, in_deg,
                                   out_deg);

        size_t n_d = n_bv.size();
        size_t n_D = _D;

        if (d == _D && _dhist[d] == 1)
        {
            n_D = 1;
            for (auto& bc : _bhist)
            {
                if (bc.first.size() == d || bc.second == 0)
                    continue;
                n_D = max(n_D, bc.first.size());
            }
        }

        n_D = max(n_D, n_d);

        double S_a = 0, S_b = 0;

        if (n_D != _D)
        {
            S_b += lbinom(_D  + _N - 1, _N);
            S_a += lbinom(n_D + _N - 1, _N);
        }

        int dB = 0;
        if (overlap_stats.virtual_remove_size(v, r, in_deg, out_deg) == 0)
            dB--;
        if (overlap_stats.get_block_size(nr) == 0)
            dB++;

        auto get_S_d = [&] (size_t d_i, int delta, int dB) -> double
            {
                int nd = int(_dhist[d_i]) + delta;
                if (nd == 0)
                    return 0.;
                double x = lbinom_fast(_actual_B + dB, d_i);
                double S = lbinom(exp(x) + nd - 1, nd); // not fast
                if (std::isinf(S) || std::isnan(S))
                    S = nd * x - lgamma_fast(nd + 1);
                return S;
            };

        if (dB == 0)
        {
            if (n_d != d)
            {
                S_b += get_S_d(d,  0, 0) + get_S_d(n_d, 0, 0);
                S_a += get_S_d(d, -1, 0) + get_S_d(n_d, 1, 0);
            }
        }
        else
        {
            for (size_t di = 0; di < min(_actual_B + abs(dB) + 1, _dhist.size()); ++di)
            {
                if (d != n_d)
                {
                    if (di == d)
                    {
                        S_b += get_S_d(d,  0, 0);
                        S_a += get_S_d(d, -1, dB);
                        continue;
                    }
                    if (di == n_d)
                    {
                        S_b += get_S_d(n_d, 0, 0);
                        S_a += get_S_d(n_d, 1, dB);
                        continue;
                    }
                }
                if (_dhist[di] == 0)
                    continue;
                S_b += get_S_d(di, 0, 0);
                S_a += get_S_d(di, 0, dB);
            }

            if (_edges_dl)
            {
                auto get_x = [](size_t B) -> size_t
                {
                    if (is_directed::apply<Graph>::type::value)
                        return B * B;
                    else
                        return (B * (B + 1)) / 2;
                };

                size_t E = num_vertices(g) / 2;
                S_b += lbinom(get_x(_actual_B) + E - 1, E);
                S_a += lbinom(get_x(_actual_B + dB) + E - 1, E);
            }

        }

        size_t bv_count = get_bv_count(bv);
        size_t n_bv_count = get_bv_count(n_bv);

        auto get_S_b = [&] (bool is_bv, int delta) -> double
            {
                if (is_bv)
                    return -lgamma_fast(bv_count + delta + 1);
                return -lgamma_fast(n_bv_count + delta + 1);
            };

        if (!is_same_bv)
        {
            S_b += get_S_b(true,  0) + get_S_b(false, 0);
            S_a += get_S_b(true, -1) + get_S_b(false, 1);
        }

        return S_a - S_b;
    }


    template <class Graph, class EWeight>
    double get_delta_deg_dl(size_t v, size_t r, size_t nr, EWeight&,
                            overlap_stats_t& overlap_stats, Graph& g,
                            size_t in_deg = 0, size_t out_deg = 0)
    {
        if (r == nr)
            return 0;

        double S_b = 0, S_a = 0;

        size_t u = overlap_stats.get_node(v);
        auto& bv = _bvs[u];
        auto& n_bv = _nbvs[u];

        cdeg_t& deg = _degs[u];
        cdeg_t& n_deg = _ndegs[u];

        bool is_same_bv = get_n_bv(v, r, nr, bv, deg, n_bv, n_deg, g, in_deg,
                                   out_deg);

        size_t bv_count = get_bv_count(bv);
        size_t n_bv_count = get_bv_count(n_bv);

        auto get_S_bv = [&] (bool is_bv, int delta) -> double
            {
                if (is_bv)
                    return lgamma_fast(bv_count + delta + 1);
                return lgamma_fast(n_bv_count + delta + 1);
            };

        auto get_S_e = [&] (bool is_bv, int bdelta, int deg_delta) -> double
            {
                size_t bv_c = ((is_bv) ? bv_count : n_bv_count) + bdelta;
                if (bv_c == 0)
                    return 0.;

                const cdeg_t& deg_i = (is_bv) ? deg : n_deg;
                const auto& bv_i = (is_bv) ? bv : n_bv;

                double S = 0;
                if (((is_bv) ? bv_count : n_bv_count) > 0)
                {
                    const auto& bmh = _embhist.find(bv_i)->second;
                    const auto& bph = _epbhist.find(bv_i)->second;

                    assert(bmh.size() == bv_i.size());
                    assert(bph.size() == bv_i.size());

                    for (size_t i = 0; i < bv_i.size(); ++i)
                    {
                        S += get_xi_fast(bv_c, bmh[i] + deg_delta * int(get<0>(deg_i[i])));
                        S += get_xi_fast(bv_c, bph[i] + deg_delta * int(get<1>(deg_i[i])));
                    }
                }
                else
                {
                    for (size_t i = 0; i < bv_i.size(); ++i)
                    {
                        S += get_xi_fast(bv_c, deg_delta * int(get<0>(deg_i[i])));
                        S += get_xi_fast(bv_c, deg_delta * int(get<1>(deg_i[i])));
                    }
                }

                return S;
            };

        auto get_S_e2 = [&] (int deg_delta, int ndeg_delta) -> double
            {
                double S = 0;
                const auto& bmh = _embhist.find(bv)->second;
                const auto& bph = _epbhist.find(bv)->second;

                for (size_t i = 0; i < bv.size(); ++i)
                {
                    S += get_xi_fast(bv_count, bmh[i] +
                                     deg_delta * int(get<0>(deg[i])) +
                                     ndeg_delta * int(get<0>(n_deg[i])));
                    S += get_xi_fast(bv_count, bph[i] +
                                     deg_delta * int(get<1>(deg[i])) +
                                     ndeg_delta * int(get<1>(n_deg[i])));
                }
                return S;
            };

        if (!is_same_bv)
        {
            S_b += get_S_bv(true,  0) + get_S_bv(false, 0);
            S_a += get_S_bv(true, -1) + get_S_bv(false, 1);

            S_b += get_S_e(true,  0,  0) + get_S_e(false, 0, 0);
            S_a += get_S_e(true, -1, -1) + get_S_e(false, 1, 1);
        }
        else
        {
            S_b += get_S_e2( 0, 0);
            S_a += get_S_e2(-1, 1);
        }

        size_t deg_count = get_deg_count(bv, deg);
        size_t n_deg_count = get_deg_count(n_bv, n_deg);

        auto get_S_deg = [&] (bool is_deg, int delta) -> double
            {
                if (is_deg)
                    return -lgamma_fast(deg_count + delta + 1);
                return -lgamma_fast(n_deg_count + delta + 1);
            };

        S_b += get_S_deg(true,  0) + get_S_deg(false, 0);
        S_a += get_S_deg(true, -1) + get_S_deg(false, 1);

        auto is_in = [&] (const bv_t& bv, size_t r) -> bool
            {
                auto iter = lower_bound(bv.begin(), bv.end(), r);
                if (iter == bv.end())
                    return false;
                if (size_t(*iter) != r)
                    return false;
                return true;
            };

        for (size_t s : bv)
        {
            S_b += lbinom_fast(_r_count[s] + _emhist[s] - 1, _emhist[s]);
            S_b += lbinom_fast(_r_count[s] + _ephist[s] - 1, _ephist[s]);
        }

        for (size_t s : n_bv)
        {
            if (is_in(bv, s))
                continue;
            S_b += lbinom_fast(_r_count[s] + _emhist[s] - 1, _emhist[s]);
            S_b += lbinom_fast(_r_count[s] + _ephist[s] - 1, _ephist[s]);
        }

#ifdef HAVE_SPARSEHASH
        dense_hash_map<size_t, pair<int, int>, std::hash<size_t>> deg_delta;
        dense_hash_map<size_t, int, std::hash<size_t>> r_count_delta;
        deg_delta.set_empty_key(numeric_limits<size_t>::max());
        r_count_delta.set_empty_key(numeric_limits<size_t>::max());
#else
        unordered_map<size_t, pair<int, int>> deg_delta;
        unordered_map<size_t, int> r_count_delta;
#endif

        if (bv != n_bv)
        {
            if (n_bv_count == 0)
            {
                for (auto s : n_bv)
                    r_count_delta[s] += 1;
            }

            if (bv_count == 1)
            {
                for (auto s : bv)
                   r_count_delta[s] -= 1;
            }
        }

        if (r != nr)
        {
            size_t kin = (in_deg + out_deg == 0) ? in_degreeS()(v, g) : in_deg;
            size_t kout = (in_deg + out_deg == 0) ? out_degreeS()(v, g) : out_deg;

            auto& d_r = deg_delta[r];
            d_r.first -= kin;
            d_r.second -= kout;

            auto& d_nr = deg_delta[nr];
            d_nr.first += kin;
            d_nr.second += kout;
        }

        for (size_t s : bv)
        {
            S_a += lbinom_fast(_r_count[s] + r_count_delta[s] + _emhist[s] + deg_delta[s].first - 1,
                               _emhist[s] + deg_delta[s].first);
            S_a += lbinom_fast(_r_count[s] + r_count_delta[s] + _ephist[s] + deg_delta[s].second - 1,
                               _ephist[s] + deg_delta[s].second);
        }

        for (size_t s : n_bv)
        {
            if (!is_in(bv, s))
            {
                S_a += lbinom_fast(_r_count[s] + r_count_delta[s] + _emhist[s] + deg_delta[s].first - 1,
                                   _emhist[s] + deg_delta[s].first);
                S_a += lbinom_fast(_r_count[s] + r_count_delta[s] + _ephist[s] + deg_delta[s].second - 1,
                                   _ephist[s] + deg_delta[s].second);
            }
        }

        return S_a - S_b;
    }

    template <class Graph, class EWeight>
    void move_vertex(size_t v, size_t r, size_t nr, bool,
                     overlap_stats_t& overlap_stats, Graph& g, EWeight&,
                     size_t in_deg = 0, size_t out_deg = 0)
    {
        if (r == nr)
            return;

        size_t u = overlap_stats.get_node(v);
        auto& bv = _bvs[u];
        auto& n_bv = _nbvs[u];
        cdeg_t& deg = _degs[u];
        cdeg_t& n_deg = _ndegs[u];
        size_t d = bv.size();

        bool is_same_bv = get_n_bv(v, r, nr, bv, deg, n_bv, n_deg, g, in_deg,
                                   out_deg);
        size_t n_d = n_bv.size();

        if (!is_same_bv)
        {
            _dhist[d] -= 1;
            auto& bv_count = _bhist[bv];
            bv_count -= 1;

            if (bv_count == 0)
            {
                _bhist.erase(bv);
                for (auto s : bv)
                {
                    _r_count[s]--;
                    if (_r_count[s] == 0)
                        _actual_B--;
                }
            }

            if (d == _D && _dhist[d] == 0)
            {
                _D = 1;
                for (auto& bc : _bhist)
                {
                    if (bc.second == 0)
                        continue;
                    _D = max(_D, bc.first.size());
                }
            }

            _dhist[n_d] += 1;
            auto& n_bv_count = _bhist[n_bv];
            n_bv_count += 1;

            if (n_bv_count == 1)
            {
                for (auto s : n_bv)
                {
                    if (_r_count[s] == 0)
                        _actual_B++;
                    _r_count[s]++;
                }
            }

            if (n_d > _D)
                _D = n_d;
        }

        auto& deg_h = _deg_hist[bv];
        auto& deg_count = deg_h[deg];
        deg_count -= 1;
        if (deg_count == 0)
            deg_h.erase(deg);
        auto& bmh = _embhist[bv];
        auto& bph = _epbhist[bv];
        assert(bmh.size() == bv.size());
        assert(bph.size() == bv.size());
        for (size_t i = 0; i < bv.size(); ++i)
        {
            bmh[i] -= get<0>(deg[i]);
            bph[i] -= get<1>(deg[i]);
        }

        if (deg_h.empty())
        {
            _deg_hist.erase(bv);
            _embhist.erase(bv);
            _epbhist.erase(bv);
        }

        size_t kin = (in_deg + out_deg == 0) ? in_degreeS()(v, g) : in_deg;
        size_t kout = (in_deg + out_deg == 0) ? out_degreeS()(v, g) : out_deg;
        _emhist[r] -= kin;
        _ephist[r] -= kout;

        auto& hist = _deg_hist[n_bv];
#ifdef HAVE_SPARSEHASH
        if (hist.empty())
        {
            cdeg_t empty = {deg_t(numeric_limits<int>::max(),
                                  numeric_limits<int>::max())};
            cdeg_t deleted = {deg_t(numeric_limits<int>::max() - 1,
                                    numeric_limits<int>::max() - 1)};
            hist.set_empty_key(empty);
            hist.set_deleted_key(deleted);
        }
#endif
        hist[n_deg] += 1;
        auto& n_bmh = _embhist[n_bv];
        auto& n_bph = _epbhist[n_bv];
        n_bmh.resize(n_bv.size());
        n_bph.resize(n_bv.size());
        for (size_t i = 0; i < n_bv.size(); ++i)
        {
            n_bmh[i] += get<0>(n_deg[i]);
            n_bph[i] += get<1>(n_deg[i]);
        }

        _emhist[nr] += kin;
        _ephist[nr] += kout;

        _bvs[u].swap(n_bv);
        _degs[u].swap(n_deg);
    }

    bool is_enabled() const { return _enabled; }

private:
    bool _enabled;
public:
    size_t _N;
private:
    size_t _B;
    size_t _actual_B;
    size_t _D;
    vector<int> _dhist;        // d-histogram
    vector<int> _r_count;      // m_r
    bhist_t _bhist;            // b-histogram
    vector<size_t> _emhist;    // e-_r histogram
    vector<size_t> _ephist;    // e+_r histogram
    ebhist_t _embhist;         // e+^r_b histogram
    ebhist_t _epbhist;         // e+^r_b histogram
    deg_hist_t _deg_hist;      // n_k^b histogram
    vector<bv_t> _bvs;         // bv node map
    vector<bv_t> _nbvs;
    vector<cdeg_t> _degs;      // deg node map
    vector<cdeg_t> _ndegs;

    vector<pair<int, int>> _deg_delta;
    bool _edges_dl;
};

struct entropy_parallel_edges_overlap
{
    template <class Graph>
    void operator()(Graph&, overlap_stats_t& overlap_stats, double& S) const
    {
        S = 0;
        for(auto& h : overlap_stats.get_parallel_bundles())
        {
            for (auto& kc : h)
                S += lgamma_fast(kc.second + 1);
        }
    }
};

template <class Graph>
struct half_edge_neighbour_policy
{
    typedef typename graph_traits<Graph>::edge_descriptor edge_t;
    typedef edge_t* iter_t;

    half_edge_neighbour_policy(Graph& g)
        : _in_edges(num_vertices(g), edge_t()),
          _out_edges(num_vertices(g), edge_t())
    {
        for (auto v : vertices_range(g))
        {
            for (auto e : out_edges_range(v, g))
                _out_edges[v] = e;
            for (auto e : in_edges_range(v, g))
                _in_edges[v] = e;
        }
    }

    vector<edge_t> _in_edges;
    vector<edge_t> _out_edges;

    struct iter_range_t
    {
        iter_range_t(const edge_t* e): _e(e), _e_end(e)
        {
            if (_e != nullptr)
                _e_end++;
        }
        iter_t begin() {return const_cast<edge_t*>(_e);}
        iter_t end()   {return const_cast<edge_t*>(_e_end);}
        const edge_t* _e;
        const edge_t* _e_end;
    };

    template <class Vertex>
    iter_range_t
    get_out_edges(Vertex v, Graph&) const
    {
        const auto& e = _out_edges[v];
        if (e != edge_t())
            return iter_range_t(&e);
        else
            return iter_range_t(nullptr);
    }

    template <class Vertex>
    iter_range_t
    get_in_edges(Vertex v, Graph&) const
    {
        const auto& e = _in_edges[v];
        if (e != edge_t())
            return iter_range_t(&e);
        else
            return iter_range_t(nullptr);
    }

    template <class Vertex, class Weight>
    int get_out_degree(Vertex& v, Graph& g, Weight&) const
    {
        return out_degreeS()(v, g);
    }

    template <class Vertex, class Weight>
    int get_in_degree(Vertex& v, Graph& g, Weight&) const
    {
        return in_degreeS()(v, g);
    }

};

template <class Graph>
class SingleEntrySet
{
public:
    SingleEntrySet() : _pos(0) {}

    void set_move(size_t, size_t) {}

    __attribute__((always_inline))
    void insert_delta(size_t r, size_t s, int delta, bool source)
    {
        if (source)
            _entries[_pos] = make_pair(s, r);
        else
            _entries[_pos] = make_pair(r, s);
        if (!is_directed::apply<Graph>::type::value &&
            _entries[_pos].second < _entries[_pos].first)
            std::swap(_entries[_pos].first, _entries[_pos].second);
        _delta[_pos] = delta;
        ++_pos;
    }

    __attribute__((always_inline))
    int get_delta(size_t t, size_t s)
    {
        if (make_pair(t, s) == _entries[0])
            return _delta[0];
        return 0;
    }

    void clear() { _pos = 0; }

    std::array<pair<size_t, size_t>,2>& get_entries() { return _entries; }
    std::array<int, 2>& get_delta() { return _delta; }

private:
    size_t _pos;
    std::array<pair<size_t, size_t>, 2> _entries;
    std::array<int, 2> _delta;
};

template <class RNG>
size_t get_lateral_half_edge(size_t v,
                             overlap_stats_t& overlap_stats,
                             RNG& rng)
{
    size_t vv = overlap_stats.get_node(v);
    size_t w = overlap_stats.sample_half_edge(vv, rng);
    return w;
}

//A single Monte Carlo Markov chain sweep
template <class Graph, class Vprop, class VVprop, class VLprop,
          class RNG, class BlockState, class MEntries>
void move_sweep_overlap(vector<BlockState>& states,
                        vector<MEntries>& m_entries_r,
                        overlap_stats_t& overlap_stats, Vprop wr, Vprop b,
                        VLprop cv, VVprop vmap, Vprop clabel,
                        vector<int>& vlist, bool deg_corr, bool dense,
                        bool multigraph, double beta, Vprop vweight, Graph& g,
                        bool sequential, bool parallel, bool random_move,
                        double c, size_t niter, size_t B, bool verbose,
                        RNG& rng, double& S, size_t& nmoves)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

    nmoves = 0;
    S = 0;

    if (vlist.size() < 100)
        parallel = false;

    vector<pair<vertex_t, double> > best_move;
    vector<rng_t*> rngs;
    if (parallel)
    {
        best_move.resize(num_vertices(g), make_pair(vertex_t(0), numeric_limits<double>::max()));

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

    half_edge_neighbour_policy<Graph> npolicy(g);

    vector<MEntries> m_entries = m_entries_r;

    for (size_t iter = 0; iter < niter; ++iter)
    {
        std::shuffle(vlist.begin(), vlist.end(), rng);

        int i = 0, N = vlist.size();
        #pragma omp parallel for default(shared) private(i) \
            firstprivate(m_entries) schedule(runtime) if (parallel)
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
            std::uniform_int_distribution<size_t> s_rand(0, B - 1);

            vertex_t v;
            if (sequential)
            {
                v = vertex(vlist[i], g);
            }
            else
            {
                v = vertex(uniform_sample(vlist, *rngs[tid]), g);
            }

            vertex_t r = b[v];

            if (vweight[v] == 0)
                continue;

            // attempt random block
            vertex_t s = s_rand(*rngs[tid]);

            if (!random_move && total_degreeS()(v, g) > 0)
            {
                auto& state = states[0];
                size_t w = get_lateral_half_edge(v, overlap_stats, *rngs[tid]);

                int u = overlap_stats.get_out_neighbour(w);
                if (u == -1)
                    u = overlap_stats.get_in_neighbour(w);

                vertex_t t = b[u];

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

            if (s == r)
                continue;

            if (wr[s] == 0 && std::isinf(beta)) // don't populate empty blocks
                continue;

            if (clabel[s] != clabel[r])
                continue;

            if (overlap_stats.virtual_remove_size(v, r) == 0) // no empty groups
                continue;

            double dS = virtual_move(v, s, b, cv, vmap, states, m_entries,
                                     dense, deg_corr, multigraph, npolicy);

            assert (!std::isinf(dS) && !std::isnan(dS));

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

            if (overlap_stats.virtual_remove_size(v, r) == 0)
                accept = false;

            if (accept)
            {
                if (!parallel)
                {

                    assert(b[v] == int(r));
                    move_vertex(v, s, b, cv, vmap, deg_corr, states,
                                not random_move, npolicy);
                    S += dS;
                    ++nmoves;

                    assert(b[v] == int(s));
                    assert(wr[r] > 0);
                    if (verbose)
                        cout << v << ": " << r << " -> " << s << " " << S << " "
                             << vlist.size() << " " << wr[r] << " " << wr[s]
                             << " " << overlap_stats.is_enabled() << endl;
                }
                else
                {
                    best_move[v].first = s;
                    best_move[v].second = dS;
                }
            }
        }

        if (parallel)
        {
            for (vertex_t v : vlist)
            {
                vertex_t r = b[v];
                vertex_t s = best_move[v].first;
                double dS = best_move[v].second;
                if (s != r && dS != numeric_limits<double>::max())
                {
                    double dS = virtual_move(v, s, b, cv, vmap, states,
                                             m_entries, dense, deg_corr,
                                             multigraph, npolicy);
                    assert (!std::isinf(dS) && !std::isnan(dS));
                    if (dS > 0 && std::isinf(beta))
                        continue;
                    move_vertex(v, s, b, cv, vmap, deg_corr, states,
                                not random_move, npolicy);
                    S += dS;
                    ++nmoves;
                }
            }
        }
    }

    if (parallel)
    {
        for (auto r : rngs)
            delete r;
        rngs.clear();
    }
}

template <class Graph, class CEMap>
size_t get_ce(size_t v, CEMap& ce, Graph& g)
{
    for (auto e : all_edges_range(v, g))
        return ce[e];
    return 0;
}

template <class Graph, class Vprop, class VVprop, class VLprop,
          class RNG, class BlockState, class MEntries>
void coherent_move_sweep_overlap(vector<BlockState>& states,
                                 vector<MEntries>& m_entries,
                                 overlap_stats_t& overlap_stats, Vprop wr,
                                 Vprop b, VLprop cv, VVprop vmap, Vprop clabel,
                                 vector<int>& vlist, bool deg_corr, bool dense,
                                 bool multigraph, double beta, Vprop vweight,
                                 Graph& g, bool sequential, bool random_move,
                                 double c, bool confine_layers, size_t niter,
                                 size_t B, RNG& rng, double& S, size_t& nmoves)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

    typedef std::uniform_real_distribution<> rdist_t;
    auto rand_real = std::bind(rdist_t(), std::ref(rng));

    std::uniform_int_distribution<size_t> s_rand(0, B - 1);

    nmoves = 0;
    S = 0;

    half_edge_neighbour_policy<Graph> npolicy(g);

    vector<uint8_t> visited(num_vertices(g), false);

    vector<size_t> vs;

    for (size_t iter = 0; iter < niter; ++iter)
    {
        std::shuffle(vlist.begin(), vlist.end(), rng);
        std::fill(visited.begin(), visited.end(), false);

        int i = 0, N = vlist.size();
        for (i = 0; i < N; ++i)
        {
            vertex_t v;
            if (sequential)
            {
                v = vertex(vlist[i], g);
            }
            else
            {
                v = vertex(uniform_sample(vlist, rng), g);
            }

            if (visited[v])
                continue;

            vertex_t r = b[v];

            if (vweight[v] == 0)
                continue;

            // attempt random block
            vertex_t s = s_rand(rng);

            if (!random_move && total_degreeS()(v, g) > 0)
            {
                auto& state = states[0];
                size_t w = get_lateral_half_edge(v, overlap_stats, rng);

                int u = state.overlap_stats.get_out_neighbour(w);
                if (u == -1)
                    u = state.overlap_stats.get_in_neighbour(w);

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
                    const auto& e = egroups_manage::sample_edge(state.egroups[t], rng);
                    s = state.b[target(e, state.g)];
                    if (s == t)
                        s = state.b[source(e, state.g)];
                }
            }

            if (s == r)
                continue;

            if (wr[s] == 0) // don't populate empty blocks
                continue;

            if (clabel[s] != clabel[r])
                continue;


            // move all half-edges of the same node that have the same color
            size_t w = overlap_stats.get_node(v);
            int l_v = *cv[v].rbegin();
            size_t kin = 0, kout = 0;
            vs.clear();
            for (vertex_t u : overlap_stats.get_half_edges(w))
            {
                // if (visited[u])
                //     continue;

                if (b[u] != int(r))                    // skip half-edges of
                    continue;                          // different colors

                int l_u = *cv[u].rbegin();

                if (confine_layers && l_u != l_v)      // skip half-edges of
                    continue;                          // different layers

                kin += in_degreeS()(u, g);
                kout += out_degreeS()(u, g);

                vs.push_back(u);
                visited[u] = true;
            }

            if (!std::isinf(beta))
            {
                // preserve reversibility by avoiding local merges
                bool merge = false;
                for (vertex_t u : overlap_stats.get_half_edges(w))
                {
                    int l_u = *cv[u].rbegin();
                    if (size_t(b[u]) == s && (!confine_layers || l_u == l_v))
                    {
                        merge = true;
                        break;
                    }
                }
                if (merge)
                    continue;
            }

            if (overlap_stats.virtual_remove_size(v, r, kin, kout) == 0) // no empty groups
                continue;

            double dS = virtual_move(vs, s, b, cv, vmap, states, m_entries,
                                     dense, deg_corr, multigraph, npolicy);

            assert (!std::isinf(dS) && !std::isnan(dS));

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
                    double p = get_move_prob(v, r, s, c, state.b, state.mrs,
                                             state.mrp, state.mrm, state.emat,
                                             state.eweight, state.g, state.bg,
                                             m_entries[0], false,
                                             state.overlap_stats);
                    pf += p;

                    p = get_move_prob(v, s, r, c, state.b, state.mrs, state.mrp,
                                      state.mrm, state.emat, state.eweight,
                                      state.g, state.bg, m_entries[0], true,
                                      state.overlap_stats, kin, kout);
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

            if (accept && overlap_stats.get_block_size(r) > 0)
            {
                S += dS;
                nmoves += vs.size();
                move_vertex(vs, s, b, cv, vmap, deg_corr, states,
                            not random_move, npolicy);
            }
        }
    }
}


//A single Monte Carlo Markov chain sweep
template <class Graph, class Vprop, class VVprop, class VLprop, class EVprop,
          class RNG, class BlockState, class MEntries>
void merge_sweep_overlap(vector<BlockState>& states,
                         vector<MEntries>& m_entries,
                         overlap_stats_t& overlap_stats, Vprop wr, Vprop b,
                         EVprop ce, VLprop cv, VVprop vmap, Vprop clabel,
                         vector<int>& vlist, bool deg_corr, bool dense,
                         bool multigraph, Graph& g, bool random_move,
                         size_t confine_layers, size_t nmerges, size_t ntries,
                         size_t B, RNG& rng, double& S, size_t& nmoves)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

    Vprop last_b = b.copy();
    Vprop best_move = b.copy();
    vector<double> best_move_dS;
    vector<vector<vertex_t>> groups(B);
    vector<unordered_map<size_t, vector<size_t>>> bundles(B);
    best_move_dS.resize(B, numeric_limits<double>::max());

    std::shuffle(vlist.begin(), vlist.end(), rng);

    vector<vector<vertex_t>> blocks;
    for (size_t r = 0; r < B; ++r)
    {
        if (wr[r] > 0)
        {
            size_t x = clabel[r];
            if (x >= blocks.size())
                blocks.resize(x + 1);
            blocks[x].push_back(r);
        }
    }

    for (auto v : vlist)
    {
        groups[b[v]].push_back(v);

        if (!confine_layers)
        {
            //separate into subgroups according to opposing sides
            int u = overlap_stats.get_out_neighbour(v);
            if (u == -1)
                u = overlap_stats.get_in_neighbour(v);
            bundles[b[v]][b[u]].push_back(v);
        }
        else
        {
            //separate into subgroups according to edge layer
            bundles[b[v]][get_ce(v, ce, g)].push_back(v);
        }
    }

    std::uniform_int_distribution<size_t> s_rand(0, B - 1);

    nmoves = 0;
    S = 0;

    half_edge_neighbour_policy<Graph> npolicy(g);

    for (size_t r = 0; r < B; ++r)
    {
        if (groups[r].empty())
            continue;

        double dS = 0;

        for (size_t j = 0; j < ntries; ++j)
        {
            for (auto& tb : bundles[r])
            {
                auto& bundle = tb.second;

                vertex_t s;

                if (!random_move)
                {
                    vertex_t v = uniform_sample(bundle, rng);
                    size_t w  = get_lateral_half_edge(v, overlap_stats, rng);
                    auto& state = states[0];

                    // neighbour sampler points to the *block graph*
                    vertex_t t = last_b[w];
                    s = state.neighbour_sampler[t].sample(rng);
                    if (s == r)
                        s = state.cavity_neighbour_sampler[s].sample(rng);
                    else
                        s = state.neighbour_sampler[s].sample(rng);
                }
                else
                {
                    // attempt random block
                    s = uniform_sample(blocks[clabel[r]], rng);
                }

                if (s == r)
                    continue;

                if (wr[s] == 0) // don't populate empty blocks
                    continue;

                if (clabel[s] != clabel[r])
                    continue;

                double ddS = 0;
                for (auto u : bundle)
                {
                    ddS += virtual_move(u, s, b, cv, vmap, states, m_entries,
                                        dense, deg_corr, multigraph, npolicy);
                    move_vertex(u, s, b, cv, vmap, deg_corr, states, false,
                                npolicy);
                }


                if (ddS < 0 || best_move[bundle[0]] == int(r))
                {
                    for (auto u : bundle)
                        best_move[u] = s;
                    dS += ddS;
                }
                else
                {
                    for (auto u : bundle)
                    {
                        move_vertex(u, best_move[u], b, cv, vmap, deg_corr,
                                    states, false, npolicy);
                    }
                }

            }
        }

        // reset
        for (auto v : groups[r])
        {
            move_vertex(v, r, b, cv, vmap, deg_corr, states, false, npolicy);
        }

        best_move_dS[r] = dS;
    }

    auto cmp = [&](vertex_t r, vertex_t s) -> bool {return best_move_dS[r] > best_move_dS[s];};

    // top is the merge with _smallest_ dS
    priority_queue<vertex_t, vector<vertex_t>, decltype(cmp)> move_heap(cmp);

    for (size_t r = 0; r < B; ++r)
    {
        if (best_move_dS[r] == numeric_limits<double>::max() || groups[r].empty())
            continue;
        move_heap.push(r);
    }

    vector<bool> touched(B, false);
    while (!move_heap.empty() && nmoves < nmerges)
    {
        vertex_t r = move_heap.top();
        move_heap.pop();

        if (touched[r])
            continue;

        bool skip = false;
        for (vertex_t v : groups[r])
        {
            vertex_t s = best_move[v];
            if (wr[s] == 0)
            {
                skip = true;
                break;
            }
        }

        if (skip)
            continue;

        touched[r] = true;

        for (auto v : groups[r])
        {
            vertex_t s = best_move[v];
            touched[s] = true;

            //assert(s != r);
            double dS = virtual_move(v, s, b, cv, vmap, states, m_entries,
                                     dense, deg_corr, multigraph, npolicy);
            move_vertex(v, s, b, cv, vmap, deg_corr, states, false, npolicy);
            S += dS;
        }

        //assert (wr[r] == 0);
        if (wr[r] == 0)
            ++nmoves;
    }
}

} // namespace graph_tool

#endif // GRAPH_BLOCKMODEL_OVERLAP_HH
