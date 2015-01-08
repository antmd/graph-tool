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
    typedef property_map_type::apply<vector<int32_t>,
                                     GraphInterface::vertex_index_map_t>::type::unchecked_t
        vvmap_t;

    overlap_stats_t(): _enabled(false) {}

    template <class Graph>
    overlap_stats_t(Graph& g, vmap_t b, vvmap_t half_edges, vmap_t node_index,
                    size_t B)
        : _half_edges(half_edges), _node_index(node_index),
          _out_neighbours(num_vertices(g), -1), _in_neighbours(num_vertices(g), -1)
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

        for (auto v : vertices_range(g))
        {
            size_t vi = node_index[v];
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
        vector<bool> visited(num_vertices(g), false);
        for (auto v : vertices_range(g))
        {
            if (visited[v])
                continue;

            size_t i = node_index[v];
            auto& he = half_edges[i];

            unordered_map<size_t, vector<size_t>> out_us;
            for (auto u : he)
            {
                auto w = _out_neighbours[u];
                if (w == -1)
                    continue;
                out_us[node_index[w]].push_back(u);
                visited[u] = true;
            }

            for (auto& uc : out_us)
            {
                if (uc.second.size() > 1)
                {
                    for (auto u : uc.second)
                        _mi[u] = _parallel_bundles.size();
                    _parallel_bundles.push_back(uc.second);
                }
            }
        }
    }

    template <class Graph>
    void add_half_edge(size_t v, size_t r, Graph& g)
    {
        size_t u = _node_index[v];
        size_t kin = in_degreeS()(v, g);
        size_t kout = out_degreeS()(v, g);
        auto& k = _block_nodes[r][u];
        k.first += kin;
        k.second += kout;
    }

    template <class Graph>
    void remove_half_edge(size_t v, size_t r, Graph& g)
    {
        size_t u = _node_index[v];
        size_t kin = in_degreeS()(v, g);
        size_t kout = out_degreeS()(v, g);

        auto& k = _block_nodes[r][u];
        k.first -= kin;
        k.second -= kout;

        if (k.first + k.second == 0)
            _block_nodes[r].erase(u);
    }

    template <class Graph>
    void move_half_edge(size_t v, size_t r, size_t nr, Graph& g)
    {
        remove_half_edge(v, r, g);
        add_half_edge(v, nr, g);
    }

    size_t get_block_size(size_t r)
    {
        return _block_nodes[r].size();
    }

    template <class Graph>
    size_t virtual_remove_size(size_t v, size_t r, Graph& g)
    {
        size_t nr = _block_nodes[r].size();
        size_t u = _node_index[v];
        size_t kin = in_degreeS()(v, g);
        size_t kout = out_degreeS()(v, g);
        auto& deg = _block_nodes[r][u];
        if (deg.first == kin && deg.second == kout)
            nr--;
        return nr;
    }

    template <class Graph>
    size_t virtual_add_size(size_t v, size_t r, Graph& g)
    {
        size_t nr = _block_nodes[r].size();
        size_t u = _node_index[v];
        auto& bnodes = _block_nodes[r];
        if (bnodes.find(u) == bnodes.end())
            nr++;
        return nr;
    }

    template <class Graph>
    double virtual_move_dS(size_t v, size_t r, size_t nr, Graph& g)
    {
        double dS = 0;

        size_t u = _node_index[v];
        size_t u_kin = in_degreeS()(v, g);
        size_t u_kout = out_degreeS()(v, g);

        auto deg = _block_nodes[r][u];
        auto ndeg = deg;
        ndeg.first -= u_kin;
        ndeg.second -= u_kout;

        dS -= lgamma_fast(ndeg.first + 1) + lgamma_fast(ndeg.second + 1);
        dS += lgamma_fast(deg.first + 1) + lgamma_fast(deg.second + 1);

        auto iter = _block_nodes[nr].find(u);
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
    double virtual_move_parallel_dS(size_t v, size_t r, size_t nr, VProp& b, Graph& g)
    {
        int m = _mi[v];
        if (_out_neighbours[v] == -1)
            m = _mi[_in_neighbours[v]];
        if (m == -1)
            return 0;

        double S = 0;
        typedef std::tuple<size_t, size_t, size_t> key_t; // u, r, s

#ifdef HAVE_SPARSEHASH
        dense_hash_map<key_t, int, std::hash<key_t>> out_us;
        out_us.set_empty_key(key_t(numeric_limits<size_t>::max(),
                                   numeric_limits<size_t>::max(),
                                   numeric_limits<size_t>::max()));
#else
        unordered_map<key_t, int> out_us;
#endif

        for (auto u : _parallel_bundles[m])
        {
            size_t r_u = b[u];

            auto w = _out_neighbours[u];
            assert(w != -1);

            size_t s = b[w];
            size_t ni_w = _node_index[w];
            size_t ni_u = _node_index[u];

            if (ni_w == ni_u && !is_directed::apply<Graph>::type::value && s < r_u)
                out_us[key_t(ni_w, s, r_u)]++;
            else
                out_us[key_t(ni_w, r_u, s)]++;
        }

        size_t vi = _node_index[v];

        auto get_node = [](const key_t& k) -> size_t {return get<0>(k);};
        S -= get_parallel_neighbours_entropy(vi, out_us, g, get_node);

        out_us.clear();
        for (auto u : _parallel_bundles[m])
        {
            size_t r_u = b[u];
            if (u == v)
                r_u = nr;

            auto w = _out_neighbours[u];
            assert(w != -1);

            size_t s = b[w];

            if (w == int(v))
                s = nr;

            size_t ni_w = _node_index[w];
            size_t ni_u = _node_index[u];
            if (ni_w == ni_u && !is_directed::apply<Graph>::type::value && s < r_u)
                out_us[key_t(ni_w, s, r_u)]++;
            else
                out_us[key_t(ni_w, r_u, s)]++;
        }

        S += get_parallel_neighbours_entropy(vi, out_us, g, get_node);

        return S;
    }

    // sample another half-edge adjacent to the same node of half-edge w
    template <class RNG>
    size_t sample_half_edge(size_t w, RNG& rng)
    {
        size_t wi = _node_index[w];
        auto& half_edges = _half_edges[wi];
        return uniform_sample(half_edges, rng);
    }

    size_t get_node(size_t v) { return _node_index[v]; }
    vector<int>& get_half_edges(size_t v) { return _half_edges[v]; }

    int get_out_neighbour(size_t v) { return _out_neighbours[v]; }
    int get_in_neighbour(size_t v) { return _in_neighbours[v]; }


    bool is_enabled() { return _enabled; }

private:
    bool _enabled;
    vvmap_t _half_edges;     // half-edges to each node
    vmap_t  _node_index;      // node to each half edges

#ifdef HAVE_SPARSEHASH
        typedef dense_hash_map<size_t, deg_t> node_map_t;
#else
        typedef unordered_map<size_t, deg_t> node_map_t;
#endif

    vector<node_map_t> _block_nodes; // nodes (and degrees) in each block

    vector<int> _out_neighbours;
    vector<int> _in_neighbours;

    vector<int> _mi;
    vector<vector<size_t>> _parallel_bundles; // parallel edge bundles
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
                              Eprop eweight, size_t N, size_t B)
        : _enabled(true), _N(N), _B(B), _D(0),
          _dhist(B + 1), _r_count(B), _bhist(N), _emhist(B), _ephist(B),
          _embhist(N), _epbhist(N), _deg_hist(N), _bvs(N), _nbvs(N), _degs(N), _ndegs(N)
    {
#ifdef HAVE_SPARSEHASH
        _bhist.set_empty_key(bv_t());
        _deg_hist.set_empty_key(bv_t());
        _embhist.set_empty_key(bv_t());
        _epbhist.set_empty_key(bv_t());
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
                cdh.set_empty_key(cdeg_t());
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
                _r_count[r]++;
        }

    }

    template <class Graph, class Vprop, class Eprop>
    void get_bv_deg(size_t v, Vprop& b, Eprop& eweight,
                    overlap_stats_t& overlap_stats, Graph& g, set<size_t>& rs,
                    dmap_t& in_hist, dmap_t& out_hist)
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


    double get_partition_dl()
    {
        double S = 0;
        for (size_t d = 1; d < _dhist.size(); ++d)
        {
            size_t nd = _dhist[d];
            if (nd == 0)
                continue;
            double x = lbinom_fast(_B, d);
            double ss = lbinom((exp(x) + nd) - 1, nd); // no fast
            if (std::isinf(ss) || std::isnan(ss))
                ss = nd * x - lgamma_fast(nd + 1);
            assert(!std::isinf(ss));
            assert(!std::isnan(ss));
            S += ss;
        }

        S += lbinom(_D + _N - 1, _N) + lgamma(_N + 1);

        for (auto& bh : _bhist)
            S -= lgamma(bh.second + 1);

        // double S1 = S;
        // S = 0;
        // vector<int> rhist(_B, 0);
        // for (auto& bh : _bhist)
        //     for (auto r : bh.first)
        //         rhist[r] += bh.second;

        // for (size_t r = 0; r < _B; ++r)
        // {
        //     S += lbinom(2 + _N - 1, _N);
        //     S += lgamma(_N + 1) - lgamma(rhist[r] + 1) - lgamma(_N - rhist[r] + 1);
        // }

        // double S2 = S;

        // // cout << "P:" << S1 << " " << S2 << " " << S1 - S2 << endl;

        // S = min(S1, S2);

        return S;
    }

    double get_deg_dl(bool ent, bool deg_alt, bool xi_fast)
    {
        double S = 0;
        if (ent)
        {
            for (auto& ch : _deg_hist)
            {
                auto& bv = ch.first;
                auto& cdeg_hist = ch.second;

                size_t n_bv = _bhist[bv];

                S += xlogx(n_bv);
                for (auto& dh : cdeg_hist)
                    S -= xlogx(dh.second);
            }

            // unordered_map<size_t,
            //               unordered_map<deg_t, size_t>> bdeg_hist;

            // for (auto& bv_h : _deg_hist)
            // {
            //     auto& bv = bv_h.first;
            //     auto& cdeg_hist = bv_h.second;
            //     size_t total = 0;
            //     for (auto& d_c : cdeg_hist)
            //     {
            //         assert(d_c.first.size() == bv.size());
            //         for (size_t i = 0; i < bv.size(); ++i)
            //         {
            //             size_t r = bv[i];
            //             bdeg_hist[r][d_c.first[i]] += d_c.second;
            //         }
            //         total += d_c.second;
            //     }
            // }

            // for (auto& r_h : bdeg_hist)
            // {
            //     //auto r = r_h.first;
            //     auto& hist = r_h.second;

            //     size_t total = 0;
            //     for (auto& d_c : hist)
            //         total += d_c.second;

            //     for (auto& d_c : hist)
            //     {
            //         double c = d_c.second;
            //         S -= c * log(c / total);
            //     }
            // }

        }
        else
        {
            S = 0;

            for (auto& ch : _deg_hist)
            {
                auto& bv = ch.first;
                auto& cdeg_hist = ch.second;

                size_t n_bv = _bhist[bv];

                if (n_bv == 0)
                    continue;

                auto& bmh = _embhist[bv];
                auto& bph = _epbhist[bv];

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

                S1 += lgamma(n_bv + 1);

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

            //const double z2 = boost::math::zeta(2);

            // unordered_map<size_t,
            //               unordered_map<deg_t, size_t>>> bdeg_hist;

            // for (auto& bv_h : _deg_hist)
            // {
            //     auto& bv = bv_h.first;
            //     auto& cdeg_hist = bv_h.second;
            //     size_t total = 0;
            //     for (auto& d_c : cdeg_hist)
            //     {
            //         assert(d_c.first.size() == bv.size());
            //         for (size_t i = 0; i < bv.size(); ++i)
            //         {
            //             size_t r = bv[i];
            //             bdeg_hist[r][d_c.first[i]] += d_c.second;
            //         }
            //         total += d_c.second;
            //     }
            // }

            // for (auto& r_h : bdeg_hist)
            // {
            //     auto r = r_h.first;
            //     auto& hist = r_h.second;

            //     // S += 2 * sqrt(z2 * _emhist[r]);
            //     // S += 2 * sqrt(z2 * _ephist[r]);

            //     size_t total = 0;
            //     for (auto& d_c : hist)
            //     {
            //         total += d_c.second;
            //         S -= lgamma(d_c.second + 1);
            //     }
            //     S += lgamma(total + 1);

            //     S += get_xi(total, _emhist[r]);
            //     S += get_xi(total, _ephist[r]);

            //     // cout << total << " " <<  _ephist[r] << " " << 2 * sqrt(z2 * _ephist[r]) << " " << get_xi(total,_ephist[r]) << endl;
            //     // cout << total << " " <<  _emhist[r] << " " << 2 * sqrt(z2 * _emhist[r]) << " " << get_xi(total,_emhist[r]) << endl;

            // }

            // double S1 = S;

            // double S2 = S;

            // // cout << S1 << " " << S2 << " " << S1 - S2 << endl;
            // S = min(S1, S2);

        }
        return S;
    }


    template <class Graph>
    bool get_n_bv(const vector<size_t>& vs, const vector<size_t>& rs,
                  const vector<size_t>& nrs, const bv_t& bv, const cdeg_t& deg,
                  bv_t& n_bv, cdeg_t& n_deg, Graph& g)
    {
#ifdef HAVE_SPARSEHASH
        dense_hash_map<size_t, pair<int, int>> deg_delta(deg.size());
        deg_delta.set_empty_key(numeric_limits<size_t>::max());
        deg_delta.set_deleted_key(numeric_limits<size_t>::max() - 1);
#else
        unordered_map<size_t, pair<int, int>> deg_delta(deg.size());
#endif

        for (size_t i = 0; i < vs.size(); ++i)
        {
            auto v = vs[i];
            auto r = rs[i];
            auto nr = nrs[i];
            if (r == nr)
                continue;
            size_t kin = in_degreeS()(v, g);
            size_t kout = out_degreeS()(v, g);

            auto& d_r = deg_delta[r];
            d_r.first -= kin;
            d_r.second -= kout;

            auto& d_nr = deg_delta[nr];
            d_nr.first += kin;
            d_nr.second += kout;
        }

        n_deg.clear();
        n_bv.clear();
        bool is_same_bv = true;
        for (size_t i = 0; i < bv.size(); ++i)
        {
            size_t s = bv[i];
            auto k_s = deg[i];

            auto iter = deg_delta.find(s);
            if (iter != deg_delta.end())
            {
                get<0>(k_s) += iter->second.first;
                get<1>(k_s) += iter->second.second;
                deg_delta.erase(iter);
            }

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

        if (!deg_delta.empty())
        {
            is_same_bv = false;
            for (auto& dd : deg_delta)
            {
                auto nr = dd.first;
                assert(dd.second.first + dd.second.second > 0);
                size_t kin = dd.second.first;
                size_t kout = dd.second.second;
                auto pos = std::lower_bound(n_bv.begin(), n_bv.end(), nr);
                auto dpos = n_deg.begin();
                std::advance(dpos, pos - n_bv.begin());
                n_bv.insert(pos, nr);
                n_deg.insert(dpos, make_pair(kin, kout));
            }
        }
        return is_same_bv;
    }

    // get deg counts without increasing the container
    size_t get_deg_count(bv_t& bv, cdeg_t& deg)
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
    size_t get_bv_count(bv_t& bv)
    {
        auto iter = _bhist.find(bv);
        if (iter == _bhist.end())
            return 0;
        return iter->second;
    }

    template <class Graph>
    double get_delta_dl(size_t v, size_t r, size_t nr, bool deg_corr,
                        overlap_stats_t& overlap_stats, Graph& g)
    {
        const vector<size_t> vs = {v}, rs = {r}, nrs = {nr};
        return get_delta_dl(vs, rs, nrs, deg_corr, overlap_stats, g);
    }

    template <class Graph>
    double get_delta_dl(const vector<size_t>& vs, const vector<size_t>& rs,
                        const vector<size_t>& nrs, bool deg_corr,
                        overlap_stats_t& overlap_stats, Graph& g)
    {
        if (rs == nrs)
            return 0;

        size_t u = overlap_stats.get_node(vs[0]);
        auto& bv = _bvs[u];
        auto& n_bv = _nbvs[u];
        size_t d = bv.size();
        cdeg_t& deg = _degs[u];
        cdeg_t& n_deg = _ndegs[u];


        bool is_same_bv = get_n_bv(vs, rs, nrs, bv, deg, n_bv, n_deg, g);

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

        double S = 0;
        if (n_D != _D)
        {
            S -= lbinom(_D + _N - 1, _N);
            S += lbinom(n_D + _N - 1, _N);
        }

        double S_a = 0, S_b = 0;

        auto get_S_d = [&] (size_t d_i, int delta) -> double
            {
                int nd = int(_dhist[d_i]) + delta;
                if (nd == 0)
                    return 0.;
                double x = lbinom_fast(_B, d_i);
                double S = lbinom(exp(x) + nd - 1, nd); // no fast
                if (std::isinf(S) || std::isnan(S))
                    S = nd * x - lgamma_fast(nd + 1);
                return S;
            };

        if (n_d != d)
        {
            S_b += get_S_d(d,  0) + get_S_d(n_d, 0);
            S_a += get_S_d(d, -1) + get_S_d(n_d, 1);
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

        S += S_a - S_b;

        if (deg_corr)
        {
            double S_b = 0, S_a = 0;

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

                    cdeg_t& deg_i = (is_bv) ? deg : n_deg;
                    auto& bv_i = (is_bv) ? bv : n_bv;

                    double S = 0;
                    if (((is_bv) ? bv_count : n_bv_count) > 0)
                    {
                        auto& bmh = _embhist[bv_i];
                        auto& bph = _epbhist[bv_i];

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
                    auto& bmh = _embhist[bv];
                    auto& bph = _epbhist[bv];

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

            auto is_in = [&] (bv_t& bv, size_t r) -> bool
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
            dense_hash_map<size_t, pair<int, int>> deg_delta;
            dense_hash_map<size_t, int> r_count_delta;
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

            for (size_t i = 0; i < vs.size(); ++i)
            {
                auto v = vs[i];
                auto r = rs[i];
                auto nr = nrs[i];
                if (r == nr)
                    continue;
                size_t kin = in_degreeS()(v, g);
                size_t kout = out_degreeS()(v, g);

                auto& d_r = deg_delta[r];
                d_r.first -= kin;
                d_r.second -= kout;

                auto& d_nr = deg_delta[nr];
                d_nr.first += kin;
                d_nr.second += kout;
            }

            for (size_t s : bv)
            {
                S_a += lbinom_fast(_r_count[s] + r_count_delta[s] + _emhist[s] + deg_delta[s].first - 1, _emhist[s] + deg_delta[s].first);
                S_a += lbinom_fast(_r_count[s] + r_count_delta[s] + _ephist[s] + deg_delta[s].second - 1, _ephist[s] + deg_delta[s].second);
            }

            for (size_t s : n_bv)
            {
                if (!is_in(bv, s))
                {
                    S_a += lbinom_fast(_r_count[s] + r_count_delta[s] + _emhist[s] + deg_delta[s].first - 1, _emhist[s] + deg_delta[s].first);
                    S_a += lbinom_fast(_r_count[s] + r_count_delta[s] + _ephist[s] + deg_delta[s].second - 1, _ephist[s] + deg_delta[s].second);
                }
            }

            S += S_a - S_b;
        }

        return S;
    }

    template <class Graph>
    void move_vertex(size_t v, size_t r, size_t nr, bool deg_corr,
                     overlap_stats_t& overlap_stats, Graph& g)
    {
        const vector<size_t> vs = {v}, rs = {r}, nrs = {nr};
        move_vertex(vs, rs, nrs, deg_corr, overlap_stats, g);

    }

    template <class Graph>
    void move_vertex(const vector<size_t>& vs, const vector<size_t>& rs,
                     const vector<size_t>& nrs, bool deg_corr,
                     overlap_stats_t& overlap_stats, Graph& g)
    {
        if (rs == nrs)
            return;

        size_t u = overlap_stats.get_node(vs[0]);
        auto& bv = _bvs[u];
        auto& n_bv = _nbvs[u];
        cdeg_t& deg = _degs[u];
        cdeg_t& n_deg = _ndegs[u];
        size_t d = bv.size();

        bool is_same_bv = get_n_bv(vs, rs, nrs, bv, deg, n_bv, n_deg, g);
        size_t n_d = n_bv.size();

        if (!is_same_bv)
        {
            _dhist[d] -= 1;
            auto& bv_count = _bhist[bv];
            bv_count -= 1;

            if (bv_count == 0)
            {
                for (auto s : bv)
                    _r_count[s]--;
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
                    _r_count[s]++;
            }

            if (n_d > _D)
                _D = n_d;
        }

        _deg_hist[bv][deg] -= 1;
        auto& bmh = _embhist[bv];
        auto& bph = _epbhist[bv];
        assert(bmh.size() == bv.size());
        assert(bph.size() == bv.size());
        for (size_t i = 0; i < bv.size(); ++i)
        {
            bmh[i] -= get<0>(deg[i]);
            bph[i] -= get<1>(deg[i]);
        }

        for (size_t i = 0; i < vs.size(); ++i)
        {
            auto v = vs[i];
            auto r = rs[i];
            size_t kin = in_degreeS()(v, g);
            size_t kout = out_degreeS()(v, g);
            _emhist[r] -= kin;
            _ephist[r] -= kout;
        }

        auto& hist = _deg_hist[n_bv];
#ifdef HAVE_SPARSEHASH
        if (hist.empty())
            hist.set_empty_key(cdeg_t());
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

        for (size_t i = 0; i < vs.size(); ++i)
        {
            auto v = vs[i];
            auto nr = nrs[i];
            size_t kin = in_degreeS()(v, g);
            size_t kout = out_degreeS()(v, g);
            _emhist[nr] += kin;
            _ephist[nr] += kout;
        }

        _bvs[u].swap(n_bv);
        _degs[u].swap(n_deg);
    }

    bool is_enabled() { return _enabled; }

private:
    bool _enabled;
public:
    size_t _N;
private:
    size_t _B;
    size_t _D;
    vector<int> _dhist; // d-histogram
    vector<int> _r_count; // m_r
    bhist_t _bhist;     // b-histogram
    vector<size_t> _emhist;    // e-_r histogram
    vector<size_t> _ephist;    // e+_r histogram
    ebhist_t _embhist;  // e+^r_b histogram
    ebhist_t _epbhist;  // e+^r_b histogram
    deg_hist_t _deg_hist; // n_k^b histogram
    vector<bv_t> _bvs;    // bv node map
    vector<bv_t> _nbvs;
    vector<cdeg_t> _degs;  // deg node map
    vector<cdeg_t> _ndegs;

};

struct entropy_parallel_edges_overlap
{
    template <class Graph, class VProp>
    void operator()(Graph& g, VProp b, overlap_stats_t& overlap_stats, double& S) const
    {
        S = 0;
        vector<bool> visited(num_vertices(g), false);
        for (auto v : vertices_range(g))
        {
            if (visited[v])
                continue;

            size_t i = overlap_stats.get_node(v);
            auto& he = overlap_stats.get_half_edges(i);

            typedef std::tuple<size_t, size_t, size_t> key_t; // u, r, s
            unordered_map<key_t, int> out_us;
            for (auto u : he)
            {
                if (visited[u])
                    continue;
                visited[u] = true;

                size_t r = b[u];
                int w = overlap_stats.get_out_neighbour(u);
                if (w == -1)
                    continue;

                if (!is_directed::apply<Graph>::type::value)
                    visited[w] = true;

                size_t ni_w = overlap_stats.get_node(w);
                size_t ni_u = overlap_stats.get_node(u);

                size_t s = b[w];

                if (ni_w == ni_u && !is_directed::apply<Graph>::type::value && s < r)
                    out_us[key_t(ni_w, s, r)]++;
                else
                    out_us[key_t(ni_w, r, s)]++;
            }

            auto get_node = [](const key_t& k) -> size_t {return get<0>(k);};
            S += get_parallel_neighbours_entropy(i, out_us, g, get_node);
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
    get_out_edges(Vertex v, Graph& g) const
    {
        const auto& e = _out_edges[v];
        if (e != edge_t())
            return iter_range_t(&e);
        else
            return iter_range_t(nullptr);
    }

    template <class Vertex>
    iter_range_t
    get_in_edges(Vertex v, Graph& g) const
    {
        const auto& e = _in_edges[v];
        if (e != edge_t())
            return iter_range_t(&e);
        else
            return iter_range_t(nullptr);
    }

    template <class Vertex, class Weight>
    int get_out_degree(Vertex& v, Graph& g, Weight& eweight) const
    {
        return out_degreeS()(v, g);
    }

    template <class Vertex, class Weight>
    int get_in_degree(Vertex& v, Graph& g, Weight& eweight) const
    {
        return in_degreeS()(v, g);
    }

};


//A single Monte Carlo Markov chain sweep
template <class Graph, class BGraph, class EMprop, class Eprop, class Vprop,
          class EMat, class EVprop, class VEprop, class RNG>
void move_sweep_overlap(EMprop mrs, Vprop mrp, Vprop mrm, Vprop wr, Vprop b,
                        Vprop clabel, vector<int>& vlist, bool deg_corr,
                        bool dense, bool multigraph, bool parallel_edges,
                        double beta, Eprop eweight, Vprop vweight,
                        EVprop egroups, VEprop esrcpos, VEprop etgtpos,
                        Graph& g, BGraph& bg, EMat& emat, bool sequential,
                        bool parallel, bool random_move, double c,
                        overlap_stats_t& overlap_stats,
                        overlap_partition_stats_t& partition_stats,
                        bool verbose, RNG& rng, double& S, size_t& nmoves)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

    size_t B = num_vertices(bg);

    std::shuffle(vlist.begin(), vlist.end(), rng);

    nmoves = 0;
    S = 0;

    if (vlist.size() < 100)
        parallel = false;

    EntrySet<Graph> m_entries(B);

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
            // attempt "lateral" moves across overlaps
            vertex_t w = overlap_stats.sample_half_edge(v, *rngs[tid]);

            int u = overlap_stats.get_out_neighbour(w);
            if (u == -1)
                u = overlap_stats.get_in_neighbour(w);

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

        if (s == r)
            continue;

        if (wr[s] == 0 && std::isinf(beta)) // don't populate empty blocks
            continue;

        if (clabel[s] != clabel[r])
            continue;

        double dS = virtual_move(v, s, dense, mrs, mrp, mrm, wr, b, deg_corr,
                                 eweight, vweight, g, bg, emat, m_entries,
                                 overlap_stats, multigraph, partition_stats,
                                 npolicy);

        bool accept = false;
        if (std::isinf(beta))
        {
            accept = dS < 0;
        }
        else
        {
            double pf = random_move ? 1 :
                get_move_prob(v, r, s, c, b, mrs, mrp, mrm, emat, eweight, g,
                              bg, m_entries, false, overlap_stats);

            double pb = random_move ? 1 :
                get_move_prob(v, s, r, c, b, mrs, mrp, mrm, emat, eweight, g,
                              bg, m_entries, true, overlap_stats);

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

        if (overlap_stats.virtual_remove_size(v, r, g) == 0)
            accept = false;

        if (accept)
        {
            if (!parallel)
            {

                assert(b[v] == int(r));
                move_vertex(v, s, mrs, mrp, mrm, wr, b, deg_corr, eweight,
                            vweight, g, bg, emat, overlap_stats,
                            partition_stats, npolicy);

                //partition_stats = overlap_partition_stats_t(g, b, overlap_stats, eweight, partition_stats._N, B);

                if (!random_move)
                    egroups_manage::update_egroups(v, r, s, eweight, egroups,
                                                   esrcpos, etgtpos, g);

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
                                         multigraph, partition_stats, npolicy);
                if (dS > 0 && std::isinf(beta))
                    continue;
                move_vertex(v, s, mrs, mrp, mrm, wr, b, deg_corr, eweight,
                            vweight, g, bg, emat, overlap_stats,
                            partition_stats, npolicy);
                //partition_stats = overlap_partition_stats_t(g, b, overlap_stats, eweight, partition_stats._N, B);

                if (!random_move)
                    egroups_manage::update_egroups(v, r, s, eweight, egroups,
                                                   esrcpos, etgtpos, g);
                S += dS;
                ++nmoves;
            }
        }
    }

}


template <class Graph, class BGraph, class EMprop, class Eprop, class Vprop,
          class EMat, class EVprop, class VEprop, class RNG>
void coherent_move_sweep_overlap(EMprop mrs, Vprop mrp, Vprop mrm, Vprop wr,
                                 Vprop b, Vprop clabel, vector<int>& vlist,
                                 bool deg_corr, bool dense, bool multigraph,
                                 bool parallel_edges, Eprop eweight,
                                 Vprop vweight, EVprop egroups, VEprop esrcpos,
                                 VEprop etgtpos, Graph& g, BGraph& bg,
                                 EMat& emat, bool sequential, bool parallel,
                                 bool random_move, double c,
                                 overlap_stats_t& overlap_stats,
                                 overlap_partition_stats_t& partition_stats,
                                 bool verbose, RNG& rng, double& S,
                                 size_t& nmoves)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

    size_t B = num_vertices(bg);

    typedef std::uniform_real_distribution<> rdist_t;
    auto rand_real = std::bind(rdist_t(), std::ref(rng));

    std::shuffle(vlist.begin(), vlist.end(), rng);
    std::uniform_int_distribution<size_t> s_rand(0, B - 1);

    nmoves = 0;
    S = 0;

    EntrySet<Graph> m_entries(B);

    half_edge_neighbour_policy<Graph> npolicy(g);
    overlap_partition_stats_t dummy_pstats;

    vector<bool> visited(num_vertices(g), false);

    Vprop last_b = b.copy();

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

        assert(last_b[v] == int(r));

        if (vweight[v] == 0)
            continue;

        // attempt random block
        vertex_t s = s_rand(rng);

        if (!random_move && total_degreeS()(v, g) > 0)
        {
            // attempt "lateral" moves across overlaps
            vertex_t w = overlap_stats.sample_half_edge(v, rng);

            int u = overlap_stats.get_out_neighbour(w);
            if (u == -1)
                u = overlap_stats.get_in_neighbour(w);

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
                const auto& e = egroups_manage::sample_edge(egroups[t], rng);
                s = b[target(e, g)];
                if (s == t)
                    s = b[source(e, g)];
            }
        }

        if (s == r)
            continue;

        if (wr[s] == 0) // don't populate empty blocks
            continue;

        if (clabel[s] != clabel[r])
            continue;

        // move all half-edges of the same node with the same color

        double dS = 0;
        size_t w = overlap_stats.get_node(v);
        size_t nm = 0;
        vector<size_t> vs, rs, nrs;

        if (partition_stats.is_enabled())
        {
            for (vertex_t u : overlap_stats.get_half_edges(w))
            {
                if (last_b[u] != int(r))
                    continue;
                vs.push_back(u);
                rs.push_back(r);
                nrs.push_back(s);
            }

            dS += partition_stats.get_delta_dl(vs, rs, nrs, deg_corr,
                                               overlap_stats, g);
            partition_stats.move_vertex(vs, rs, nrs, deg_corr, overlap_stats, g);
        }

        for (vertex_t u : overlap_stats.get_half_edges(w))
        {
            if (last_b[u] != int(r))
                continue;
            assert(last_b[u] == int(r));
            assert(!visited[u]);
            visited[u] = true;

            dS += virtual_move(u, s, dense, mrs, mrp, mrm, wr, b,
                               deg_corr, eweight, vweight, g, bg, emat,
                               m_entries, overlap_stats, multigraph,
                               dummy_pstats, npolicy);

            move_vertex(u, s, mrs, mrp, mrm, wr, b, deg_corr, eweight, vweight,
                        g, bg, emat, overlap_stats, dummy_pstats, npolicy);
            //partition_stats = overlap_partition_stats_t(g, b, overlap_stats, eweight, partition_stats._N, B);

            if (!random_move)
                egroups_manage::update_egroups(u, r, s, eweight, egroups,
                                               esrcpos, etgtpos, g);
            nm++;

        }


        if (dS < 0 && overlap_stats.get_block_size(r) > 0)
        {
            S += dS;
            nmoves += nm;
        }
        else
        {
            if (partition_stats.is_enabled())
            {
                dS += partition_stats.get_delta_dl(vs, nrs, rs, deg_corr, overlap_stats, g);
                partition_stats.move_vertex(vs, nrs, rs, deg_corr, overlap_stats, g);
            }

            for (vertex_t u : overlap_stats.get_half_edges(w))
            {
                if (last_b[u] != int(r))
                    continue;
                move_vertex(u, r, mrs, mrp, mrm, wr, b, deg_corr, eweight, vweight,
                            g, bg, emat, overlap_stats, dummy_pstats, npolicy);
                //partition_stats = overlap_partition_stats_t(g, b, overlap_stats, eweight, partition_stats._N, B);
                if (!random_move)
                    egroups_manage::update_egroups(u, s, r, eweight, egroups,
                                                   esrcpos, etgtpos, g);
            }
        }
    }
}


//A single Monte Carlo Markov chain sweep
template <class Graph, class BGraph, class EMprop, class Eprop, class Vprop,
          class EMat, class EVprop, class VEprop, class SamplerMap, class RNG>
void merge_sweep_overlap(EMprop mrs, Vprop mrp, Vprop mrm, Vprop wr, Vprop b,
                         Vprop clabel, vector<int>& vlist, bool deg_corr,
                         bool dense, bool multigraph, bool parallel_edges,
                         Eprop eweight, Vprop vweight, EVprop egroups,
                         VEprop esrcpos, VEprop etgtpos, Graph& g, BGraph& bg,
                         EMat& emat, SamplerMap neighbour_sampler,
                         SamplerMap cavity_neighbour_sampler,
                         bool uniform_moves, bool sequential, bool parallel,
                         bool random_move, double c, size_t nmerges,
                         size_t ntries, Vprop merge_map,
                         overlap_stats_t& overlap_stats,
                         overlap_partition_stats_t& partition_stats,
                         bool verbose, RNG& rng, double& S, size_t& nmoves)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

    size_t B = num_vertices(bg);

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

        //separate into subgroups according to opposing sides
        int u = overlap_stats.get_out_neighbour(v);
        if (u == -1)
            u = overlap_stats.get_in_neighbour(v);
        bundles[b[v]][b[u]].push_back(v);
    }

    std::uniform_int_distribution<size_t> s_rand(0, B - 1);

    nmoves = 0;
    S = 0;

    EntrySet<Graph> m_entries(B);

    half_edge_neighbour_policy<Graph> npolicy(g);
    overlap_partition_stats_t dummy_pstats;

    for (size_t r = 0; r < B; ++r)
    {
        if (groups[r].empty())
            continue;

        // "explode"
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

                    assert(last_b[v] == int(r));

                    // neighbour sampler points to the *block graph*
                    vertex_t w = overlap_stats.sample_half_edge(v, rng);
                    vertex_t t = last_b[w];
                    s = neighbour_sampler[t].sample(rng);
                    if (s == r)
                        s = cavity_neighbour_sampler[s].sample(rng);
                    else
                        s = neighbour_sampler[s].sample(rng);
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
#ifdef HAVE_SPARSEHASH
                dense_hash_map<size_t, vector<size_t>> vs;
                dense_hash_map<size_t, vector<size_t>> rs;
                dense_hash_map<size_t, vector<size_t>> nrs;
                vs.set_empty_key(numeric_limits<size_t>::max());
                rs.set_empty_key(numeric_limits<size_t>::max());
                nrs.set_empty_key(numeric_limits<size_t>::max());
#else
                unordered_map<size_t, vector<size_t>> vs;
                unordered_map<size_t, vector<size_t>> rs;
                unordered_map<size_t, vector<size_t>> nrs;
#endif

                double ddS = 0;
                if (partition_stats.is_enabled())
                {
                    for (auto u : bundle)
                    {
                        size_t vi = overlap_stats.get_node(u);
                        vs[vi].push_back(u);
                        rs[vi].push_back(best_move[u]);
                        nrs[vi].push_back(s);
                    }

                    for (auto& vi : vs)
                    {
                        ddS += partition_stats.get_delta_dl(vi.second,
                                                            rs[vi.first],
                                                            nrs[vi.first],
                                                            deg_corr,
                                                            overlap_stats, g);
                        partition_stats.move_vertex(vi.second,
                                                    rs[vi.first],
                                                    nrs[vi.first],
                                                    deg_corr,
                                                    overlap_stats, g);
                    }
                }

                for (auto u : bundle)
                {
                    ddS += virtual_move(u, s, dense, mrs, mrp, mrm, wr, b,
                                        deg_corr, eweight, vweight, g, bg, emat,
                                        m_entries, overlap_stats, multigraph,
                                        dummy_pstats, npolicy);

                    move_vertex(u, s, mrs, mrp, mrm, wr, b, deg_corr, eweight,
                                vweight, g, bg, emat, overlap_stats,
                                dummy_pstats, npolicy);
                }


                if (ddS < 0 || best_move[bundle[0]] == int(r))
                {
                    for (auto u : bundle)
                        best_move[u] = s;
                    dS += ddS;
                }
                else
                {
                    if (partition_stats.is_enabled())
                    {
                        for (auto& vi : vs)
                        {
                            partition_stats.move_vertex(vi.second,
                                                        nrs[vi.first],
                                                        rs[vi.first],
                                                        deg_corr,
                                                        overlap_stats, g);
                        }
                    }

                    for (auto u : bundle)
                    {
                        move_vertex(u, best_move[u], mrs, mrp, mrm, wr, b,
                                    deg_corr, eweight, vweight, g, bg, emat,
                                    overlap_stats, dummy_pstats, npolicy);
                    }
                }

            }
        }

        for (auto v : groups[r])
        {
            move_vertex(v, r, mrs, mrp, mrm, wr, b, deg_corr, eweight, vweight,
                        g, bg, emat, overlap_stats, partition_stats, npolicy);
        }

        best_move_dS[r] = dS;

        // uniform move
        if (uniform_moves)
        {
            double best_dS = numeric_limits<double>::max();
            for (size_t j = 0; j < ntries; ++j)
            {
                vertex_t s;

                if (!random_move)
                {
                    vertex_t v = uniform_sample(groups[r], rng);

                    // neighbour sampler points to the *block graph*
                    vertex_t w = overlap_stats.sample_half_edge(v, rng);
                    vertex_t t = b[w];
                    s = neighbour_sampler[t].sample(rng);
                    if (s == r)
                        s = cavity_neighbour_sampler[s].sample(rng);
                    else
                        s = neighbour_sampler[s].sample(rng);
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

                if (s == r || wr[s] == 0)
                    continue;

                dS = 0;
                for (auto v: groups[r])
                {
                    if (vweight[v] == 0)
                        continue;

                    dS += virtual_move(v, s, dense, mrs, mrp, mrm, wr, b,
                                       deg_corr, eweight, vweight, g, bg, emat,
                                       m_entries, overlap_stats, multigraph,
                                       partition_stats, npolicy);

                    move_vertex(v, s, mrs, mrp, mrm, wr, b, deg_corr, eweight,
                                vweight, g, bg, emat, overlap_stats,
                                partition_stats, npolicy);
                }

                if (dS < best_move_dS[r])
                {
                    best_move_dS[r] = dS;
                    for (auto v: groups[r])
                        best_move[v] = s;

                }
            }

            for (auto v: groups[r])
            {
                move_vertex(v, r, mrs, mrp, mrm, wr, b, deg_corr, eweight,
                            vweight, g, bg, emat, overlap_stats,
                            partition_stats, npolicy);
            }

            best_dS = min(best_dS, dS);
        }
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

            double dS = virtual_move(v, s, dense, mrs, mrp, mrm, wr, b,
                                     deg_corr, eweight, vweight, g, bg, emat,
                                     m_entries, overlap_stats, multigraph,
                                     partition_stats, npolicy);

            move_vertex(v, s, mrs, mrp, mrm, wr, b, deg_corr, eweight, vweight,
                        g, bg, emat, overlap_stats, partition_stats, npolicy);
            S += dS;
        }

        //assert (wr[r] == 0);
        if (wr[r] == 0)
            ++nmoves;
    }
}

} // namespace graph_tool

#endif // GRAPH_BLOCKMODEL_OVERLAP_HH
