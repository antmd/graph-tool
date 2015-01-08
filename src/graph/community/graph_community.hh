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

#ifndef GRAPH_COMMUNITY_HH
#define GRAPH_COMMUNITY_HH

#include <unordered_set>
#include <tuple>
#include <iostream>
#include <fstream>
#include <iomanip>

#define BOOST_DISABLE_ASSERTS
#include "boost/multi_array.hpp"

#include "graph_util.hh"
#include "graph_properties.hh"
#include "random.hh"

namespace graph_tool
{

using namespace std;
using namespace boost;

using std::unordered_map;
using std::unordered_set;

// computes the community structure through a spin glass system with
// simulated annealing

template <template <class G, class CommunityMap> class NNKS>
struct get_communities
{
    template <class Graph, class VertexIndex, class WeightMap,
              class CommunityMap>
    void operator()(const Graph& g, VertexIndex vertex_index, WeightMap weights,
                    CommunityMap s, double gamma, size_t n_iter,
                    pair<double, double> Tinterval, size_t n_spins, rng_t& rng,
                    pair<bool, string> verbose) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;
        typedef typename property_traits<WeightMap>::key_type weight_key_t;


        auto random = std::bind(std::uniform_real_distribution<>(), std::ref(rng));

        stringstream out_str;
        ofstream out_file;
        if (verbose.second != "")
        {
            out_file.open(verbose.second.c_str());
            if (!out_file.is_open())
                throw IOException("error opening file " + verbose.second +
                                  " for writing");
            out_file.exceptions (ifstream::eofbit | ifstream::failbit |
                                 ifstream::badbit);
        }

        double Tmin = Tinterval.first;
        double Tmax = Tinterval.second;

        unordered_map<size_t, size_t> Ns; // spin histogram
        CommunityMap temp_s(vertex_index, num_vertices(g));

        // init spins from [0,N-1] and global info
        uniform_int_distribution<size_t> sample_spin(0, n_spins-1);
        typename graph_traits<Graph>::vertex_iterator vi,vi_end;
        for (tie(vi,vi_end) = vertices(g); vi != vi_end; ++vi)
        {
            s[*vi] = temp_s[*vi] = sample_spin(rng);
            Ns[s[*vi]]++;
        }

        NNKS<Graph,CommunityMap> Nnnks(g, s); // this will retrieve the expected
                                              // number of neighbours with given
                                              // spin, as a function of degree

        // define cooling rate so that temperature starts at Tmax at temp_count
        // == 0 and reaches Tmin at temp_count == n_iter - 1
        if (Tmin < numeric_limits<double>::epsilon())
            Tmin = numeric_limits<double>::epsilon();
        double cooling_rate = -(log(Tmin)-log(Tmax))/(n_iter-1);

        // start the annealing
        for (size_t temp_count = 0; temp_count < n_iter; ++temp_count)
        {
            double T = Tmax*exp(-cooling_rate*temp_count);
            double E = 0;

            vector<std::tuple<size_t, size_t, size_t> > updates;

            // sample a new spin for every vertex
            int NV = num_vertices(g),i;
            #pragma omp parallel for default(shared) private(i)\
                reduction(+:E) schedule(runtime) if (NV > 100)
            for (i = 0; i < NV; ++i)
            {
                vertex_t v = vertex(i, g);
                if (v == graph_traits<Graph>::null_vertex())
                    continue;

                size_t new_s;
                {
                    #pragma omp critical
                    new_s = sample_spin(rng);
                }

                unordered_map<size_t, double> ns; // number of neighbours with a
                                                  // given spin 's' (weighted)

                // neighborhood spins info
                typename graph_traits<Graph>::out_edge_iterator e, e_end;
                for (tie(e,e_end) = out_edges(v,g); e != e_end; ++e)
                {
                    vertex_t t = target(*e, g);
                    if (t != v)
                        ns[s[t]] += get(weights, weight_key_t(*e));
                }

                size_t k = out_degree_no_loops(v, g);

                double curr_e = gamma*Nnnks(k,s[v]) - ns[s[v]];
                double new_e = gamma*Nnnks(k,new_s) - ns[new_s];

                double r;
                {
                    #pragma omp critical
                    r = random();
                }

                if (new_e < curr_e || r < exp(-(new_e - curr_e)/T))
                {
                    temp_s[v] = new_s;
                    curr_e = new_e;
                    {
                        #pragma omp critical
                        updates.push_back(std::make_tuple(k, size_t(s[v]),
                                                          new_s));
                        Ns[s[v]]--;
                        Ns[new_s]++;
                    }
                }
                else
                {
                    temp_s[v] = s[v];
                }
                E += curr_e;
            }
            swap(s, temp_s);

            for (typeof(updates.begin()) iter = updates.begin();
                 iter != updates.end(); ++iter)
                Nnnks.Update(std::get<0>(*iter), std::get<1>(*iter),
                             std::get<2>(*iter));

            if (verbose.first)
            {
                for (size_t j = 0; j < out_str.str().length(); ++j)
                    cout << "\b";
                out_str.str("");
                size_t ns = 0;
                for (typeof(Ns.begin()) iter = Ns.begin(); iter != Ns.end();
                     ++iter)
                    if (iter->second > 0)
                        ns++;
                out_str << setw(lexical_cast<string>(n_iter).size())
                        << temp_count << " of " << n_iter
                        << " (" << setw(2) << (temp_count+1)*100/n_iter
                        << "%) " << "temperature: " << setw(14)
                        << setprecision(10) << T << " spins: "
                        << ns << " energy: " << E;
                cout << out_str.str() << flush;
            }
            if (verbose.second != "")
            {
                try
                {
                    size_t ns = 0;
                    for (typeof(Ns.begin()) iter = Ns.begin(); iter != Ns.end();
                         ++iter)
                        if (iter->second > 0)
                            ns++;
                    out_file << temp_count << "\t" << setprecision(10) << T
                             << "\t" << ns << "\t" << E << endl;
                }
                catch (ifstream::failure e)
                {
                    throw IOException("error writing to file " +
                                      verbose.second + ": " + e.what());
                }
            }
        }

        if (n_iter % 2 != 0)
        {
            int NV = num_vertices(g), i;
            #pragma omp parallel for default(shared) private(i)\
                schedule(runtime) if (NV > 100)
            for (i = 0; i < NV; ++i)
            {
                vertex_t v = vertex(i, g);
                if (v == graph_traits<Graph>::null_vertex())
                    continue;
                temp_s[v] = s[v];
            }
        }

        // rename spins, starting from zero
        unordered_map<size_t,size_t> spins;
        for (tie(vi,vi_end) = vertices(g); vi != vi_end; ++vi)
        {
            if (spins.find(s[*vi]) == spins.end())
                spins[s[*vi]] = spins.size() - 1;
            s[*vi] = spins[s[*vi]];
        }

    }
};

template <class Graph, class CommunityMap>
class NNKSErdosReyni
{
public:
    NNKSErdosReyni(const Graph &g, CommunityMap s)
    {
        size_t N = 0;
        double _avg_k = 0.0;
        typename graph_traits<Graph>::vertex_iterator v,v_end;
        for (tie(v,v_end) = vertices(g); v != v_end; ++v)
        {
            size_t k = out_degree_no_loops(*v,g);
            _avg_k += k;
            N++;
            _Ns[s[*v]]++;
        }
        _p = _avg_k/(N*N);
    }

    void Update(size_t, size_t old_s, size_t s)
    {
        _Ns[old_s]--;
        if (_Ns[old_s] == 0)
            _Ns.erase(old_s);
        _Ns[s]++;
    }

    double operator()(size_t, size_t s) const
    {
        size_t ns = 0;
        typeof(_Ns.begin()) iter = _Ns.find(s);
        if (iter != _Ns.end())
            ns = iter->second;
        return _p*ns;
    }

private:
    double _p;
    unordered_map<size_t,size_t> _Ns;
};

template <class Graph, class CommunityMap>
class NNKSUncorr
{
public:
    NNKSUncorr(const Graph &g, CommunityMap s): _g(g), _K(0)
    {
        typename graph_traits<Graph>::vertex_iterator v,v_end;
        for (tie(v,v_end) = vertices(_g); v != v_end; ++v)
        {
            size_t k = out_degree_no_loops(*v, _g);
            _K += k;
            _Ks[s[*v]] += k;
        }
    }

    void Update(size_t k, size_t old_s, size_t s)
    {
        _Ks[old_s] -= k;
        if (_Ks[old_s] == 0)
            _Ks.erase(old_s);
        _Ks[s] += k;
    }

    double operator()(size_t k, size_t s) const
    {
        size_t ks = 0;
        typeof(_Ks.begin()) iter = _Ks.find(s);
        if (iter != _Ks.end())
            ks = iter->second;
        return k*ks/double(_K);
    }

private:
    const Graph& _g;
    size_t _K;
    unordered_map<size_t,size_t> _Ks;
};

template <class Graph, class CommunityMap>
class NNKSCorr
{
public:
    NNKSCorr(const Graph &g, CommunityMap s): _g(g)
    {
        unordered_set<size_t> spins;

        typename graph_traits<Graph>::vertex_iterator v,v_end;
        for (tie(v,v_end) = vertices(_g); v != v_end; ++v)
        {
            size_t k = out_degree_no_loops(*v, _g);
            _Nk[k]++;
            _Nks[k][s[*v]]++;
            spins.insert(s[*v]);
        }

        size_t E = 0;
        typename graph_traits<Graph>::edge_iterator e,e_end;
        for (tie(e,e_end) = edges(_g); e != e_end; ++e)
        {
            typename graph_traits<Graph>::vertex_descriptor src, tgt;

            src = source(*e,g);
            tgt = target(*e,g);
            if (src != tgt)
            {
                size_t k1 = out_degree_no_loops(src, g);
                size_t k2 = out_degree_no_loops(tgt, g);
                _Pkk[k1][k2]++;
                _Pkk[k2][k1]++;
                E++;
            }
        }

        for (typeof(_Pkk.begin()) iter1 = _Pkk.begin(); iter1 != _Pkk.end();
             ++iter1)
        {
            double sum = 0;
            for (typeof(iter1->second.begin()) iter2 = iter1->second.begin();
                 iter2 != iter1->second.end(); ++iter2)
                sum += iter2->second;
            for (typeof(iter1->second.begin()) iter2 = iter1->second.begin();
                 iter2 != iter1->second.end(); ++iter2)
                iter2->second /= sum;
        }

        for (typeof(_Nk.begin()) k_iter = _Nk.begin(); k_iter != _Nk.end();
             ++k_iter)
        {
            size_t k1 = k_iter->first;
            _degs.push_back(k1);
            for (typeof(spins.begin()) s_iter = spins.begin();
                 s_iter != spins.end(); ++s_iter)
                for (typeof(_Nk.begin()) k_iter2 = _Nk.begin();
                     k_iter2 != _Nk.end(); ++k_iter2)
                {
                    size_t k2 = k_iter2->first;
                    if (_Nks[k2].find(*s_iter) != _Nks[k2].end())
                        _NNks[k1][*s_iter] +=
                            k1*_Pkk[k1][k2] * _Nks[k2][*s_iter]/double(_Nk[k2]);
                }
        }
    }

    void Update(size_t k, size_t old_s, size_t s)
    {
        int i, NK = _degs.size();
        #pragma omp parallel for default(shared) private(i) schedule(runtime) if (NK > 100)
        for (i = 0; i < NK; ++i)
        {
            size_t k1 = _degs[i], k2 = k;
            if (_Pkk.find(k1) == _Pkk.end())
                continue;
            if (_Pkk.find(k1)->second.find(k2) == _Pkk.find(k1)->second.end())
                continue;
            unordered_map<size_t,double>& NNks_k1 = _NNks[k1];
            double Pk1k2 = _Pkk[k1][k2];
            unordered_map<size_t,size_t>& Nksk2 = _Nks[k2];
            double Nk2 = _Nk[k2];
            NNks_k1[old_s] -=  k1*Pk1k2 * Nksk2[old_s]/Nk2;
            if (NNks_k1[old_s] == 0.0)
                NNks_k1.erase(old_s);
            if (Nksk2.find(s) != Nksk2.end())
                NNks_k1[s] -=  k1*Pk1k2 * Nksk2[s]/Nk2;
            if (NNks_k1[s] == 0.0)
                NNks_k1.erase(s);
        }

        _Nks[k][old_s]--;
        if (_Nks[k][old_s] == 0)
            _Nks[k].erase(old_s);
        _Nks[k][s]++;

        #pragma omp parallel for default(shared) private(i) schedule(runtime) if (NK > 100)
        for (i = 0; i < NK; ++i)
        {
            size_t k1 = _degs[i], k2 = k;
            if (_Pkk.find(k1) == _Pkk.end())
                continue;
            if (_Pkk.find(k1)->second.find(k2) == _Pkk.find(k1)->second.end())
                continue;
            unordered_map<size_t,double>& NNks_k1 = _NNks[k1];
            double Pk1k2 = _Pkk[k1][k2];
            unordered_map<size_t,size_t>& Nksk2 = _Nks[k2];
            double Nk2 = _Nk[k2];
            NNks_k1[old_s] +=  k1*Pk1k2 * Nksk2[old_s]/Nk2;
            if (NNks_k1[old_s] == 0.0)
                NNks_k1.erase(old_s);
            NNks_k1[s] +=  k1*Pk1k2 * Nksk2[s]/Nk2;
        }

    }

    double operator()(size_t k, size_t s) const
    {
        const unordered_map<size_t,double>& nnks = _NNks.find(k)->second;
        const typeof(nnks.begin()) iter = nnks.find(s);
        if (iter != nnks.end())
            return iter->second;
        return 0.0;
    }

private:
    const Graph& _g;
    vector<size_t> _degs;
    unordered_map<size_t,size_t> _Nk;
    unordered_map<size_t,unordered_map<size_t,double> > _Pkk;
    unordered_map<size_t,unordered_map<size_t,size_t> > _Nks;
    unordered_map<size_t,unordered_map<size_t,double> > _NNks;
};

enum comm_corr_t
{
    ERDOS_REYNI,
    UNCORRELATED,
    CORRELATED
};

struct get_communities_selector
{
    get_communities_selector(comm_corr_t corr,
                             GraphInterface::vertex_index_map_t index)
        : _corr(corr), _index(index) {}
    comm_corr_t _corr;
    GraphInterface::vertex_index_map_t _index;

    template <class Graph, class WeightMap, class CommunityMap>
    void operator()(const Graph& g, WeightMap weights, CommunityMap s,
                    double gamma, size_t n_iter, pair<double, double> Tinterval,
                    size_t Nspins, rng_t& rng, pair<bool, string> verbose)
        const
    {
        switch (_corr)
        {
        case ERDOS_REYNI:
            get_communities<NNKSErdosReyni>()(g, _index, weights, s, gamma,
                                              n_iter, Tinterval, Nspins, rng,
                                              verbose);
            break;
        case UNCORRELATED:
            get_communities<NNKSUncorr>()(g, _index, weights, s, gamma, n_iter,
                                          Tinterval, Nspins, rng, verbose);
            break;
        case CORRELATED:
            get_communities<NNKSCorr>()(g, _index, weights, s, gamma, n_iter,
                                        Tinterval, Nspins, rng, verbose);
            break;
        }
    }
};

// get Newman's modularity of a given community partition
struct get_modularity
{
    template <class Graph, class WeightMap, class CommunityMap>
    void operator()(const Graph& g, WeightMap weights, CommunityMap b, double& Q) const
    {
        typedef typename property_traits<WeightMap>::key_type weight_key_t;
        typedef typename property_traits<WeightMap>::value_type weight_val_t;
        typedef typename property_traits<CommunityMap>::value_type s_val_t;

        vector<double> er, err;
        double W = 0;

        typename graph_traits<Graph>::edge_iterator e, e_end;
        for (tie(e, e_end) = edges(g); e != e_end; ++e)
        {
            size_t r = get(b, source(*e,g));
            size_t s = get(b, target(*e,g));

            double w = get(weights, *e);
            W += 2 * w;

            if (er.size() <= r)
                er.resize(r + 1);
            er[r] += w;

            if (er.size() <= s)
                er.resize(s + 1);
            er[s] += w;

            if (r == s)
            {
                if (err.size() <= r)
                    err.resize(r + 1);
                err[r] += 2 * w;
            }
        }

        Q = 0;
        for (size_t r = 0; r < er.size(); ++r)
        {
            if (err.size() <= r)
                err.resize(r + 1);

            Q += err[r] - (er[r] * er[r]) / W;
        }
        Q /= W;

    }
};

} // graph_tool namespace

#endif //GRAPH_COMMUNITY_HH
