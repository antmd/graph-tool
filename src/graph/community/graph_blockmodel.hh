// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2014 Tiago de Paula Peixoto <tiago@skewed.de>
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

#include "config.h"
#include <unordered_set>
#include <unordered_map>
#include <tuple>

#ifdef HAVE_SPARSEHASH
#include <sparsehash/dense_hash_set>
#include <sparsehash/dense_hash_map>
#endif

#include "../generation/sampler.hh"
#include "../generation/dynamic_sampler.hh"

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

template <class Type>
__attribute__((always_inline))
inline double safelog(Type x)
{
    if (x == 0)
        return 0;
    return log(x);
}

__attribute__((always_inline))
inline double safelog(size_t x);

__attribute__((always_inline))
inline double xlogx(size_t x);

__attribute__((always_inline))
inline double lgamma_fast(size_t x);


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
    void operator()(Eprop mrs, Vprop mrp, Vprop mrm, Vprop wr, bool deg_corr,
                    Graph& g, double& S) const
    {
        S = 0;
        typename graph_traits<Graph>::edge_iterator e, e_end;
        for (tie(e, e_end) = edges(g); e != e_end; ++e)
            S += eterm(source(*e, g), target(*e, g), mrs[*e], g);
        typename graph_traits<Graph>::vertex_iterator v, v_end;
        for (tie(v, v_end) = vertices(g); v != v_end; ++v)
            S += vterm(*v, mrp[*v], mrm[*v], wr[*v], deg_corr, g);
    }
};


//
// Dense entropy
//

double lbinom(double N, double k)
{
    return lgamma(N + 1) - lgamma(N - k + 1) - lgamma(k + 1);
}

double lbinom_fast(int N, int k)
{
    return lgamma_fast(N + 1) - lgamma_fast(N - k + 1) - lgamma_fast(k + 1);
}


// "edge" term of the entropy
template <class Graph>
__attribute__((always_inline))
inline double eterm_dense(size_t r, size_t s, int mrs, int wr_r,
                          int wr_s, bool multigraph, const Graph& g)
{
    int nrns;
    int ers = mrs;

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
        S = lbinom_fast(nrns + ers - 1, ers);
    else
        S = lbinom_fast(nrns, ers);
    return S;
}

struct entropy_dense
{
    template <class Graph, class Eprop, class Vprop>
    void operator()(Eprop mrs, Vprop wr, bool multigraph, Graph& g, double& S) const
    {
        S = 0;
        typename graph_traits<Graph>::edge_iterator e, e_end;
        for (tie(e, e_end) = edges(g); e != e_end; ++e)
        {
            typename graph_traits<Graph>::vertex_descriptor r, s;
            r = source(*e, g);
            s = target(*e, g);
            S += eterm_dense(r, s, mrs[*e], wr[r], wr[s], multigraph, g);
        }
    }
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
        typename graph_traits<Graph>::edge_iterator e, e_end;
        for (tie(e, e_end) = edges(g); e != e_end; ++e)
        {
            if (source(*e, g) >= B || target(*e, g) >= B)
                throw GraphException("incorrect number of blocks when creating emat!");
            emat[source(*e, g)][target(*e, g)] = make_pair(*e, true);
            if (!is_directed::apply<Graph>::type::value)
                emat[target(*e, g)][source(*e, g)] = make_pair(*e, true);
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
    const typename get_ehash_t::apply<Graph>::map_t& map = ehash[r];
    typeof(map.begin()) iter = map.find(s);
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
        typename graph_traits<Graph>::vertex_iterator v, v_end;
        for (tie(v, v_end) = vertices(g); v != v_end; ++v)
        {
            emat[*v].set_empty_key(numeric_limits<vertex_t>::max());
            emat[*v].set_deleted_key(num_vertices(g));
        }
#endif

        typename graph_traits<Graph>::edge_iterator e, e_end;
        for (tie(e, e_end) = edges(g); e != e_end; ++e)
            put_me(source(*e, g), target(*e, g), *e, emat, g);

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

// remove a vertex from its current block
template <class Graph, class BGraph, class Eprop, class Vprop, class EWprop,
          class VWprop, class EMat>
void remove_vertex(size_t v, Eprop& mrs, Vprop& mrp, Vprop& mrm, Vprop& wr,
                   Vprop& b, const EWprop& eweight, const VWprop& vweight,
                   Graph& g, BGraph& bg, EMat& emat)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    vertex_t r = b[v];

    bool self_count = false;
    typename graph_traits<Graph>::out_edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = out_edges(v, g); ei != ei_end; ++ei)
    {
        typename graph_traits<Graph>::edge_descriptor e = *ei;
        vertex_t u = target(e, g);
        if (u == v && !is_directed::apply<Graph>::type::value)
        {
            if (self_count)
                continue;
            self_count = true;
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

    if (is_directed::apply<Graph>::type::value)
    {
        typename in_edge_iteratorS<Graph>::type ie, ie_end;
        for (tie(ie, ie_end) = in_edge_iteratorS<Graph>::get_edges(vertex(v, g), g);
             ie != ie_end; ++ie)
        {
            typename graph_traits<Graph>::edge_descriptor e = *ie;
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
    }

    wr[r] -= vweight[v];
}

// add a vertex to block rr
template <class Graph, class BGraph, class Eprop, class Vprop, class EWprop,
          class VWprop,class EMat>
void add_vertex(size_t v, size_t r, Eprop& mrs, Vprop& mrp, Vprop& mrm,
                Vprop& wr, Vprop& b, const EWprop& eweight,
                const VWprop& vweight, Graph& g, BGraph& bg, EMat& emat)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

    bool self_count = false;
    typename graph_traits<Graph>::out_edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = out_edges(vertex(v, g), g); ei != ei_end; ++ei)
    {
        typename graph_traits<Graph>::edge_descriptor e = *ei;
        vertex_t u = target(e, g);
        if (u == v && !is_directed::apply<Graph>::type::value )
        {
            if (self_count)
                continue;
            self_count = true;
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

    typename in_edge_iteratorS<Graph>::type ie, ie_end;
    for (tie(ie, ie_end) = in_edge_iteratorS<Graph>::get_edges(vertex(v, g), g);
         ie != ie_end; ++ie)
    {
        typename graph_traits<Graph>::edge_descriptor e = *ie;
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

    wr[r] += vweight[v];
    b[v] = r;
}

// move a vertex from its current block to block nr
template <class Graph, class BGraph, class Eprop, class Vprop, class EWprop,
          class VWprop, class EMat>
void move_vertex(size_t v, size_t nr, Eprop& mrs, Vprop& mrp, Vprop& mrm,
                 Vprop& wr, Vprop& b, bool deg_corr, const EWprop& eweight,
                 const VWprop& vweight, Graph& g, BGraph& bg, EMat& emat)

{
    remove_vertex(v, mrs, mrp, mrm, wr, b, eweight, vweight, g, bg, emat);
    add_vertex(v, nr, mrs, mrp, mrm, wr, b, eweight, vweight, g, bg, emat);
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
template <class Graph, class BGraph, class Vertex, class Vprop, class Eprop>
void move_entries(Vertex v, Vertex nr, Vprop b, Eprop eweights, Graph& g,
                  BGraph& bg, EntrySet<Graph>& m_entries)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    vertex_t r = b[v];

    m_entries.SetMove(r, nr);

    bool self_count = false;
    typename graph_traits<Graph>::out_edge_iterator e, e_end;
    for (tie(e, e_end) = out_edges(v, g); e != e_end; ++e)
    {
        vertex_t u = target(*e, g);
        if (u == v && !is_directed::apply<Graph>::type::value)
        {
            if (self_count)
                continue;
            self_count = true;
        }
        vertex_t s = b[u];
        int ew = eweights[*e];
        //assert(ew > 0);

        m_entries.InsertDeltaTarget(r, s, -ew);

        //insert_m_entry(r,  s, -ew, m_entries, m_entries_set, g);
        if (u == v)
            s = nr;
        m_entries.InsertDeltaTarget(nr, s, +ew);

        //insert_m_entry(nr, s, +ew, m_entries, m_entries_set, g);
    }

    typename in_edge_iteratorS<Graph>::type ie, ie_end;
    for (tie(ie, ie_end) = in_edge_iteratorS<Graph>::get_edges(v, g);
         ie != ie_end; ++ie)
    {
        vertex_t u = source(*ie, g);
        if (u == v)
            continue;
        vertex_t s = b[u];
        int ew = eweights[*ie];

        m_entries.InsertDeltaSource(s,  r, -ew);
        m_entries.InsertDeltaSource(s, nr, +ew);

        // insert_m_entry(s,  r, -ew, m_entries, m_entries_set, g);
        // insert_m_entry(s, nr, +ew, m_entries, m_entries_set, g);
    }

    //clear_entries_set(m_entries, m_entries_set);
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
        //assert(ers + delta >= 0);
        dS += eterm(er, es, ers + d, bg) - eterm(er, es, ers, bg);
    }
    return dS;
}

// compute the entropy difference of a virtual move of vertex r to block nr
template <class Graph, class BGraph, class Eprop, class Vprop, class EWprop,
          class VWprop, class EMat>
double virtual_move(size_t v, size_t nr, Eprop& mrs, Vprop& mrp, Vprop& mrm,
                    Vprop& wr, Vprop& b, bool deg_corr, const EWprop& eweight,
                    const VWprop& vweight, Graph& g, BGraph& bg, EMat& emat,
                    EntrySet<Graph>& m_entries)

{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    vertex_t r = b[v];

    if (r == nr)
        return 0.;

    m_entries.Clear();
    move_entries(v, nr, b, eweight, g, bg, m_entries);
    double dS = entries_dS(m_entries, mrs, bg, emat);
    int kout = out_degreeS()(v, g, eweight);
    int kin = kout;
    if (is_directed::apply<Graph>::type::value)
        kin = in_degreeS()(v, g, eweight);
    //assert(mrm[r]  - kin >= 0);
    //assert(mrp[r]  - kout >= 0);
    dS += vterm(r,  mrp[r]  - kout, mrm[r]  - kin, wr[r]  - vweight[v], deg_corr, bg);
    dS += vterm(nr, mrp[nr] + kout, mrm[nr] + kin, wr[nr] + vweight[v], deg_corr, bg);
    dS -= vterm(r,  mrp[r]        , mrm[r]       , wr[r]              , deg_corr, bg);
    dS -= vterm(nr, mrp[nr]       , mrm[nr]      , wr[nr]             , deg_corr, bg);
    return dS;
}


// compute the entropy difference of a virtual move of vertex r to block nr
template <class Graph, class BGraph, class Eprop, class Vprop, class EWprop,
          class VWprop, class EMat>
double virtual_move_dense(size_t v, size_t nr, Eprop& mrs, Vprop& mrp,
                          Vprop& mrm, Vprop& wr, Vprop& b, bool deg_corr,
                          const EWprop& eweight, const VWprop& vweight,
                          Graph& g, BGraph& bg, EMat& emat, bool multigraph)

{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    vertex_t r = b[v];

    if (r == nr)
        return 0;

    if (deg_corr)
        throw GraphException("Dense entropy for degree corrected model not implemented!");

    vector<int> deltap(num_vertices(bg), 0);
    int deltal = 0;
    bool self_count = false;
    typename graph_traits<Graph>::out_edge_iterator e, e_end;
    for (tie(e, e_end) = out_edges(v, g); e != e_end; ++e)
    {
        vertex_t u = target(*e, g);
        vertex_t s = b[u];
        if (u == v)
        {
            if (self_count)
                continue;
            self_count = true;
            deltal += eweight[*e];
        }
        else
        {
            deltap[s] += eweight[*e];
        }
    }

    vector<int> deltam(num_vertices(bg), 0);
    typename in_edge_iteratorS<Graph>::type ie, ie_end;
    for (tie(ie, ie_end) = in_edge_iteratorS<Graph>::get_edges(v, g);
         ie != ie_end; ++ie)
    {
        vertex_t u = source(*ie, g);
        if (u == v)
            continue;
        vertex_t s = b[u];
        deltam[s] += eweight[*ie];
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
                Si += eterm_dense(r,  s, ers,              wr[r],               wr[s], multigraph, bg);
                Sf += eterm_dense(r,  s, ers - deltap[s],  wr[r] - vweight[v],  wr[s], multigraph, bg);
                Si += eterm_dense(nr, s, enrs,             wr[nr],              wr[s], multigraph, bg);
                Sf += eterm_dense(nr, s, enrs + deltap[s], wr[nr] + vweight[v], wr[s], multigraph, bg);
            }

            if (s == r)
            {
                Si += eterm_dense(r, r, ers,                      wr[r],              wr[r],              multigraph, bg);
                Sf += eterm_dense(r, r, ers - deltap[r] - deltal, wr[r] - vweight[v], wr[r] - vweight[v], multigraph, bg);
            }

            if (s == nr)
            {
                Si += eterm_dense(nr, nr, enrs,                       wr[nr],              wr[nr],              multigraph, bg);
                Sf += eterm_dense(nr, nr, enrs + deltap[nr] + deltal, wr[nr] + vweight[v], wr[nr] + vweight[v], multigraph, bg);

                Si += eterm_dense(r, nr, ers,                          wr[r],              wr[nr], multigraph, bg);
                Sf += eterm_dense(r, nr, ers - deltap[nr] + deltap[r], wr[r] - vweight[v], wr[nr] + vweight[v], multigraph, bg);
            }
        }
        else
        {
            int esr = get_mrs(s, r, mrs, emat, bg);
            int esnr = get_mrs(s, nr, mrs, emat, bg);

            if (s != nr && s != r)
            {
                Si += eterm_dense(r, s, ers            , wr[r]             , wr[s]             , multigraph, bg);
                Sf += eterm_dense(r, s, ers - deltap[s], wr[r] - vweight[v], wr[s]             , multigraph, bg);
                Si += eterm_dense(s, r, esr            , wr[s]             , wr[r]             , multigraph, bg);
                Sf += eterm_dense(s, r, esr - deltam[s], wr[s]             , wr[r] - vweight[v], multigraph, bg);

                Si += eterm_dense(nr, s, enrs            , wr[nr]             , wr[s]              , multigraph, bg);
                Sf += eterm_dense(nr, s, enrs + deltap[s], wr[nr] + vweight[v], wr[s]              , multigraph, bg);
                Si += eterm_dense(s, nr, esnr            , wr[s]              , wr[nr]             , multigraph, bg);
                Sf += eterm_dense(s, nr, esnr + deltam[s], wr[s]              , wr[nr] + vweight[v], multigraph, bg);
            }

            if(s == r)
            {
                Si += eterm_dense(r, r, ers                                  , wr[r]             , wr[r]             , multigraph, bg);
                Sf += eterm_dense(r, r, ers - deltap[r]  - deltam[r] - deltal, wr[r] - vweight[v], wr[r] - vweight[v], multigraph, bg);

                Si += eterm_dense(r, nr, esnr                         , wr[r]             , wr[nr]             , multigraph, bg);
                Sf += eterm_dense(r, nr, esnr - deltap[nr] + deltam[r], wr[r] - vweight[v], wr[nr] + vweight[v], multigraph, bg);
            }

            if(s == nr)
            {
                Si += eterm_dense(nr, nr, esnr                                   , wr[nr]             , wr[nr]             , multigraph, bg);
                Sf += eterm_dense(nr, nr, esnr + deltap[nr] + deltam[nr] + deltal, wr[nr] + vweight[v], wr[nr] + vweight[v], multigraph, bg);

                Si += eterm_dense(nr, r, esr                         , wr[nr]             , wr[r]             , multigraph, bg);
                Sf += eterm_dense(nr, r, esr + deltap[r] - deltam[nr], wr[nr] + vweight[v], wr[r] - vweight[v], multigraph, bg);
            }
        }
    }
    return Sf - Si;
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
    typename graph_traits<Graph>::vertex_descriptor u, v;
    u = source(e, g);
    v = target(e, g);
    return min(u, v);
}

template <class Edge, class Graph>
inline typename graph_traits<Graph>::vertex_descriptor
get_target(const Edge& e, const Graph &g)
{
    typename graph_traits<Graph>::vertex_descriptor u, v;
    u = source(e, g);
    v = target(e, g);
    return max(u, v);
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
                      VertexIndex vertex_index, bool weighted)
    {
        if (weighted)
        {
            typedef typename get_sampler<Graph, mpl::true_>::type sampler_t;
            typedef typename property_map_type::apply<sampler_t,
                                                      VertexIndex>::type vemap_t;
            vemap_t egroups_checked(vertex_index);
            oegroups = egroups_checked;
            build_dispatch(b, egroups_checked.get_unchecked(num_vertices(g)),
                           esrcpos, etgtpos, eweight, g, vertex_index,
                           mpl::true_());
        }
        else
        {
            typedef typename get_sampler<Graph, mpl::false_>::type sampler_t;
            typedef typename property_map_type::apply<sampler_t,
                                                      VertexIndex>::type vemap_t;
            vemap_t egroups_checked(vertex_index);
            oegroups = egroups_checked;
            build_dispatch(b, egroups_checked.get_unchecked(num_vertices(g)),
                           esrcpos, etgtpos, eweight, g, vertex_index,
                           mpl::true_());
        }
    }

    template <class Eprop, class Vprop, class VEprop, class Graph, class VertexIndex, class Egroups>
    static void build_dispatch(Vprop b, Egroups egroups, VEprop esrcpos,
                               VEprop etgtpos, Eprop eweight, Graph& g,
                               VertexIndex vertex_index, mpl::true_)
    {

        typename graph_traits<Graph>::edge_iterator e, e_end;
        for (tie(e, e_end) = edges(g); e != e_end; ++e)
        {
            size_t r = b[get_source(*e, g)];
            auto& r_elist = egroups[r];
            esrcpos[*e] = insert_edge(std::make_tuple(*e, true), r_elist, eweight[*e]);

            size_t s = b[get_target(*e, g)];
            auto& s_elist = egroups[s];
            etgtpos[*e] = insert_edge(std::make_tuple(*e, false), s_elist, eweight[*e]);
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
    static void remove_egroups(Vertex v, Vertex r, Eprop eweight,
                               EVprop egroups, VEprop esrcpos, VEprop etgtpos,
                               Graph& g)
    {
        typedef Vertex vertex_t;
        bool self_count = false;

        //update the half-edge lists
        typename all_edges_iteratorS<Graph>::type e, e_end;
        for (tie(e, e_end) = all_edges_iteratorS<Graph>::get_edges(v, g);
             e != e_end; ++e)
        {
            vertex_t src = get_source(*e, g);
            vertex_t tgt = get_target(*e, g);

            bool is_src = (src == v);

            // self-loops will appear twice
            if (src == tgt)
            {
                is_src = !self_count;
                self_count = true;
            }

            auto& elist = egroups[r];
            size_t pos = (is_src) ? esrcpos[*e] : etgtpos[*e];
            remove_edge(pos, esrcpos, etgtpos, elist);
        }
    }

    template <class Vertex, class Graph, class EVprop, class Eprop, class VEprop>
    static void add_egroups(Vertex v, Vertex s, Eprop eweight, EVprop egroups,
                            VEprop esrcpos, VEprop etgtpos, Graph& g)
    {
        typedef Vertex vertex_t;
        bool self_count = false;

        //update the half-edge lists
        typename all_edges_iteratorS<Graph>::type e, e_end;
        for (tie(e, e_end) = all_edges_iteratorS<Graph>::get_edges(v, g);
             e != e_end; ++e)
        {
            vertex_t src = get_source(*e, g);
            vertex_t tgt = get_target(*e, g);

            bool is_src = (src == v);

            // self-loops will appear twice
            if (src == tgt)
            {
                is_src = !self_count;
                self_count = true;
            }

            auto& elist = egroups[s];
            typedef typename tuple_element<0, typename property_traits<EVprop>::value_type::value_type>::type e_type;
            size_t pos = insert_edge(std::make_tuple(e_type(*e), is_src),
                                     elist, size_t(eweight[*e]));
            if (is_src)
                esrcpos[*e] = pos;
            else
                etgtpos[*e] = pos;
        }
    }

    template <class Vertex, class Graph, class EVprop, class Eprop, class VEprop>
    static void update_egroups(Vertex v, Vertex r, Vertex s, Eprop eweight,
                               EVprop egroups, VEprop esrcpos, VEprop etgtpos,
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
          class BGraph>
inline double
get_move_prob(Vertex v, Vertex r, Vertex s, double c, Vprop b, Eprop mrs,
              Vprop mrp, Vprop mrm, Emat& emat, Eprop eweight, Graph& g,
              BGraph& bg, EntrySet<Graph>& m_entries, bool reverse)
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

    typename all_edges_iteratorS<Graph>::type e, e_end;
    for (tie(e, e_end) = all_edges_iteratorS<Graph>::get_edges(v, g);
         e != e_end; ++e)
    {
        vertex_t u = target(*e, g);
        if (is_directed::apply<Graph>::type::value && u == v)
            u = source(*e, g);
        vertex_t t = b[u];
        if (u == v)
            t = r;
        size_t ew = eweight[*e];
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
    return p / w;
}


//A single Monte Carlo Markov chain sweep
template <class Graph, class BGraph, class EMprop, class Eprop, class Vprop,
          class EMat, class EVprop, class VEprop, class SamplerMap, class RNG>
void move_sweep(EMprop mrs, Vprop mrp, Vprop mrm, Vprop wr, Vprop b,
                Vprop clabel, vector<int>& vlist, bool deg_corr, bool dense,
                bool multigraph, double beta, Eprop eweight, Vprop vweight,
                EVprop egroups, VEprop esrcpos, VEprop etgtpos, Graph& g,
                BGraph& bg, EMat& emat, SamplerMap neighbour_sampler,
                bool sequential, bool random_move, double c, bool verbose,
                RNG& rng, double& S, size_t& nmoves)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

    size_t B = num_vertices(bg);

    typedef std::uniform_real_distribution<> rdist_t;
    auto rand_real = std::bind(rdist_t(), std::ref(rng));

    // it is useful to shuffle the vertex order even in the parallel case, so
    // threads become more balanced
    std::shuffle(vlist.begin(), vlist.end(), rng);

    std::uniform_int_distribution<size_t> s_rand(0, B - 1);

    nmoves = 0;
    S = 0;

    EntrySet<Graph> m_entries(B);

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
            std::uniform_int_distribution<size_t> v_rand(0, N - 1);
            v = vertex(vlist[v_rand(rng)], g);
        }

        vertex_t r = b[v];

        // blocks can't become empty
        if (wr[r] == 1)
            continue;

        // attempt random block
        vertex_t s = s_rand(rng);

        if (!random_move && total_degreeS()(v, g) > 0)
        {
            vertex_t u = neighbour_sampler[v].sample(rng);
            vertex_t t = b[u];

            double p_rand = 0;
            if (c > 0)
            {
                if (is_directed::apply<Graph>::type::value)
                    p_rand = c * B / double(mrp[t] + mrm[t] + c * B);
                else
                    p_rand = c * B / double(mrp[t] + c * B);
            }

            double sample_real = rand_real();

            if (sample_real >= p_rand)
            {
                const auto& e = egroups_manage::sample_edge(egroups[t], rng);
                s = b[target(e, g)];
                if (s == t)
                    s = b[source(e, g)];
            }
        }

        if (s == r)
            continue;

        if (clabel[s] != clabel[r])
            continue;

        double dS;
        if (!dense)
        {
            dS = virtual_move(v, s, mrs, mrp, mrm, wr, b, deg_corr, eweight,
                              vweight, g, bg, emat, m_entries);
        }
        else
        {
            if (!std::isinf(beta) && !random_move)
            {
                m_entries.Clear();
                move_entries(v, s, b, eweight, g, bg, m_entries);
            }

            dS = virtual_move_dense(v, s, mrs, mrp, mrm, wr, b, deg_corr,
                                    eweight, vweight, g, bg, emat, multigraph);
        }

        bool accept = false;
        if (std::isinf(beta))
        {
            accept = dS < 0;
        }
        else
        {
            double pf = random_move ? 1 :
                get_move_prob(v, r, s, c, b, mrs, mrp, mrm, emat, eweight, g,
                              bg, m_entries, false);

            double pb = random_move ? 1 :
                get_move_prob(v, s, r, c, b, mrs, mrp, mrm, emat, eweight, g,
                              bg, m_entries, true);

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

        if (!dense || (!std::isinf(beta) && !random_move))
            m_entries.Clear();

        if (accept)
        {
            move_vertex(v, s, mrs, mrp, mrm, wr, b, deg_corr, eweight,
                        vweight, g, bg, emat);
            if (!random_move)
                egroups_manage::update_egroups(v, r, s, eweight, egroups, esrcpos, etgtpos, g);
            S += dS;
            ++nmoves;
            if (verbose)
                cout << v << ": " << r << " -> " << s << " " << S << " " << vlist.size() << endl;
        }
    }
}


template <class Vertex, class Graph, class Eprop, class SMap>
void build_neighbour_sampler(Vertex v, SMap sampler, Eprop eweight, Graph& g)
{
    vector<Vertex> neighbours;
    vector<double> probs;
    neighbours.reserve(total_degreeS()(v, g));
    probs.reserve(total_degreeS()(v, g));

    typename all_edges_iteratorS<Graph>::type e, e_end;
    for(tie(e, e_end) = all_edges_iteratorS<Graph>::get_edges(v, g);
        e != e_end; ++e)
    {
        Vertex u = target(*e, g);
        if (is_directed::apply<Graph>::type::value && u == v)
            u = source(*e, g);
        neighbours.push_back(u);
        probs.push_back(eweight[*e]); // Self-loops will be added twice, and
                                      // hence will be sampled with probability
                                      // 2 * eweight[e]
    }

    sampler[v] = Sampler<Vertex, mpl::false_>(neighbours, probs);
};

struct init_neighbour_sampler
{
    template <class Graph, class Eprop>
    void operator()(Graph& g, Eprop eweight, boost::any& asampler) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

        typedef typename property_map<Graph, vertex_index_t>::type vindex_map_t;
        typedef typename property_map_type::apply<Sampler<vertex_t, mpl::false_>,
                                                  vindex_map_t>::type::unchecked_t
            sampler_map_t;

        sampler_map_t sampler(get(vertex_index_t(), g), num_vertices(g));
        asampler = sampler;

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) schedule(static)
        for (i = 0; i < N; ++i)
        {
            vertex_t v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            build_neighbour_sampler(v, sampler, eweight, g);
        }
    }
};


struct merge_cmp_less
{
    template<class Vertex>
    double operator()(const std::tuple<Vertex, Vertex, double>& a,
                      const std::tuple<Vertex, Vertex, double>& b)
    {
        return get<2>(a) < get<2>(b);
    }
};

struct match_cmp_less
{
    template<class Vertex>
    double operator()(const pair<Vertex, double>& a,
                      const pair<Vertex, double>& b)
    {
        return a.second < b.second;
    }
};


// merge block r into s
template <class BGraph, class Eprop, class Vprop, class EMat>
void merge_blocks(size_t r, size_t s, Eprop& mrs, Vprop& mrp, Vprop& mrm,
                  Vprop& wr, BGraph& bg, EMat& emat)
{
    typedef typename graph_traits<BGraph>::vertex_descriptor vertex_t;
    bool self_count = false;
    typename graph_traits<BGraph>::out_edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = out_edges(r, bg); ei != ei_end; ++ei)
    {
        typename graph_traits<BGraph>::edge_descriptor e = *ei;
        vertex_t t = target(e, bg);
        if (t == r && !is_directed::apply<BGraph>::type::value)
        {
            if (self_count)
                continue;
            self_count = true;
        }

        remove_me(r, t, e, emat, bg, false);

        if (t == r)
            t = s;

        pair<typename graph_traits<BGraph>::edge_descriptor, bool> ne =
            get_me(s, t, emat, bg);

        if (!ne.second)
        {
            ne = add_edge(s, t, bg);
            put_me(s, t, ne.first, emat, bg);
        }
        mrs[ne.first] += mrs[e];
        mrs[e] = 0;
    }

    typename in_edge_iteratorS<BGraph>::type ie, ie_end;
    for (tie(ie, ie_end) = in_edge_iteratorS<BGraph>::get_edges(r, bg);
         ie != ie_end; ++ie)
    {
        typename graph_traits<BGraph>::edge_descriptor e = *ie;
        vertex_t t = source(e, bg);
        if (t == r)
            continue;

        remove_me(r, t, e, emat, bg, false);

        pair<typename graph_traits<BGraph>::edge_descriptor, bool> ne =
            get_me(t, s, emat, bg);

        if (!ne.second)
        {
            ne = add_edge(t, s, bg);
            put_me(t, s, ne.first, emat, bg);
        }
        mrs[ne.first] += mrs[e];
        mrs[e] = 0;
    }

    mrp[s] += mrp[r];
    if (is_directed::apply<BGraph>::type::value)
        mrm[s] += mrm[r];
    wr[s] += wr[r];
    wr[r] = 0;
    clear_vertex(r, bg);
}


template <class BGraph, class Eprop, class Vprop, class EMat, class SamplerMap,
          class RNG>
void merge_sweep(Eprop mrs, Vprop mrp, Vprop mrm, Vprop wr, Vprop b,
                 Vprop clabel, bool deg_corr, bool dense, bool multigraph,
                 BGraph& bg, EMat& emat, SamplerMap neighbour_sampler,
                 size_t nmerges, size_t nsweeps, bool random_moves,
                 bool verbose, RNG& rng, double& S, size_t& nmoves)
{
    typedef typename graph_traits<BGraph>::vertex_descriptor vertex_t;

    size_t B = num_vertices(bg);
    auto_ptr<EntrySet<BGraph> > m_entries;

    vector<set<pair<vertex_t, double>, match_cmp_less> > past_moves;
    vector<pair<vertex_t, double> > best_move;
    past_moves.resize(B);
    best_move.resize(B, make_pair(vertex_t(0), numeric_limits<double>::max()));

    nmoves = 0;
    size_t n_attempts = 0;
    while (n_attempts < nsweeps)
    {
        n_attempts++;

        EntrySet<BGraph> m_entries(B);
        int i = 0, N = B;
        for (i = 0; i < N; ++i)
        {
            vertex_t r = vertex(i, bg);

            if (b[r] != r)
                continue;

            vertex_t s;
            if (random_moves || total_degreeS()(r, bg) == 0)
            {
                std::uniform_int_distribution<size_t> srand(0, B - 1);
                s = srand(rng);
                while (b[s] != s)
                    s = srand(rng);
            }
            else
            {
                vertex_t t = neighbour_sampler[r].sample(rng);
                s = neighbour_sampler[t].sample(rng);
            }

            if (s == r)
                continue;

            if (clabel[r] != clabel[s])
                continue;

            double dS = 0;

            typeof(past_moves[r]).begin() iter = past_moves[r].find(make_pair(s, dS));
            if (iter == past_moves[r].end())
            {
                if (!dense)
                {
                    dS = virtual_move(r, s, mrs, mrp, mrm, wr, b, deg_corr, mrs,
                                      wr, bg, bg, emat, m_entries);
                    m_entries.Clear();
                }
                else
                {
                    dS = virtual_move_dense(r, s, mrs, mrp, mrm, wr, b,
                                            deg_corr, mrs, wr, bg, bg, emat,
                                            multigraph);
                }
                past_moves[r].insert(make_pair(s, dS));
            }
            else
            {
                dS = iter->second;
            }

            if (dS < best_move[r].second)
            {
                best_move[r].first = s;
                best_move[r].second = dS;
            }

        }
    }

    if (nmerges > 0)
    {
        EntrySet<BGraph> m_entries(B);

        // top is the merge with _largest_ dS
        priority_queue<std::tuple<vertex_t, vertex_t, double>,
                       vector<std::tuple<vertex_t, vertex_t, double> >,
                       merge_cmp_less> merge_heap;

        for (size_t i = 0; i < B; ++i)
        {
            vertex_t r = vertex(i, bg);
            vertex_t s = best_move[r].first;
            double dS = best_move[r].second;

            if (r != s && dS < numeric_limits<double>::max() &&
                (merge_heap.size() < nmerges || dS < get<2>(merge_heap.top())))
            {
                merge_heap.push(std::make_tuple(r, s, dS));
                if (merge_heap.size() > nmerges)
                    merge_heap.pop();
            }
        }

        // back is the merge with _smallest_ dS
        vector<pair<vertex_t, vertex_t> > best_merges;
        best_merges.reserve(best_merges.size());

        while (!merge_heap.empty())
        {
            std::tuple<vertex_t, vertex_t, double> v = merge_heap.top();
            best_merges.push_back(make_pair(get<0>(v), get<1>(v)));
            merge_heap.pop();
        }

        vector<bool> touched(B, false);

        while (!best_merges.empty())
        {
            vertex_t r = best_merges.back().first;
            vertex_t s = best_merges.back().second;
            best_merges.pop_back();

            double dS;
            if (!dense)
            {
                dS = virtual_move(r, s, mrs, mrp, mrm, wr, b, deg_corr, mrs,
                                  wr, bg, bg, emat, m_entries);
                m_entries.Clear();
            }
            else
            {
                dS = virtual_move_dense(r, s, mrs, mrp, mrm, wr, b,
                                        deg_corr, mrs, wr, bg, bg, emat,
                                        multigraph);
            }

            if (touched[r] || touched[s])
                continue;

            touched[r] = touched[s] = true;

            typename Eprop::checked_t ew = mrs.get_checked();

            merge_blocks(r, s, ew, mrp, mrm, wr, bg, emat);
            b[r] = s;

            ++nmoves;
            S += dS;
        }

        // collapse merge tree across multiple calls
        int i = 0, N = B;
        for (i = 0; i < N; ++i)
        {
            vertex_t r = vertex(i, bg);
            vertex_t s = r;
            while (b[s] != s)
                s = b[s];
            b[r] = s;
        }

        #pragma omp parallel for default(shared) private(i) schedule(static) if (N > 100)
        for (i = 0; i < N; ++i)
        {
            vertex_t r = vertex(i, bg);
            if (total_degreeS()(r, bg) > 0)
                build_neighbour_sampler(r, neighbour_sampler, mrs, bg);
            else
                neighbour_sampler[r] = Sampler<vertex_t, mpl::false_>();
        }
    }
}


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
