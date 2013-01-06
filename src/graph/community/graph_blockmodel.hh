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

#ifndef GRAPH_BLOCKMODEL_HH
#define GRAPH_BLOCKMODEL_HH

#include <cmath>
#include <iostream>

#include <ext/numeric>
using __gnu_cxx::power;

#  ifndef __clang__
#   include <tr1/unordered_set>
#   include <tr1/tuple>
#  else
#   include <boost/tr1/unordered_set.hpp>
#   include <boost/tr1/tuple.hpp>
#  endif

#include <dense_hash_set>
#include <dense_hash_map>

namespace graph_tool
{

using google::dense_hash_set;
using google::dense_hash_map;

using namespace boost;

// ====================
// Entropy calculation
// ====================

template <class Type>
__attribute__((always_inline))
inline double safelog(Type x)
{
    if (x == 0)
        return 0;
    return log(x);
}

// "edge" term of the entropy
template <class Graph, class Edge, class Eprop>
__attribute__((always_inline))
inline double eterm(const Edge& e, const Eprop& mrs, const Graph& g)
{
    typename graph_traits<Graph>::vertex_descriptor r, s;
    r = source(e, g);
    s = target(e, g);

    const size_t mrse = mrs[e];

    // Repeated computation of x*log(x) actually adds up to a lot of time. A
    // significant speedup can be made by caching pre-computed values. This
    // is doable since the values of mrse are bounded in [0, 2E], where E is
    // the total number of edges in the network. Here we simply grow the
    // cache as needed.

    static vector<double> cache;  // FIXME: This will persist throughout the
                                  // process' life. It would be better to
                                  // replace with a version which will be
                                  // cleaned up at some point.
    {
        #pragma omp critical
        if (mrse >= cache.size())
        {
            size_t old_size = cache.size();
            cache.resize(mrse + 1);
            for (size_t i = old_size; i < cache.size(); ++i)
                cache[i] = i * safelog(i);
        }
    }

    double val = cache[mrse];

    if (is_directed::apply<Graph>::type::value || r != s)
        return -val;
    else
        return -val / 2;
}

// "vertex" term of the entropy
template <class Graph, class Vertex, class Vprop>
inline double vterm(Vertex v, Vprop& mrp, Vprop& mrm, Vprop& wr, bool deg_corr,
                    Graph& g)
{
    double one = 0.5;

    if (is_directed::apply<Graph>::type::value)
        one = 1;

    if (deg_corr)
        return one * (mrm[v] * safelog(mrm[v]) + mrp[v] * safelog(mrp[v]));
    else
        return one * (mrm[v] * safelog(wr[v]) + mrp[v] * safelog(wr[v]));
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
            S += eterm(*e, mrs, g);
        typename graph_traits<Graph>::vertex_iterator v, v_end;
        for (tie(v, v_end) = vertices(g); v != v_end; ++v)
            S += vterm(*v, mrp, mrm, wr, deg_corr, g);
    }
};

// ===============================
// Block merge and merge distance
// ===============================


// obtain the "super block" which is a merger of two blocks
template <class Graph, class Vertex, class Eprop, class Vprop>
void super_block(Vertex u, Vertex v, Eprop mrs, Vprop mrp, Vprop mrm, Vprop wr,
                 bool deg_corr, Graph& g, unordered_set<Vertex>& nsp,
                 unordered_set<Vertex>& nsm, unordered_map<Vertex, size_t>& msp,
                 unordered_map<Vertex, size_t>& msm, size_t& ml, size_t& m_rp,
                 size_t& m_rm, size_t& w_r)
{
    ml = 0;
    size_t ml2 = 0;

    Vertex ws[2];
    ws[0] = u;
    ws[1] = v;

    for (size_t i = 0; i < 2; ++i)
    {
        Vertex w = ws[i];
        typename all_edges_iteratorS<Graph>::type e, e_end;
        for (tie(e, e_end) = all_edges_iteratorS<Graph>::get_edges(w, g);
             e != e_end; ++e)
        {
            Vertex t = target(*e, g);
            Vertex s = source(*e, g);
            if ((s == w && (t == u || t == v)) ||
                (t == w && (s == u || s == v)))
            {
                if (target(*e, g) == w || is_directed::apply<Graph>::type::value)
                    ml2 += mrs[*e];
                else
                    ml += mrs[*e];
            }
            else
            {
                if (s == w)
                {
                    nsp.insert(t);
                    msp[t] += mrs[*e];
                }
                else
                {
                    nsm.insert(s);
                    msm[s] += mrs[*e];
                }
            }
        }
    }
    ml += ml2 / 2;

    m_rp = mrp[u] + mrp[v];
    m_rm = mrm[u] + mrm[v];
    w_r = wr[u] + wr[v];
}

// compute the "distance" between two blocks (i.e. entropy difference when merged)
struct dist
{
    template <class Graph, class Vertex, class Eprop, class Vprop>
    void operator()(Vertex r, Vertex s, const Eprop& mrs, const Vprop& mrp,
                    const Vprop& mrm, const Vprop& wr, bool deg_corr, Graph& g,
                    double& d) const
    {
        if (r == s)
        {
            d = 0;
            return;
        }

        d = vterm(r, mrp, mrm, wr, deg_corr, g) +
            vterm(s, mrp, mrm, wr, deg_corr, g);

        typename all_edges_iteratorS<Graph>::type e, e_end;
        for (tie(e, e_end) = all_edges_iteratorS<Graph>::get_edges(r, g);
             e != e_end; ++e)
        {
            d += eterm(*e, mrs, g);
        }

        for (tie(e, e_end) = all_edges_iteratorS<Graph>::get_edges(s, g);
             e != e_end; ++e)
        {
            if ((target(*e, g) == s && source(*e, g) == r) ||
                (target(*e, g) == r && source(*e, g) == s))
                continue;
            d += eterm(*e, mrs, g);
        }

        unordered_set<Vertex> nsp, nsm;
        unordered_map<Vertex, size_t> msp, msm;
        size_t ml, m_rp, m_rm, w_r;

        super_block(r, s, mrs, mrp, mrm, wr, deg_corr, g, nsp, nsm, msp, msm,
                    ml, m_rp, m_rm, w_r);

        if (!is_directed::apply<Graph>::type::value)
        {
            msm = msp;
        }


        for (typeof(nsp.begin()) iter = nsp.begin(); iter != nsp.end(); ++iter)
        {
            Vertex t = *iter;
            if (deg_corr)
                d += msp[t] * (safelog(msp[t]) - safelog(m_rp) - safelog(mrm[t]));
            else
                d += msp[t] * (safelog(msp[t]) - safelog(w_r) - safelog(wr[t]));
        }

        for (typeof(nsm.begin()) iter = nsm.begin(); iter != nsm.end(); ++iter)
        {
            Vertex t = *iter;
            if (deg_corr)
                d += msm[t] * (safelog(msm[t]) - safelog(m_rm) - safelog(mrp[t]));
            else
                d += msm[t] * (safelog(msm[t]) - safelog(w_r) - safelog(wr[t]));
        }

        if (ml > 0)
        {
            double one = 0.5;
            if (is_directed::apply<Graph>::type::value)
                one = 1;

            if (deg_corr)
                d += one * ml * (safelog(ml) - safelog(m_rp) - safelog(m_rm));
            else
                d += one * ml * (safelog(ml) - safelog(w_r) - safelog(w_r));
        }

        d *= -1;
    }
};

// obtain the minimum distance between all block pairs
template <class RNG>
struct min_dist
{
    min_dist(RNG& rng): rng(rng) {}
    RNG& rng;

    template <class Graph, class Eprop, class Vprop>
    void operator()(int n, Eprop mrs, Vprop mrp, Vprop mrm, Vprop wr,
                    bool deg_corr, Graph& g, vector<int>& vlist,
                    boost::tuple<double, size_t, size_t>& min_d) const
    {
        get<0>(min_d) = numeric_limits<double>::max();
        std::tr1::uniform_int<size_t> rand(0, vlist.size() - 1);

        bool random = true;
        if (n <= 0)
        {
            n = vlist.size() * vlist.size();
            random = false;
        }

        int l;
        //#pragma omp parallel for default(shared) private(l)
        for (l = 0; l < n; ++l)
        {
            typename graph_traits<Graph>::vertex_descriptor r;
            if (random)
            {
                #pragma omp critical
                r = vertex(vlist[rand(rng)], g);
            }
            else
            {
                r = vlist[l / vlist.size()];
            }

            if (wr[r] == 0)
                continue;

            typename graph_traits<Graph>::vertex_descriptor s;
            if (random)
            {
                do
                {
                    #pragma omp critical
                    s = vertex(vlist[rand(rng)], g);
                }
                while (s == r);
            }
            else
            {
                s = vlist[l % vlist.size()];
                if (s <= r)
                    continue;
            }

            if (wr[s] == 0)
                continue;

            if (s != r)
            {
                double d = 0;
                dist()(r, s, mrs, mrp, mrm, wr, deg_corr, g, d);

                {
                    #pragma omp critical
                    if (d < get<0>(min_d))
                    {
                        get<0>(min_d) = d;
                        get<1>(min_d) = r;
                        get<2>(min_d) = s;
                    }
                }
            }
        }
    }
};

// merge two blocks into one
struct b_join
{
    b_join(size_t r, size_t s): rr(r), ss(s) {}
    size_t rr, ss;

    template <class Graph, class Eprop, class Vprop>
    void operator()(Eprop mrs, Vprop mrp, Vprop mrm, Vprop wr, bool deg_corr,
                    vector<int>& vlist, Graph& g) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        unordered_set<vertex_t> nsp, nsm;
        unordered_map<vertex_t, size_t> msp, msm;
        size_t ml, m_rp, m_rm, w_r;

        vertex_t r = vertex(rr, g);
        vertex_t s = vertex(ss, g);

        super_block(r, s, mrs, mrp, mrm, wr, deg_corr, g, nsp, nsm, msp, msm,
                    ml, m_rp, m_rm, w_r);

        mrp[r] = m_rp;
        mrm[r] = m_rm;
        wr[r] = w_r;

        clear_vertex(r, g);
        clear_vertex(s, g);

        for (typeof(nsp.begin()) iter = nsp.begin(); iter != nsp.end(); ++iter)
        {
            vertex_t t = *iter;
            typename graph_traits<Graph>::edge_descriptor e = add_edge(r, t, g).first;
            mrs[e] = msp[t];
        };

        for (typeof(nsm.begin()) iter = nsm.begin(); iter != nsm.end(); ++iter)
        {
            vertex_t s = *iter;
            typename graph_traits<Graph>::edge_descriptor e = add_edge(s, r, g).first;
            mrs[e] = msm[s];
        };

        if (ml > 0)
        {
            typename graph_traits<Graph>::edge_descriptor e = add_edge(r, r, g).first;
            mrs[e] = ml;
        }

        mrp[s] = 0;
        mrm[s] = 0;
        wr[s] = 0;

    }
};

// ===============================
// Block moves
// ===============================


// this structure speeds up the access to the edges between to given blocks,
// since we're using and adjacency list to store the block structure (the emat_t
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
    void operator()(Graph& g, boost::any& oemap, size_t B) const
    {
        typedef typename get_emat_t::apply<Graph>::type emat_t;
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

// remove a vertex from its current block
template <class Graph, class BGraph, class Eprop, class Vprop, class EWprop,
          class VWprop, class EMat>
void remove_vertex(size_t v, Eprop& mrs, Vprop& mrp, Vprop& mrm, Vprop& wr,
                   Vprop& b, const EWprop& eweight, const VWprop& vweight,
                   Graph& g, BGraph& bg, EMat& emat)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    vertex_t r = b[v];

    typename graph_traits<Graph>::out_edge_iterator e, e_end;
    for (tie(e, e_end) = out_edges(vertex(v, g), g); e != e_end; ++e)
    {
        vertex_t u = target(*e, g);
        vertex_t s = b[u];

        if (!emat[r][s].second)
            throw GraphException("no edge? " + lexical_cast<string>(r) +
                                 " " + lexical_cast<string>(s));

        const typename graph_traits<BGraph>::edge_descriptor& me = emat[r][s].first;

        size_t ew = eweight[*e];
        mrs[me] -= ew;

        if (u != v && r == s && !is_directed::apply<Graph>::type::value)
            mrs[me] -= ew;

        mrp[r] -= ew;
        if (u != v || is_directed::apply<Graph>::type::value)
            mrm[s] -= ew;
    }

    if (is_directed::apply<Graph>::type::value)
    {
        typename in_edge_iteratorS<Graph>::type ie, ie_end;
        for (tie(ie, ie_end) = in_edge_iteratorS<Graph>::get_edges(vertex(v, g), g);
             ie != ie_end; ++ie)
        {
            vertex_t u = source(*ie, g);
            if (u == v)
                continue;
            vertex_t s = b[u];

            if (!emat[s][r].second)
                throw GraphException("no edge? " + lexical_cast<string>(s) +
                                     " " + lexical_cast<string>(r));

            typename graph_traits<BGraph>::edge_descriptor me = emat[s][r].first;

            size_t ew = eweight[*ie];
            mrs[me] -= ew;

            mrp[s] -= ew;
            mrm[r] -= ew;
        }
    }

    wr[vertex(r, bg)] -= vweight[vertex(v, g)];

}

// add a vertex to block rr
template <class Graph, class BGraph, class Eprop, class Vprop, class EWprop,
          class VWprop,class EMat>
void add_vertex(size_t v, size_t rr, Eprop& mrs, Vprop& mrp, Vprop& mrm,
                Vprop& wr, Vprop& b, const EWprop& eweight,
                const VWprop& vweight, Graph& g, BGraph& bg, EMat& emat)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    vertex_t r = vertex(rr, bg);

    typename graph_traits<Graph>::out_edge_iterator e, e_end;
    for (tie(e, e_end) = out_edges(vertex(v, g), g); e != e_end; ++e)
    {
        vertex_t u = target(*e, g);
        vertex_t s;

        if (u != v)
            s = b[u];
        else
            s = r;

        typename graph_traits<BGraph>::edge_descriptor me;

        pair<typename graph_traits<BGraph>::edge_descriptor, bool>& mep = emat[r][s];

        if (!mep.second)
        {
            emat[r][s] = add_edge(r, s, bg);
            if (!is_directed::apply<Graph>::type::value)
                emat[s][r] = emat[r][s];
            me = emat[r][s].first;
            mrs[me] = 0;
        }
        else
        {
            me = mep.first;
        }

        size_t ew = eweight[*e];

        mrs[me] += ew;

        if (u != v && r == s && !is_directed::apply<Graph>::type::value)
            mrs[me] += ew;

        mrp[r] += eweight[*e];
        if (u != v || is_directed::apply<Graph>::type::value)
            mrm[s] += ew;
    }

    typename in_edge_iteratorS<Graph>::type ie, ie_end;
    for (tie(ie, ie_end) = in_edge_iteratorS<Graph>::get_edges(vertex(v, g), g);
         ie != ie_end; ++ie)
    {
        vertex_t u = source(*ie, g);
        if (u == v)
            continue;

        vertex_t s = b[u];

        typename graph_traits<BGraph>::edge_descriptor me;
        pair<typename graph_traits<BGraph>::edge_descriptor, bool>& mep = emat[s][r];

        if (!mep.second)
        {
            emat[s][r] = add_edge(s, r, bg);
            me = emat[s][r].first;
            mrs[me] = 0;
        }
        else
        {
            me = mep.first;
        }

        size_t ew = eweight[*ie];

        mrs[me] += ew;

        mrp[s] += ew;
        mrm[r] += ew;
    }

    wr[vertex(r, bg)] += vweight[vertex(v, g)];
    b[v] = r;
}

template <class Type1, class Type2, class Graph>
inline pair<Type1,Type2> make_ordered_pair(const Type1& v1, const Type2& v2, const Graph&)
{
    if (!is_directed::apply<Graph>::type::value)
        return make_pair(min(v1, v2), max(v1, v2));
     return make_pair(v1, v2);
}

template <class Vertex, class MEntry, class MEntrySet, class Graph>
inline void insert_m_entry(Vertex r, Vertex s, MEntry& m_entries,
                           MEntrySet& m_entries_set, Graph& g)
{
    pair<Vertex, Vertex> rs = make_ordered_pair(r, s, g);
    if (m_entries_set.find(rs) == m_entries_set.end())
    {
        m_entries.push_back(rs);
        m_entries_set.insert(rs);
    }
}

// obtain the necessary entries in the e_rs matrix which need to be modified
// after the move
template <class Graph, class BGraph, class Vertex, class Vprop, class MEntry,
          class MEntrySet>
void move_entries(Vertex v, Vertex nr, Vprop& b, Graph& g, BGraph& bg,
                  MEntry& m_entries, MEntrySet& m_entries_set)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    vertex_t r = b[v];

    typename graph_traits<Graph>::out_edge_iterator e, e_end;
    for (tie(e, e_end) = out_edges(vertex(v, g), g); e != e_end; ++e)
    {
        vertex_t s = b[target(*e, g)];
        insert_m_entry(r, s, m_entries, m_entries_set, g);
        insert_m_entry(nr, s, m_entries, m_entries_set, g);
        if (s == r || s == nr)
            insert_m_entry(nr, nr, m_entries, m_entries_set, g);
        if (s == nr)
            insert_m_entry(r, r, m_entries, m_entries_set, g);
    }

    typename in_edge_iteratorS<Graph>::type ie, ie_end;
    for (tie(ie, ie_end) = in_edge_iteratorS<Graph>::get_edges(vertex(v, g), g);
         ie != ie_end; ++ie)
    {
        vertex_t s = b[source(*ie, g)];

        insert_m_entry(s, r, m_entries, m_entries_set, g);
        insert_m_entry(s, nr, m_entries, m_entries_set, g);
        if (s == r || s == nr)
            insert_m_entry(nr, nr, m_entries, m_entries_set, g);
        if (s == nr)
            insert_m_entry(r, r, m_entries, m_entries_set, g);
    }
}


// obtain the entropy difference given a set of entries in the e_rs matrix
template <class MEntry, class Eprop, class BGraph, class EMat>
double entries_dS(MEntry& m_entries, Eprop& mrs, BGraph& bg, const EMat& emat)
{
    typedef typename graph_traits<BGraph>::vertex_descriptor vertex_t;
    double dS = 0;
    for (typeof(m_entries.begin()) iter = m_entries.begin(); iter != m_entries.end(); ++iter)
    {
        vertex_t er = iter->first;
        vertex_t es = iter->second;

        const pair<typename graph_traits<BGraph>::edge_descriptor, bool>& e = emat[er][es];
        if (!e.second)
            continue;
        dS += eterm(e.first, mrs, bg);
    }
    return dS;
}


// move a vertex from its current block to block nr, and return the entropy
// difference
template <class Graph, class BGraph, class Eprop, class Vprop, class EWprop,
          class VWprop, class EMat>
double move_vertex(size_t v, size_t nr, Eprop& mrs, Vprop& mrp, Vprop& mrm,
                   Vprop& wr, Vprop& b, bool deg_corr, const EWprop& eweight,
                   const VWprop& vweight, Graph& g, BGraph& bg, EMat& emat,
                   bool get_ds=true)

{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    vertex_t r = b[vertex(v, g)];

    if (r == nr)
        return 0;

    double Si = 0, Sf = 0;

    vector<pair<vertex_t, vertex_t> > m_entries;
    dense_hash_set<pair<vertex_t, vertex_t>, boost::hash<pair<vertex_t, vertex_t> > > m_entries_set;
    if (get_ds)
    {
        m_entries_set.set_empty_key(make_pair(graph_traits<Graph>::null_vertex(),
                                              graph_traits<Graph>::null_vertex()));

        move_entries(vertex(v, g), nr, b, g, bg, m_entries, m_entries_set);

        Si = entries_dS(m_entries, mrs, bg, emat);
        Si += vterm(r, mrp, mrm, wr, deg_corr, bg);
        Si += vterm(nr, mrp, mrm, wr, deg_corr, bg);
    }

    remove_vertex(v, mrs, mrp, mrm, wr, b, eweight, vweight, g, bg, emat);
    add_vertex(v, nr, mrs, mrp, mrm, wr, b, eweight, vweight, g, bg, emat);

    if (get_ds)
    {
        Sf = entries_dS(m_entries, mrs, bg, emat);
        Sf += vterm(r, mrp, mrm, wr, deg_corr, bg);
        Sf += vterm(nr, mrp, mrm, wr, deg_corr, bg);
    }
    return Sf - Si;
}
//============
// Main loops
//============


template <class Graph, class RNG>
typename graph_traits<Graph>::vertex_descriptor
random_neighbour(typename graph_traits<Graph>::vertex_descriptor v, Graph& g,
                 RNG& rng)
{
    if (total_degreeS()(v, g) == 0)
        return v;
    std::tr1::uniform_int<size_t> rand(0, total_degreeS()(v, g) - 1);

    // if (!is_convertible<typename iterator_traits<typename all_edges_iteratorS<Graph>::type>::iterator_category,
    //                     random_access_iterator_tag>::value)
    //     throw GraphException("not random access...");

    size_t i = rand(rng);
    typename all_edges_iteratorS<Graph>::type e, e_end;
    tie(e, e_end) = all_edges_iteratorS<Graph>::get_edges(v, g);
    std::advance(e, i);
    return target(*e, g);
}

template <class Vertex, class Eprop, class Emat>
inline size_t get_mrs(Vertex r, Vertex s, const Eprop& mrs, const Emat& emat)
{
    const typename Emat::element& me = emat[r][s];
    if (me.second)
        return mrs[me.first];
    else
        return 0;
}

//the following guarantee a stable (source, target) ordering even for undirected
//graphs
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

//computes the move proposal probability
template <class Vertex, class Graph, class Vprop, class Eprop, class Emat>
inline double
get_move_prob(Vertex v, Vertex s, Vprop b, Eprop mrs, Vprop mrp, Vprop mrm,
              Emat& emat, Graph& g, size_t B)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    double p = 0;
    typename all_edges_iteratorS<Graph>::type e, e_end;
    for (tie(e, e_end) = all_edges_iteratorS<Graph>::get_edges(v, g);
         e != e_end; ++e)
    {
        vertex_t t = b[target(*e, g)];
        if (is_directed::apply<Graph>::type::value)
            p += (get_mrs(t, s, mrs, emat) + get_mrs(s, t, mrs, emat) + 1.) / (mrp[t] + mrm[t] + B);
        else
            p += (get_mrs(t, s, mrs, emat) + 1.) / (mrp[t] + B);
    }
    return p / total_degreeS()(v, g);
}

//A single Monte Carlo Markov chain sweep
template <class Graph, class BGraph, class Eprop, class Vprop, class EMat,
          class EVprop, class VEprop, class RNG>
void move_sweep(Eprop mrs, Vprop mrp, Vprop mrm, Vprop wr, Vprop b, Vprop label,
                size_t L, vector<int>& vlist, bool deg_corr, double beta,
                Eprop eweight, Vprop vweight, EVprop egroups, VEprop esrcpos,
                VEprop etgtpos, Graph& g, BGraph& bg, EMat& emat,
                bool sequential, bool verbose, RNG& rng, double& S,
                size_t& nmoves)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    nmoves = 0;


    typedef std::tr1::uniform_real<> rdist_t;
    std::tr1::variate_generator<RNG&, rdist_t> rand_real(rng, rdist_t());

    typedef random_permutation_iterator<typename vector<int>::iterator,
                                        RNG> random_vertex_iter;
    random_vertex_iter viter(vlist.begin(), vlist.end(), rng),
        vi_end(vlist.end(), vlist.end(), rng);
    for (; viter != vi_end; ++viter)
    {
        std::tr1::uniform_int<size_t> rand(0, num_vertices(bg) / L - 1);

        vertex_t v = vertex(*viter, g);
        if (!sequential)
        {
            std::tr1::uniform_int<size_t> randr(0, vlist.size() - 1);
            v = vertex(vlist[randr(rng)], g);
        }

        vertex_t r = b[v];


        vertex_t u = random_neighbour(v, g, rng);
        if (u == v)
            continue;
        vertex_t t = b[u];

        // attempt random block first
        vertex_t s = vertex(rand(rng) * L + label[v], bg);
        size_t B = num_vertices(bg);
        double p_rand;
        if (is_directed::apply<Graph>::type::value)
            p_rand = B / double(mrp[t] + mrm[t] + B);
        else
            p_rand = B / double(mrp[t] + B);

        if (rand_real() >= p_rand)
        {
            std::tr1::uniform_int<size_t> urand(0, egroups[t].size() - 1);

            const typename graph_traits<Graph>::edge_descriptor& e = get<0>(egroups[t][urand(rng)]);
            s = b[target(e, g)];
            if (s == t)
                s = b[source(e, g)];

            if (int(s % L) != label[v])
                continue;
        }

        if (s == r)
            continue;

        double pf = isinf(beta) ? 1 : get_move_prob(v, s, b, mrs, mrp, mrm, emat, g, B);

        double dS = move_vertex(v, s, mrs, mrp, mrm, wr, b, deg_corr, eweight,
                                vweight, g, bg, emat, true);

        double pb = isinf(beta) ? 1 : get_move_prob(v, r, b, mrs, mrp, mrm, emat, g, B);

        double a = isinf(beta) ? -dS : -beta * dS + log(pb) - log(pf);

        bool accept = false;
        if (a > 0)
            accept = true;
        else if (!isinf(beta))
            accept = rand_real() < exp(a);

        if (!accept)
        {
            move_vertex(v, r, mrs, mrp, mrm, wr, b, deg_corr, eweight,
                        vweight, g, bg, emat, false);
        }
        else
        {
            S += dS;
            ++nmoves;

            size_t self_count = 0;

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
                    is_src = (self_count == 0);
                    ++self_count;
                }

                for (size_t i = 0; i < size_t(eweight[*e]); ++i)
                {
                    size_t pos = (is_src) ? esrcpos[*e][i] : etgtpos[*e][i];

                    const boost::tuple<typename graph_traits<Graph>::edge_descriptor, bool, size_t>& eback = egroups[r].back();

                    vertex_t u = get<1>(eback) ? get_source(get<0>(eback), g) : get_target(get<0>(eback), g);

                    if (get<1>(eback))
                        esrcpos[get<0>(eback)][get<2>(eback)] = pos;
                    else
                        etgtpos[get<0>(eback)][get<2>(eback)] = pos;

                    egroups[r][pos] = eback;
                    egroups[r].pop_back();
                }
            }

            self_count = 0;
            for (tie(e, e_end) = all_edges_iteratorS<Graph>::get_edges(v, g);
                 e != e_end; ++e)
            {
                vertex_t src = get_source(*e, g);
                vertex_t tgt = get_target(*e, g);

                bool is_src = (src == v);

                // self-loops will appear twice
                if (src == tgt)
                {
                    is_src = (self_count == 0);
                    ++self_count;
                }

                for (size_t i = 0; i < size_t(eweight[*e]); ++i)
                {
                    egroups[s].push_back(boost::make_tuple(*e, is_src, i));
                    if (is_src)
                        esrcpos[*e][i] = egroups[s].size() - 1;
                    else
                        etgtpos[*e][i] = egroups[s].size() - 1;
                }
            }

            if (verbose)
                cout << v << ": " << r << " -> " << s << " " << S << " " << vlist.size() << endl;
        }
    }
}


// construct half-edge lists
struct build_egroups
{

    template <class Eprop, class Vprop, class VEprop, class Graph, class VertexIndex>
    void operator()(Vprop b, boost::any& oegroups, VEprop esrcpos, VEprop etgtpos,
                    Eprop eweight, Graph& g, VertexIndex vertex_index) const
    {
        typedef typename property_map_type::apply<vector<boost::tuple<typename graph_traits<Graph>::edge_descriptor, bool, size_t> >,
                                                  GraphInterface::vertex_index_map_t>::type vemap_t;
        vemap_t egroups(vertex_index);
        oegroups = egroups;

        typename graph_traits<Graph>::edge_iterator e, e_end;
        for (tie(e, e_end) = edges(g); e != e_end; ++e)
        {
            size_t r = b[get_source(*e, g)];
            esrcpos[*e].resize(eweight[*e]);
            for (size_t i = 0; i < size_t(eweight[*e]); ++i)
            {
                egroups[r].push_back(boost::make_tuple(*e, true, i));
                esrcpos[*e][i] = egroups[r].size() - 1;
            }
            size_t s = b[get_target(*e, g)];
            etgtpos[*e].resize(eweight[*e]);
            for (size_t i = 0; i < size_t(eweight[*e]); ++i)
            {
                egroups[s].push_back(boost::make_tuple(*e, false, i));
                etgtpos[*e][i] = egroups[s].size() - 1;
            }
        }
    }

};

// Sampling marginal probabilities on the edges
template <class Graph, class Vprop, class MEprop>
void collect_edge_marginals(size_t B, Vprop b, MEprop p, Graph& g, Graph& bg)
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
