// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2010 Tiago de Paula Peixoto <tiago@forked.de>
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

#ifndef SUBGRAPH_ISOMORPHISM_HPP
#define SUBGRAPH_ISOMORPHISM_HPP

#include <boost/graph/graph_traits.hpp>
#include <utility>
#include <tr1/unordered_set>
#include <tr1/random>

namespace boost
{
using namespace std;

typedef tr1::mt19937 rng_t;

namespace detail {

//sparse matrix
typedef vector<tr1::unordered_set<size_t> > matrix_t;

struct check_adjacency
{
    template <class Graph>
    typename graph_traits<Graph>::vertex_descriptor
    get_vertex(const typename graph_traits<Graph>::edge_descriptor& e, Graph& g,
               mpl::true_ out_edges)
    {
        return target(e, g);
    }

    template <class Graph>
    typename graph_traits<Graph>::vertex_descriptor
    get_vertex(const typename graph_traits<Graph>::edge_descriptor& e, Graph& g,
               mpl::false_ out_edges)
    {
        return source(e, g);
    }

    template <class Graph, class IsOut>
    struct get_edge_iterator
    {
        typedef typename graph_traits<Graph>::out_edge_iterator type;

        static pair<type,type>
        edges(typename graph_traits<Graph>::vertex_descriptor v, Graph& g)
        {
            return out_edges(v, g);
        }
    };

    template <class Graph>
    struct get_edge_iterator<Graph, mpl::false_>
    {
        typedef typename graph_traits<Graph>::in_edge_iterator type;

        static pair<type,type>
        edges(typename graph_traits<Graph>::vertex_descriptor v, Graph& g)
        {
            return in_edges(v, g);
        }
    };

    template <class Graph1, class Graph2, class EdgeLabelling>
    bool operator()(typename graph_traits<Graph1>::vertex_descriptor k,
                    typename graph_traits<Graph2>::vertex_descriptor l,
                    matrix_t& M, EdgeLabelling& edge_labelling,
                    Graph1& g1, Graph2& g2, vector<size_t>& vindex,
                    mpl::true_ directed)
    {
        return do_check(k, l, M, edge_labelling, g1, g2, vindex,
                        mpl::true_()) &&
            do_check(k, l, M, edge_labelling, g1, g2, vindex, mpl::false_());
    }

    template <class Graph1, class Graph2, class EdgeLabelling>
    bool operator()(typename graph_traits<Graph1>::vertex_descriptor k,
                    typename graph_traits<Graph2>::vertex_descriptor l,
                    matrix_t& M, EdgeLabelling& edge_labelling, Graph1& g1,
                    Graph2& g2, vector<size_t>& vindex, mpl::false_ directed)
    {
        return do_check(k, l, M, edge_labelling, g1, g2, vindex, mpl::true_());
    }

    template <class Graph1, class Graph2, class EdgeLabelling, class IsOut>
    bool do_check(typename graph_traits<Graph1>::vertex_descriptor k,
                  typename graph_traits<Graph2>::vertex_descriptor l,
                  matrix_t& M, EdgeLabelling& edge_labelling, Graph1& g1,
                  Graph2& g2, vector<size_t>& vindex, IsOut)
    {
        bool valid = true;
        typename get_edge_iterator<Graph1, IsOut>::type e1, e1_end;
        for (tie(e1, e1_end) =
                 get_edge_iterator<Graph1, IsOut>::edges(k, g1);
             e1 != e1_end; ++e1)
        {
            typename graph_traits<Graph1>::vertex_descriptor v1 =
                get_vertex(*e1, g1, IsOut());

            bool is_adjacent = false;
            typename get_edge_iterator<Graph2, IsOut>::type e2, e2_end;
            for (tie(e2, e2_end) =
                     get_edge_iterator<Graph2, IsOut>::edges(l, g2);
                 e2 != e2_end; ++e2)
            {
                typename graph_traits<Graph2>::vertex_descriptor v2 =
                    get_vertex(*e2, g2, IsOut());

                if (M[v1].find(vindex[v2]) != M[v1].end() &&
                    edge_labelling(*e1, *e2))
                {
                    is_adjacent = true;
                    break;
                }
            }

            if (!is_adjacent)
            {
                valid = false;
                break;
            }
        }
        return valid;
    }
};


template <class Graph1, class Graph2, class EdgeLabelling>
bool refine_check(const Graph1& sub, const Graph2& g, matrix_t& M, size_t count,
                  tr1::unordered_set<size_t>& already_mapped,
                  EdgeLabelling edge_labelling, vector<size_t>& vlist,
                  vector<size_t>& vindex)
{
    matrix_t M_temp(num_vertices(sub));

    int k = 0, N = num_vertices(sub);
    #pragma omp parallel for default(shared) private(k) schedule(dynamic)
    for (k = 0; k < int(count); ++k)
        M_temp[k] = M[k];

    size_t n_mod = 1;
    while (n_mod > 0)
    {
        n_mod = 0;
        bool abort = false;
        #pragma omp parallel for default(shared) private(k) schedule(dynamic) \
            reduction(+:n_mod)
        for (k = count; k < N; ++k)
        {
            if (abort)
                continue;
            if (vertex(k, sub) == graph_traits<Graph1>::null_vertex())
                continue;
            tr1::unordered_set<size_t> m_new;
            for (typeof(M[k].begin()) li = M[k].begin(); li != M[k].end(); ++li)
            {
                size_t l = *li;
                if (already_mapped.find(l) != already_mapped.end())
                    continue;
                bool valid = check_adjacency()
                    (vertex(k, sub), vertex(vlist[l], g), M, edge_labelling,
                     sub, g, vindex,
                     typename is_directed::apply<Graph1>::type());
                if (valid)
                    m_new.insert(l);
            }
            if (m_new.empty())
            {
                abort = true;
                continue;
            }
            M_temp[k].swap(m_new);
            if (M_temp[k].size() < M[k].size())
                n_mod++;
        }
        if (abort)
            return false;
        M.swap(M_temp);
    }
    return true;
}


template <class Graph1, class Graph2, class EdgeLabelling, class Mapping>
void find_mappings(const Graph1& sub, const Graph2& g, matrix_t& M0,
                   vector<Mapping>& FF, EdgeLabelling edge_labelling,
                   vector<size_t>& vlist, vector<size_t>& vindex, size_t max_n)
{
    size_t i = 0;
    for (i = 0; i < num_vertices(sub); ++i)
        if (vertex(i, sub) != graph_traits<Graph1>::null_vertex())
            break;
    int last_i = 0;
    for (last_i = num_vertices(sub) - 1; last_i >= 0; --last_i)
        if (vertex(i, sub) != graph_traits<Graph1>::null_vertex())
            break;
    for (; i < vlist.size(); ++i)
        if (vertex(vlist[i], g) != graph_traits<Graph2>::null_vertex())
            break;

    Mapping F;
    list<tuple<matrix_t, size_t,
               typename matrix_t::value_type::const_iterator> > Mstack;
    Mstack.push_back(make_tuple(M0, i, M0[i].begin()));
    get<2>(Mstack.back()) = get<0>(Mstack.back())[i].begin();
    tr1::unordered_set<size_t> already_mapped;

    // perform depth-first search of combination space
    while (!Mstack.empty() && (max_n == 0 || FF.size() < max_n))
    {
        const matrix_t& M = get<0>(Mstack.back());
        size_t& i = get<1>(Mstack.back());
        typename matrix_t::value_type::const_iterator& mi =
            get<2>(Mstack.back());

        if (mi == M[i].end())
        {
            // dead end
            Mstack.pop_back();
            if (!F.empty())
            {
                already_mapped.erase(F.back().second);
                F.pop_back();
            }
            continue;
        }

        matrix_t M_prime(M);
        M_prime[i].clear();
        M_prime[i].insert(*mi);
        already_mapped.insert(*mi);
        size_t c_mi = *mi;

        // move locally to next child
        ++mi;

        size_t ni = i + 1;
        for (; ni < num_vertices(sub); ++ni)
            if (vertex(ni, sub) != graph_traits<Graph1>::null_vertex())
                break;

        // refine search tree
        if (refine_check(sub, g, M_prime, ni, already_mapped, edge_labelling,
                         vlist, vindex))
        {
            // store the current mapping so far
            F.push_back(std::make_pair(i, c_mi));

            if (ni < size_t(last_i))
            {
                // proceed with search at a higher depth
                Mstack.push_back(make_tuple(M_prime, ni, M_prime[ni].begin()));
                get<2>(Mstack.back()) = get<0>(Mstack.back())[ni].begin();
            }
            else
            {
                // maximum depth reached: visit all end leafs
                for (typeof(M_prime[ni].begin()) iter = M_prime[ni].begin();
                     iter != M_prime[ni].end(); ++iter)
                {
                    F.push_back(std::make_pair(ni, *iter));
                    FF.push_back(F);
                    F.pop_back();
                }

                // we are done which this tree node
                mi = M[i].end();
                F.pop_back();
                already_mapped.erase(c_mi);
            }
        }
        else
        {
            already_mapped.erase(c_mi);
        }
    }
}

}  // namespace detail


template <class Graph1, class Graph2, class VertexLabelling,
          class EdgeLabelling, class Mapping>
void subgraph_isomorphism(const Graph1& sub, const Graph2& g,
                          VertexLabelling vertex_labelling,
                          EdgeLabelling edge_labelling, vector<Mapping>& F,
                          vector<size_t>& vlist, size_t max_n)
{
    // initial mapping candidates
    detail::matrix_t M0(num_vertices(sub));
    vector<size_t> vindex(num_vertices(g));
    for (size_t j = 0; j < num_vertices(g); ++j)
        vindex[vlist[j]] = j;

    bool abort = false;
    int i, N = num_vertices(sub);
    #pragma omp parallel for default(shared) private(i) schedule(dynamic)
    for (i = 0; i < N; ++i)
    {
        if (vertex(i, sub) == graph_traits<Graph1>::null_vertex() || abort)
            continue;

        for (size_t j = 0; j < num_vertices(g); ++j)
        {
            if (vertex(vlist[j], g) == graph_traits<Graph1>::null_vertex())
                continue;
            if (vertex_labelling(vertex(i, sub), vertex(vlist[j], g)))
                M0[i].insert(j);
        }
        if (M0[i].empty())
            abort = true;
    }
    if (abort)
        return;
    detail::find_mappings(sub, g, M0, F, edge_labelling, vlist, vindex, max_n);
}

}  // namespace boost


#endif // SUBGRAPH_ISOMORPHISM_HPP
