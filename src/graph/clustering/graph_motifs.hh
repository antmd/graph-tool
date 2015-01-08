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

#ifndef GRAPH_MOTIFS_HH
#define GRAPH_MOTIFS_HH

#include <boost/graph/isomorphism.hpp>
#include <unordered_set>
#include <boost/functional/hash.hpp>
#include <algorithm>
#include <vector>

#include "random.hh"

namespace graph_tool
{
using namespace boost;

template <class Value>
void insert_sorted(std::vector<Value>& v, const Value& val)
{
    typeof(v.begin()) iter = lower_bound(v.begin(), v.end(), val);
    if (iter != v.end() && *iter == val)
        return; // no repetitions
    v.insert(iter, val);
}

template <class Value>
bool has_val(std::vector<Value>& v, const Value& val)
{
    typeof(v.begin()) iter = lower_bound(v.begin(), v.end(), val);
    if (iter == v.end())
        return false;
    return *iter == val;
}

// gets all the subgraphs starting from vertex v and store it in subgraphs.
template <class Graph, class Sampler>
void get_subgraphs(Graph& g, typename graph_traits<Graph>::vertex_descriptor v,
                   size_t n,
                   std::vector<std::vector<typename graph_traits<Graph>::vertex_descriptor> >& subgraphs,
                   Sampler sampler)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

    // extension and subgraph stack
    std::vector<std::vector<vertex_t> > ext_stack(1);
    std::vector<std::vector<vertex_t> > sub_stack(1);
    std::vector<std::vector<vertex_t> > sub_neighbours_stack(1);

    sub_stack[0].push_back(v);
    typename graph_traits<Graph>::out_edge_iterator e, e_end;
    for (tie(e, e_end) = out_edges(v, g); e != e_end; ++e)
    {
        typename graph_traits<Graph>::vertex_descriptor u = target(*e, g);
        if (u > v && !has_val(ext_stack[0], u))
        {
            insert_sorted(ext_stack[0], u);
            insert_sorted(sub_neighbours_stack[0],u);
        }
    }

    while (!sub_stack.empty())
    {
        std::vector<vertex_t>& ext = ext_stack.back();
        std::vector<vertex_t>& sub = sub_stack.back();
        std::vector<vertex_t>& sub_neighbours = sub_neighbours_stack.back();

        if (sub.size() == n)
        {
            // found a subgraph of the desired size; put it in the list and go
            // back a level
            subgraphs.push_back(sub);
            sub_stack.pop_back();
            ext_stack.pop_back();
            sub_neighbours_stack.pop_back();
            continue;
        }

        if (ext.empty())
        {
            // no where else to go
            ext_stack.pop_back();
            sub_stack.pop_back();
            sub_neighbours_stack.pop_back();
            continue;
        }
        else
        {
            // extend subgraph
            std::vector<vertex_t> new_ext, new_sub = sub,
                new_sub_neighbours = sub_neighbours;

            // remove w from ext
            vertex_t w = ext.back();
            ext.pop_back();

            // insert w in subgraph
            insert_sorted(new_sub, w);

            // update new_ext
            new_ext = ext;
            for (tie(e, e_end) = out_edges(w, g); e != e_end; ++e)
            {
                vertex_t u = target(*e,g);
                if (u > v)
                {
                    if (!has_val(sub_neighbours, u))
                        insert_sorted(new_ext, u);
                    insert_sorted(new_sub_neighbours, u);
                }
            }

            sampler(new_ext, ext_stack.size());

            ext_stack.push_back(new_ext);
            sub_stack.push_back(new_sub);
            sub_neighbours_stack.push_back(new_sub_neighbours);
        }
    }
}

// sampling selectors

struct sample_all
{
    template <class val_type>
    void operator()(std::vector<val_type>&, size_t) {}
};

struct sample_some
{
    sample_some(std::vector<double>& p, rng_t& rng): _p(&p), _rng(&rng) {}
    sample_some() {}

    template <class val_type>
    void operator()(std::vector<val_type>& extend, size_t d)
    {
        typedef std::uniform_real_distribution<double> rdist_t;
        auto random = std::bind(rdist_t(), std::ref(*_rng));

        double pd = (*_p)[d+1];
        size_t nc = extend.size();
        double u = nc*pd - floor(nc*pd);
        size_t n;
        double r;
        {
            #pragma omp critical
            r = random();
        }
        if (r < u)
            n = size_t(ceil(nc*pd));
        else
            n = size_t(floor(nc*pd));

        if (n == extend.size())
            return;
        if (n == 0)
        {
            extend.clear();
            return;
        }

        typedef std::uniform_int_distribution<size_t> idist_t;
        for (size_t i = 0; i < n; ++i)
        {
            auto random_v = std::bind(idist_t(0, extend.size()-i-1),
                                      std::ref(*_rng));
            size_t j;
            {
                #pragma omp critical
                j = i + random_v();
            }
            swap(extend[i], extend[j]);
        }
        extend.resize(n);
    }

    std::vector<double>* _p;
    rng_t* _rng;
};


// build the actual induced subgraph from the vertex list
template <class Graph, class GraphSG>
void make_subgraph
    (std::vector<typename graph_traits<Graph>::vertex_descriptor>& vlist,
     Graph& g, GraphSG& sub)
{
    for (size_t i = 0; i < vlist.size(); ++i)
        add_vertex(sub);
    for (size_t i = 0; i < vlist.size(); ++i)
    {
        typename graph_traits<Graph>::vertex_descriptor ov = vlist[i], ot;
        typename graph_traits<GraphSG>::vertex_descriptor nv = vertex(i,sub);
        typename graph_traits<Graph>::out_edge_iterator e, e_end;
        for (tie(e, e_end) = out_edges(ov, g); e != e_end; ++e)
        {
            ot = target(*e, g);
            typeof(vlist.begin()) viter =
                lower_bound(vlist.begin(), vlist.end(), ot);
            size_t ot_index = viter - vlist.begin();
            if (viter != vlist.end() && vlist[ot_index] == ot &&
                (is_directed::apply<Graph>::type::value || ot < ov))
                add_edge(nv, vertex(ot_index, sub), sub);
        }
    }
}

// compare two graphs for labeled exactness (not isomorphism)
template <class Graph>
bool graph_cmp(Graph& g1, Graph& g2)
{
    if (num_vertices(g1) != num_vertices(g2) || num_edges(g1) != num_edges(g2))
        return false;

    typename graph_traits<Graph>::vertex_iterator v1, v1_end;
    typename graph_traits<Graph>::vertex_iterator v2, v2_end;
    tie(v2, v2_end) = vertices(g2);
    for (tie(v1, v1_end) = vertices(g1); v1 != v1_end; ++v1)
    {
        if (out_degree(*v1, g1) != out_degree(*v2, g2))
            return false;
        if (in_degreeS()(*v1, g1) != in_degreeS()(*v2, g2))
            return false;

        std::vector<typename graph_traits<Graph>::vertex_descriptor> out1, out2;
        typename graph_traits<Graph>::out_edge_iterator e, e_end;
        for (tie(e, e_end) = out_edges(*v1, g1); e != e_end; ++e)
            out1.push_back(target(*e, g1));
        for (tie(e, e_end) = out_edges(*v2, g2); e != e_end; ++e)
            out2.push_back(target(*e, g2));
        sort(out1.begin(), out1.end());
        sort(out2.begin(), out2.end());
        if (out1 != out2)
            return false;
    }
    return true;
}

// short hand for both types of subgraphs
typedef adj_list<size_t> d_graph_t;
typedef adj_list<size_t> u_graph_t;

// we need this wrap to use the UndirectedAdaptor only on directed graphs
struct wrap_undirected
{
    template <class Graph>
    struct apply
    {
        typedef typename mpl::if_<typename is_directed::apply<Graph>::type,
                                  UndirectedAdaptor<Graph>,
                                  Graph&>::type type;
    };
};

// get the signature of the graph: sorted degree sequence
template <class Graph>
void get_sig(Graph& g, std::vector<size_t>& sig)
{
    sig.clear();
    size_t N = num_vertices(g);
    if (N > 0)
        sig.resize(is_directed::apply<Graph>::type::value ? 2 * N : N);
    for (size_t i = 0; i < N; ++i)
    {
        typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
        sig[i] = out_degree(v, g);
        if(is_directed::apply<Graph>::type::value)
            sig[i + N] = in_degreeS()(v, g);
    }
    sort(sig.begin(), sig.end());
}

// gets (or samples) all the subgraphs in graph g
struct get_all_motifs
{
    get_all_motifs(bool collect_vmaps, double p, bool comp_iso, bool fill_list,
                   rng_t& rng)
        : collect_vmaps(collect_vmaps), p(p),
          comp_iso(comp_iso), fill_list(fill_list), rng(rng) {}
    bool collect_vmaps;
    double p;
    bool comp_iso;
    bool fill_list;
    rng_t& rng;

    template <class Graph, class Sampler, class VMap>
    void operator()(Graph& g, size_t k, boost::any& list,
                    std::vector<size_t>& hist, std::vector<std::vector<VMap> >& vmaps,
                    Sampler sampler) const
    {
        typedef typename mpl::if_<typename is_directed::apply<Graph>::type,
                                  d_graph_t,
                                  u_graph_t>::type graph_sg_t;

        // the main subgraph lists
        std::vector<graph_sg_t>& subgraph_list =
            any_cast<std::vector<graph_sg_t>&>(list);

        // this hashes subgraphs according to their signature
        std::unordered_map<std::vector<size_t>,
                           std::vector<pair<size_t, graph_sg_t> >,
                           std::hash<std::vector<size_t>>> sub_list;
        std::vector<size_t> sig; // current signature

        for (size_t i = 0; i < subgraph_list.size(); ++i)
        {
            get_sig(subgraph_list[i], sig);
            sub_list[sig].push_back(make_pair(i,subgraph_list[i]));
        }

        // the subgraph count
        hist.resize(subgraph_list.size());

        typedef std::uniform_real_distribution<double> rdist_t;
        auto random = std::bind(rdist_t(), std::ref(rng));

        // the set of vertices V to be sampled (filled only if p < 1)
        std::vector<size_t> V;
        if (p < 1)
        {
            typename graph_traits<Graph>::vertex_iterator v, v_end;
            for (tie(v, v_end) = vertices(g); v != v_end; ++v)
                V.push_back(*v);

            size_t n;
            if (random() < p)
                n = size_t(ceil(V.size()*p));
            else
                n = size_t(floor(V.size()*p));

            typedef std::uniform_int_distribution<size_t> idist_t;
            for (size_t i = 0; i < n; ++i)
            {
                auto random_v = std::bind(idist_t(0, V.size()-i-1),
                                          std::ref(rng));
                size_t j = i + random_v();
                swap(V[i], V[j]);
            }
            V.resize(n);
        }

        int i, N = (p < 1) ? V.size() : num_vertices(g);
        #pragma omp parallel for default(shared) private(i, sig) \
            schedule(runtime) if (N > 100)
        for (i = 0; i < N; ++i)
        {
            std::vector<std::vector<typename graph_traits<Graph>::vertex_descriptor> >
                subgraphs;
            typename graph_traits<Graph>::vertex_descriptor v =
                (p < 1) ? V[i] : vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            typename wrap_undirected::apply<Graph>::type ug(g);
            get_subgraphs(ug, v, k, subgraphs, sampler);

            #pragma omp critical

            for (size_t j = 0; j < subgraphs.size(); ++j)
            {
                graph_sg_t sub;
                make_subgraph(subgraphs[j], g, sub);
                get_sig(sub, sig);

                typeof(sub_list.begin()) iter = sub_list.find(sig);
                if(iter == sub_list.end())
                {
                    if (!fill_list)
                        continue; // avoid inserting an element in sub_list
                    sub_list[sig].clear();
                }

                bool found = false;
                size_t pos;
                typeof(sub_list.begin()) sl = sub_list.find(sig);
                if (sl != sub_list.end())
                {
                    for (size_t l = 0; l < sl->second.size(); ++l)
                    {
                        graph_sg_t& motif = sl->second[l].second;
                        if (comp_iso)
                        {
                            if (isomorphism(motif, sub,
                                            vertex_index1_map(get(vertex_index, motif)).
                                            vertex_index2_map(get(vertex_index, sub))))
                                found = true;
                        }
                        else
                        {
                            if (graph_cmp(motif, sub))
                                found = true;
                        }
                        if (found)
                        {
                            pos = sl->second[l].first;
                            hist[pos]++;
                            break;
                        }
                    }
                }

                if (found == false && fill_list)
                {
                    subgraph_list.push_back(sub);
                    sub_list[sig].push_back(make_pair(subgraph_list.size() - 1,
                                                      sub));
                    hist.push_back(1);
                    pos = hist.size() - 1;
                    found = true;
                }

                if (found && collect_vmaps)
                {
                    if (pos >= vmaps.size())
                        vmaps.resize(pos + 1);
                    vmaps[pos].push_back(VMap(get(boost::vertex_index,sub)));
                    for (size_t vi = 0; vi < num_vertices(sub); ++vi)
                        vmaps[pos].back()[vertex(vi, sub)] = subgraphs[j][vi];
                }
            }
        }
    }
};

} //graph-tool namespace

#endif // GRAPH_MOTIFS_HH
