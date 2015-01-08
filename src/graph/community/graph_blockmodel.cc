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

#define BOOST_PYTHON_MAX_ARITY 40
#include <boost/python.hpp>
#include <cmath>
#include <iostream>

#include "numpy_bind.hh"

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "graph_filtering.hh"

#include "graph.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"
#include "graph_util.hh"

#include "random.hh"

#include "config.h"

#include "graph_blockmodel_overlap.hh"
#include "graph_blockmodel.hh"

using namespace boost;
using namespace graph_tool;



// ====================
// Entropy calculation
// ====================


// Repeated computation of x*log(x) and log(x) actually adds up to a lot of
// time. A significant speedup can be made by caching pre-computed values. This
// is doable since the values of mrse are bounded in [0, 2E], where E is the
// total number of edges in the network. Here we simply grow the cache as
// needed.

namespace graph_tool
{

vector<double> __safelog_cache;
vector<double> __xlogx_cache;
vector<double> __lgamma_cache;

void init_safelog(size_t x)
{
    size_t old_size = __safelog_cache.size();
    if (x >= old_size)
    {
        __safelog_cache.resize(x + 1);
        for (size_t i = old_size; i < __safelog_cache.size(); ++i)
            __safelog_cache[i] = safelog(double(i));
    }
}

void clear_safelog()
{
    vector<double>().swap(__safelog_cache);
}


void init_xlogx(size_t x)
{
    size_t old_size = __xlogx_cache.size();
    if (x >= old_size)
    {
        __xlogx_cache.resize(x + 1);
        for (size_t i = old_size; i < __xlogx_cache.size(); ++i)
            __xlogx_cache[i] = i * safelog(i);
    }
}

void clear_xlogx()
{
    vector<double>().swap(__xlogx_cache);
}

void init_lgamma(size_t x)
{
    size_t old_size = __lgamma_cache.size();
    if (x >= old_size)
    {
        __lgamma_cache.resize(x + 1);
        for (size_t i = old_size; i < __lgamma_cache.size(); ++i)
            __lgamma_cache[i] = lgamma(i);
    }
}

void clear_lgamma()
{
    vector<double>().swap(__lgamma_cache);
}


}

double do_get_ent(GraphInterface& gi, boost::any omrs, boost::any omrp,
                  boost::any omrm, boost::any owr, bool deg_corr)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        vmap_t;
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::edge_index_map_t>::type
        emap_t;
    emap_t mrs = any_cast<emap_t>(omrs);
    vmap_t mrp = any_cast<vmap_t>(omrp);
    vmap_t mrm = any_cast<vmap_t>(omrm);
    vmap_t wr = any_cast<vmap_t>(owr);

    double S = 0;
    run_action<>()
        (gi, std::bind(entropy(), mrs, mrp, mrm, wr, deg_corr,
                       placeholders::_1,
                       std::ref(S)))();
    return S;
}

double do_get_ent_dense(GraphInterface& gi, boost::any omrs, boost::any owr,
                        bool multigraph)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        vmap_t;
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::edge_index_map_t>::type
        emap_t;
    emap_t mrs = any_cast<emap_t>(omrs);
    vmap_t wr = any_cast<vmap_t>(owr);

    double S = 0;
    run_action<>()
        (gi, std::bind(entropy_dense(), mrs, wr, multigraph,
                       placeholders::_1, std::ref(S)))();
    return S;
}

double do_get_ent_parallel(GraphInterface& gi, boost::any oweight)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::edge_index_map_t>::type
        emap_t;
    emap_t weight = any_cast<emap_t>(oweight);

    double S = 0;
    run_action<>()
        (gi, std::bind(entropy_parallel_edges(),
                       placeholders::_1, weight, std::ref(S)))();
    return S;
}


boost::any do_create_emat(GraphInterface& gi)
{
    boost::any emat;
    run_action<>()(gi, std::bind<void>(create_emat(), placeholders::_1,
                                       std::ref(emat)))();
    return emat;
}

boost::any do_create_ehash(GraphInterface& gi)
{
    boost::any emat;
    run_action<>()(gi, std::bind<void>(create_ehash(), placeholders::_1,
                                       std::ref(emat)))();
    return emat;
}

//============
// Main loops
//============


template <class Eprop, class Vprop, class VEprop>
struct move_sweep_dispatch
{
    move_sweep_dispatch(Eprop eweight, Vprop vweight, boost::any egroups,
                        VEprop esrcpos, VEprop etgtpos, Vprop label,
                        vector<int>& vlist, bool deg_corr, bool dense,
                        bool multigraph, double beta, bool sequential,
                        bool parallel, bool random_move, double c, bool verbose,
                        size_t max_edge_index, size_t nmerges, size_t ntries,
                        Vprop merge_map, partition_stats_t& partition_stats,
                        rng_t& rng, double& S, size_t& nmoves,
                        GraphInterface& bgi)

        : eweight(eweight), vweight(vweight), oegroups(egroups), esrcpos(esrcpos),
          etgtpos(etgtpos), label(label), vlist(vlist),
          deg_corr(deg_corr), dense(dense), multigraph(multigraph), beta(beta),
          sequential(sequential), parallel(parallel), random_move(random_move),
          c(c), verbose(verbose), max_edge_index(max_edge_index),
          nmerges(nmerges), ntries(ntries), merge_map(merge_map),
          partition_stats(partition_stats), rng(rng), S(S),
          nmoves(nmoves), bgi(bgi)
    {}

    Eprop eweight;
    Vprop vweight;
    boost::any oegroups;
    VEprop esrcpos;
    VEprop etgtpos;
    Vprop label;
    size_t n;
    vector<int>& vlist;
    bool deg_corr;
    bool dense;
    bool multigraph;
    double beta;
    bool sequential;
    bool parallel;
    bool random_move;
    double c;
    bool verbose;
    size_t max_edge_index;
    size_t nmerges;
    size_t ntries;
    Vprop merge_map;
    partition_stats_t& partition_stats;
    rng_t& rng;
    double& S;
    size_t& nmoves;
    GraphInterface& bgi;

    template <class Graph>
    void operator()(Eprop mrs, Vprop mrp, Vprop mrm, Vprop wr, Vprop b,
                    Graph& g, boost::any& emat, boost::any sampler,
                    boost::any cavity_sampler, bool weighted) const
    {
        if (is_directed::apply<Graph>::type::value)
        {
            dispatch(mrs, mrp, mrm, wr, b, g, emat, sampler, cavity_sampler,
                     bgi.GetGraph(), weighted);
        }
        else
        {
            UndirectedAdaptor<GraphInterface::multigraph_t> ug(bgi.GetGraph());
            dispatch(mrs, mrp, mrm, wr, b, g, emat, sampler, cavity_sampler, ug,
                     weighted);
        }
    }

    template <class Graph, class BGraph>
    void dispatch(Eprop mrs, Vprop mrp, Vprop mrm, Vprop wr, Vprop b, Graph& g,
                  boost::any& aemat, boost::any asampler,
                  boost::any acavity_sampler, BGraph& bg, bool weighted) const
    {
        if (weighted)
        {
            typedef typename property_map_type::apply<DynamicSampler<std::tuple<typename graph_traits<Graph>::edge_descriptor, bool> >,
                                                      GraphInterface::vertex_index_map_t>::type vemap_t;
            vemap_t egroups = any_cast<vemap_t>(oegroups);
            dispatch(mrs, mrp, mrm, wr, b, g, aemat, asampler, acavity_sampler, bg, egroups);
        }
        else
        {
            typedef typename property_map_type::apply<vector<std::tuple<typename graph_traits<Graph>::edge_descriptor, bool> >,
                                                      GraphInterface::vertex_index_map_t>::type vemap_t;
            vemap_t egroups = any_cast<vemap_t>(oegroups);
            dispatch(mrs, mrp, mrm, wr, b, g, aemat, asampler, acavity_sampler, bg, egroups);
        }
    }

    template <class Graph, class BGraph, class Egroups>
    void dispatch(Eprop mrs, Vprop mrp, Vprop mrm, Vprop wr, Vprop b, Graph& g,
                  boost::any& aemat, boost::any asampler,
                  boost::any acavity_sampler, BGraph& bg, Egroups egroups) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

        size_t B = num_vertices(bg);
        size_t max_BE = is_directed::apply<Graph>::type::value ?
            B * B : (B * (B + 1)) / 2;

        size_t eidx = random_move ? 1 : max_edge_index;

        typedef typename property_map<Graph, vertex_index_t>::type vindex_map_t;
        typedef typename property_map_type::apply<Sampler<vertex_t, boost::mpl::false_>,
                                                  vindex_map_t>::type::unchecked_t
            sampler_map_t;
        sampler_map_t sampler = any_cast<sampler_map_t>(asampler);
        sampler_map_t cavity_sampler = any_cast<sampler_map_t>(acavity_sampler);

        try
        {
            typedef typename get_emat_t::apply<BGraph>::type emat_t;
            emat_t& emat = any_cast<emat_t&>(aemat);

            // make sure the properties are _unchecked_, since otherwise it
            // affects performance

            move_sweep(mrs.get_unchecked(max_BE),
                       mrp.get_unchecked(num_vertices(bg)),
                       mrm.get_unchecked(num_vertices(bg)),
                       wr.get_unchecked(num_vertices(bg)),
                       b.get_unchecked(num_vertices(g)),
                       label.get_unchecked(num_vertices(bg)), vlist, deg_corr,
                       dense, multigraph, beta,
                       eweight.get_unchecked(max_edge_index),
                       vweight.get_unchecked(num_vertices(g)),
                       egroups.get_unchecked(num_vertices(bg)),
                       esrcpos.get_unchecked(eidx),
                       etgtpos.get_unchecked(eidx), g, bg, emat, sampler,
                       cavity_sampler, sequential, parallel, random_move, c,
                       nmerges, ntries,
                       merge_map.get_unchecked(num_vertices(g)),
                       partition_stats, verbose, rng, S, nmoves,
                       overlap_stats_t());
        }
        catch (bad_any_cast&)
        {
            typedef typename get_ehash_t::apply<BGraph>::type emat_t;
            emat_t& emat = any_cast<emat_t&>(aemat);
            move_sweep(mrs.get_unchecked(num_edges(g)),
                       mrp.get_unchecked(num_vertices(bg)),
                       mrm.get_unchecked(num_vertices(bg)),
                       wr.get_unchecked(num_vertices(bg)),
                       b.get_unchecked(num_vertices(g)),
                       label.get_unchecked(num_vertices(bg)), vlist, deg_corr,
                       dense, multigraph, beta,
                       eweight.get_unchecked(max_edge_index),
                       vweight.get_unchecked(num_vertices(g)),
                       egroups.get_unchecked(num_vertices(bg)),
                       esrcpos.get_unchecked(eidx),
                       etgtpos.get_unchecked(eidx), g, bg, emat, sampler,
                       cavity_sampler, sequential, parallel, random_move, c,
                       nmerges, ntries,
                       merge_map.get_unchecked(num_vertices(g)),
                       partition_stats, verbose, rng, S, nmoves,
                       overlap_stats_t());
        }
    }
};


boost::python::object do_move_sweep(GraphInterface& gi, GraphInterface& bgi,
                                    boost::any& emat, boost::any sampler,
                                    boost::any cavity_sampler, boost::any omrs,
                                    boost::any omrp, boost::any omrm,
                                    boost::any owr, boost::any ob,
                                    boost::any olabel, vector<int>& vlist,
                                    bool deg_corr, bool dense, bool multigraph,
                                    boost::any oeweight, boost::any ovweight,
                                    boost::any oegroups, boost::any oesrcpos,
                                    boost::any oetgtpos, double beta,
                                    bool sequential, bool parallel,
                                    bool random_move, double c, bool weighted,
                                    size_t nmerges, size_t ntries,
                                    boost::any omerge_map,
                                    partition_stats_t& partition_stats,
                                    bool verbose, rng_t& rng)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        vmap_t;
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::edge_index_map_t>::type
        emap_t;
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::edge_index_map_t>::type
        vemap_t;
    emap_t mrs = any_cast<emap_t>(omrs);
    vmap_t mrp = any_cast<vmap_t>(omrp);
    vmap_t mrm = any_cast<vmap_t>(omrm);
    vmap_t wr = any_cast<vmap_t>(owr);
    vmap_t b = any_cast<vmap_t>(ob);
    vmap_t label = any_cast<vmap_t>(olabel);
    emap_t eweight = any_cast<emap_t>(oeweight);
    vmap_t vweight = any_cast<vmap_t>(ovweight);

    vemap_t esrcpos = any_cast<vemap_t>(oesrcpos);
    vemap_t etgtpos = any_cast<vemap_t>(oetgtpos);

    double S = 0;
    size_t nmoves = 0;

    vmap_t merge_map = any_cast<vmap_t>(omerge_map);

    run_action<graph_tool::detail::all_graph_views, boost::mpl::true_>()
        (gi, std::bind(move_sweep_dispatch<emap_t, vmap_t, vemap_t>
                       (eweight, vweight, oegroups, esrcpos, etgtpos,
                        label, vlist, deg_corr, dense, multigraph, beta,
                        sequential, parallel, random_move, c, verbose,
                        gi.GetMaxEdgeIndex(), nmerges, ntries, merge_map,
                        partition_stats, rng, S, nmoves, bgi),
                       mrs, mrp, mrm, wr, b, placeholders::_1,
                       std::ref(emat), sampler, cavity_sampler, weighted))();
    return boost::python::make_tuple(S, nmoves);
}

struct build_egroups
{
    template <class Eprop, class Vprop, class VEprop, class Graph, class VertexIndex>
    void operator()(Vprop b, boost::any& oegroups, VEprop esrcpos,
                    VEprop etgtpos, Eprop eweight, Graph& g,
                    VertexIndex vertex_index, size_t B, bool weighted, bool empty) const
    {
        egroups_manage::build(b, oegroups, esrcpos, etgtpos, eweight, g,
                              vertex_index, B, weighted, empty);
    }
};

boost::any do_build_egroups(GraphInterface& gi, GraphInterface& bgi,
                            boost::any ob, boost::any oeweights,
                            boost::any oesrcpos, boost::any oetgtpos,
                            bool weighted, bool empty)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        vmap_t;
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::edge_index_map_t>::type
        emap_t;
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::edge_index_map_t>::type
        vemap_t;
    vmap_t b = any_cast<vmap_t>(ob);

    vemap_t esrcpos = any_cast<vemap_t>(oesrcpos);
    vemap_t etgtpos = any_cast<vemap_t>(oetgtpos);
    emap_t eweights = any_cast<emap_t>(oeweights);

    boost::any oegroups;
    run_action<graph_tool::detail::all_graph_views, boost::mpl::true_>()
        (gi, std::bind<void>(build_egroups(), b, std::ref(oegroups),
                             esrcpos.get_unchecked(gi.GetMaxEdgeIndex()),
                             etgtpos.get_unchecked(gi.GetMaxEdgeIndex()),
                             eweights.get_unchecked(gi.GetMaxEdgeIndex()),
                             placeholders::_1, bgi.GetVertexIndex(),
                             bgi.GetNumberOfVertices(), weighted, empty))();
    return oegroups;
}

boost::any do_init_neighbour_sampler(GraphInterface& gi, boost::any oeweights,
                                     bool self_loops, bool empty)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::edge_index_map_t>::type
        emap_t;
    emap_t eweights = any_cast<emap_t>(oeweights);

    boost::any osampler;
    run_action<graph_tool::detail::all_graph_views, boost::mpl::true_>()
        (gi, std::bind(init_neighbour_sampler(), placeholders::_1, eweights,
                       self_loops, empty, std::ref(osampler)))();
    return osampler;
}

struct collect_edge_marginals_dispatch
{
    template <class Graph, class Vprop, class MEprop>
    void operator()(Graph& g, size_t B, Vprop cb, MEprop p,
                    std::tuple<boost::any, GraphInterface&> abg) const
    {
        Graph& bg = *any_cast<Graph*>(get<0>(abg));
        collect_edge_marginals(B, cb.get_unchecked(num_vertices(bg)), p, g, bg);
    }
};

void do_collect_edge_marginals(GraphInterface& gi, GraphInterface& gbi,
                               size_t B, boost::any ob, boost::any op)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        vmap_t;
    typedef property_map_type::apply<vector<int32_t>,
                                     GraphInterface::edge_index_map_t>::type
        emap_t;
    vmap_t b = any_cast<vmap_t>(ob);
    emap_t p = any_cast<emap_t>(op);

    run_action<graph_tool::detail::all_graph_views, boost::mpl::true_>()
        (gi, std::bind<void>(collect_edge_marginals_dispatch(),
                             placeholders::_1, B, b, p,
                             std::tuple<boost::any, GraphInterface&>(gbi.GetGraphView(), gbi)))();
}

boost::python::tuple do_bethe_entropy(GraphInterface& gi, size_t B, boost::any op,
                                      boost::any opv)
{
    typedef property_map_type::apply<vector<double>,
                                     GraphInterface::vertex_index_map_t>::type
        vmap_t;
    typedef property_map_type::apply<vector<int32_t>,
                                     GraphInterface::edge_index_map_t>::type
        emap_t;
    emap_t p = any_cast<emap_t>(op);
    vmap_t pv = any_cast<vmap_t>(opv);

    double H=0, sH=0, Hmf=0, sHmf=0;
    run_action<graph_tool::detail::all_graph_views, boost::mpl::true_>()
        (gi, std::bind<void>(bethe_entropy(),
                             placeholders::_1, B, p, pv, std::ref(H), std::ref(sH),
                             std::ref(Hmf), std::ref(sHmf)))();
    return boost::python::make_tuple(H, sH, Hmf, sHmf);
}


struct collect_vertex_marginals_dispatch
{
    template <class Graph, class Vprop, class MEprop>
    void operator()(Graph& g, Vprop cb, MEprop p) const
    {
        collect_vertex_marginals(cb.get_unchecked(num_vertices(g)),
                                 p.get_unchecked(num_vertices(g)), g);
    }
};


void do_collect_vertex_marginals(GraphInterface& gi, boost::any ob,
                                 boost::any op)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        vmap_t;
    vmap_t b = any_cast<vmap_t>(ob);

    run_action<graph_tool::detail::all_graph_views, boost::mpl::true_>()
        (gi, std::bind(collect_vertex_marginals_dispatch(),
                       placeholders::_1, b, placeholders::_2),
         vertex_scalar_vector_properties())(op);
}


struct get_deg_entropy_term
{
    template <class Graph>
    void operator()(Graph& g, double& S) const
    {
        for(auto v : vertices_range(g))
        {
            S -= lgamma_fast(out_degree(v, g) + 1);
            S -= lgamma_fast(in_degreeS()(v, g) + 1);
        }
    }
};

struct get_deg_entropy_term_overlap
{

    template <class Graph, class Vprop>
    void operator()(Graph& g, Vprop b, overlap_stats_t& overlap_stats, size_t N,
                    double& S) const
    {
#ifdef HAVE_SPARSEHASH
        typedef dense_hash_map<int, int> map_t;
#else
        typedef unordered_map<int, int> map_t;
#endif

        map_t in_hist, out_hist;

#ifdef HAVE_SPARSEHASH
        in_hist.set_empty_key(numeric_limits<int>::max());
        out_hist.set_empty_key(numeric_limits<int>::max());
#endif

        for (size_t v = 0; v < N; ++v)
        {
            in_hist.clear();
            out_hist.clear();

            auto& half_edges = overlap_stats.get_half_edges(v);
            for (size_t u : half_edges)
            {
                in_hist[b[u]] += in_degreeS()(u, g);
                out_hist[b[u]] += out_degree(u, g);
            }

            for (auto& k_c : in_hist)
                S -= lgamma_fast(k_c.second + 1);
            for (auto& k_c : out_hist)
                S -= lgamma_fast(k_c.second + 1);
        }
    }
};

double do_get_deg_entropy_term(GraphInterface& gi, boost::any ob,
                               overlap_stats_t& overlap_stats, size_t N)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        vmap_t;

    double S = 0;

    if (!ob.empty())
    {
        vmap_t b = any_cast<vmap_t>(ob);
        run_action<>()
            (gi, std::bind(get_deg_entropy_term_overlap(),
                           placeholders::_1, b, std::ref(overlap_stats), N,
                           std::ref(S)))();
    }
    else
    {
        run_action<>()
            (gi, std::bind(get_deg_entropy_term(),
                           placeholders::_1, std::ref(S)))();
    }

    return S;
}


vector<int32_t> get_vector(size_t n)
{
    return vector<int32_t>(n);
}

void vector_map(boost::python::object ovals, boost::python::object omap)
{

    multi_array_ref<int32_t,1> vals = get_array<int32_t,1>(ovals);
    multi_array_ref<int32_t,1> map = get_array<int32_t,1>(omap);

    size_t pos = 0;
    for (size_t i = 0; i < vals.size(); ++i)
    {
        int32_t v = vals[i];
        if (map[v] == -1)
            map[v] = pos++;
        vals[i] = map[v];
    }
}

void vector_continuous_map(boost::python::object ovals)
{

    multi_array_ref<int32_t,1> vals = get_array<int32_t,1>(ovals);
    unordered_map<int32_t, size_t> map;

    for (size_t i = 0; i < vals.size(); ++i)
    {
        int32_t v = vals[i];
        auto iter = map.find(v);
        if (iter == map.end())
            iter = map.insert(make_pair(v, map.size())).first;
        vals[i] = iter->second;
    }
}

void vector_rmap(boost::python::object ovals, boost::python::object omap)
{

    multi_array_ref<int32_t,1> vals = get_array<int32_t,1>(ovals);
    multi_array_ref<int32_t,1> map = get_array<int32_t,1>(omap);

    for (size_t i = 0; i < vals.size(); ++i)
    {
        map[vals[i]] = i;
    }
}

python::tuple python_get_mu_l(double N, double E, double epsilon=1e-8)
{
    double mu, l;
    get_mu_l(N, E, mu, l, epsilon);
    return python::make_tuple(mu, l);
}

struct get_partition_stats
{
    template <class Graph, class Vprop, class Eprop>
    void operator()(Graph& g, Vprop b, Eprop eweight, size_t N, size_t B,
                    partition_stats_t& partition_stats) const
    {
        partition_stats = partition_stats_t(g, b, eweight, N, B);
    }
};

partition_stats_t
do_get_partition_stats(GraphInterface& gi, boost::any ob, boost::any aeweight,
                       size_t N, size_t B)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        vmap_t;
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::edge_index_map_t>::type
        emap_t;

    partition_stats_t partition_stats;

    vmap_t b = any_cast<vmap_t>(ob);
    emap_t eweight = any_cast<emap_t>(aeweight);

    run_action<>()(gi, std::bind(get_partition_stats(),
                                 placeholders::_1, b, eweight, N, B,
                                 std::ref(partition_stats)))();
    return partition_stats;
}

void export_blockmodel()
{
    using namespace boost::python;

    class_<partition_stats_t>("partition_stats")
        .def("is_enabled", &partition_stats_t::is_enabled)
        .def("get_partition_dl", &partition_stats_t::get_partition_dl)
        .def("get_deg_dl", &partition_stats_t::get_deg_dl);

    def("init_partition_stats", do_get_partition_stats);

    def("init_safelog", init_safelog);
    def("clear_safelog", clear_safelog);
    def("init_xlogx", init_xlogx);
    def("clear_xlogx", clear_xlogx);
    def("init_lgamma", init_lgamma);
    def("clear_lgamma", clear_lgamma);
    def("get_xi", get_xi<double,double>);
    def("get_xi_fast", get_xi_fast<double,double>);
    def("get_mu_l", python_get_mu_l);
    def("polylog", polylog<double>);

    def("poisson_entropy", poisson_entropy<double>);
    def("lpoisson", lpoisson<double>);
    def("poisson", poisson<double>);

    def("get_vector", get_vector);
    def("vector_map", vector_map);
    def("vector_rmap", vector_rmap);
    def("vector_continuous_map", vector_continuous_map);

    def("create_emat", do_create_emat);
    def("create_ehash", do_create_ehash);
    def("build_egroups", do_build_egroups);
    def("init_neighbour_sampler", do_init_neighbour_sampler);

    def("move_sweep", do_move_sweep);

    def("entropy", do_get_ent);
    def("entropy_dense", do_get_ent_dense);
    def("entropy_parallel", do_get_ent_parallel);
    def("deg_entropy_term", do_get_deg_entropy_term);

    def("edge_marginals", do_collect_edge_marginals);
    def("bethe_entropy", do_bethe_entropy);
    def("vertex_marginals", do_collect_vertex_marginals);
}
