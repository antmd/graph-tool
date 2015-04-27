// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2016 Tiago de Paula Peixoto <tiago@skewed.de>
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

#define BOOST_PYTHON_MAX_ARITY 46
#include <boost/python.hpp>
#include <cmath>
#include <iostream>
#include <functional>

#include "numpy_bind.hh"

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "graph_filtering.hh"

#include "graph.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"
#include "graph_util.hh"
#include "graph_python_interface.hh"

#include "random.hh"

#include "config.h"

#include "graph_blockmodel_covariates.hh"
#include "graph_blockmodel.hh"
#include "graph_blockmodel_overlap.hh"

using namespace graph_tool;
using namespace std;


template <class Eprop, class Vprop, class VVprop, class VEprop, class BMap>
struct cov_move_sweep_dispatch
{
    cov_move_sweep_dispatch(Eprop ce, VVprop cv, VVprop vmap,
                            vector<std::reference_wrapper<Eprop>>& eweight,
                            vector<std::reference_wrapper<Vprop>>& vweight,
                            vector<std::reference_wrapper<boost::any>>& egroups,
                            vector<std::reference_wrapper<VEprop>>& esrcpos,
                            vector<std::reference_wrapper<VEprop>>& etgtpos,
                            Vprop& label, vector<int>& vlist, bool deg_corr,
                            bool dense, bool multigraph, double beta,
                            bool sequential, bool parallel, bool random_move,
                            bool node_coherent, bool confine_layers, double c,
                            bool verbose, size_t meindex,
                            vector<size_t> max_edge_index, size_t nmerges,
                            size_t niter, Vprop merge_map,
                            vector<std::reference_wrapper<partition_stats_t>>& partition_stats,
                            vector<std::reference_wrapper<overlap_partition_stats_t>>& overlap_partition_stats,
                            vector<std::reference_wrapper<overlap_stats_t>>& overlap_stats,
                            vector<bool>& master, vector<bool>& slave,
                            rng_t& rng, double& S, size_t& nmoves,
                            vector<std::reference_wrapper<GraphInterface>>& bgi,
                            BMap& block_map,
                            vector<std::reference_wrapper<Vprop>>& block_rmap,
                            vector<std::reference_wrapper<vector<size_t>>>& free_blocks,
                            size_t B)

        : ce(ce), cv(cv), vmap(vmap), eweight(eweight), vweight(vweight),
          oegroups(egroups), esrcpos(esrcpos), etgtpos(etgtpos), label(label), vlist(vlist),
          deg_corr(deg_corr), dense(dense), multigraph(multigraph), beta(beta),
          sequential(sequential), parallel(parallel), random_move(random_move),
          node_coherent(node_coherent), confine_layers(confine_layers),
          c(c), verbose(verbose), meindex(meindex),
          max_edge_index(max_edge_index),
          nmerges(nmerges), niter(niter), merge_map(merge_map),
          partition_stats(partition_stats),
          overlap_partition_stats(overlap_partition_stats),
          overlap_stats(overlap_stats),
          master(master), slave(slave),
          rng(rng), S(S), nmoves(nmoves), bgi(bgi), block_map(block_map),
          block_rmap(block_rmap), free_blocks(free_blocks), B(B)
    {}

    Eprop ce;
    VVprop cv;
    VVprop vmap;
    vector<std::reference_wrapper<Eprop>>& eweight;
    vector<std::reference_wrapper<Vprop>>& vweight;
    vector<std::reference_wrapper<boost::any>> oegroups;
    vector<std::reference_wrapper<VEprop>>& esrcpos;
    vector<std::reference_wrapper<VEprop>>& etgtpos;
    Vprop& label;
    size_t n;
    vector<int>& vlist;
    bool deg_corr;
    bool dense;
    bool multigraph;
    double beta;
    bool sequential;
    bool parallel;
    bool random_move;
    bool node_coherent;
    bool confine_layers;
    double c;
    bool verbose;
    size_t meindex;
    vector<size_t> max_edge_index;
    size_t nmerges;
    size_t niter;
    Vprop merge_map;
    vector<std::reference_wrapper<partition_stats_t>>& partition_stats;
    vector<std::reference_wrapper<overlap_partition_stats_t>>& overlap_partition_stats;
    vector<std::reference_wrapper<overlap_stats_t>>& overlap_stats;
    vector<bool> master;
    vector<bool> slave;
    rng_t& rng;
    double& S;
    size_t& nmoves;
    vector<std::reference_wrapper<GraphInterface>>& bgi;
    BMap& block_map;
    vector<std::reference_wrapper<Vprop>>& block_rmap;
    vector<std::reference_wrapper<vector<size_t>>>& free_blocks;
    size_t B;

    template <class Graph>
    void operator()(vector<std::reference_wrapper<Eprop>>& mrs,
                    vector<std::reference_wrapper<Vprop>>& mrp,
                    vector<std::reference_wrapper<Vprop>>& mrm,
                    vector<std::reference_wrapper<Vprop>>& wr,
                    Vprop& b, vector<std::reference_wrapper<Vprop>>& bs,
                    vector<std::reference_wrapper<GraphInterface>>& bgis, Graph& g,
                    vector<std::reference_wrapper<GraphInterface>>& ags,
                    vector<std::reference_wrapper<boost::any>>& emat,
                    vector<std::reference_wrapper<boost::any>>& sampler,
                    vector<std::reference_wrapper<boost::any>>& cavity_sampler,
                    bool weighted) const
    {
        vector<std::reference_wrapper<Graph>> gs;
        for (GraphInterface& ag : ags)
            gs.push_back(*any_cast<Graph*>(ag.GetGraphView()));

        if (is_directed::apply<Graph>::type::value)
        {
            vector<std::reference_wrapper<GraphInterface::multigraph_t>> bgs;
            for (GraphInterface& bgi : bgis)
                bgs.push_back(bgi.GetGraph());
            dispatch(mrs, mrp, mrm, wr, b, bs, g, gs, emat, sampler,
                     cavity_sampler, bgs, weighted);
        }
        else
        {
            vector<UndirectedAdaptor<GraphInterface::multigraph_t>> ubgs;
            for (GraphInterface& bgi : bgis)
                ubgs.push_back(UndirectedAdaptor<GraphInterface::multigraph_t>(bgi.GetGraph()));
            vector<std::reference_wrapper<UndirectedAdaptor<GraphInterface::multigraph_t>>> rubgs;
            for (auto& bg : ubgs)
                rubgs.push_back(bg);
            dispatch(mrs, mrp, mrm, wr, b, bs, g, gs, emat, sampler,
                     cavity_sampler, rubgs, weighted);
        }
    }

    template <class Graph, class BGraph>
    void dispatch(vector<std::reference_wrapper<Eprop>>& mrs,
                  vector<std::reference_wrapper<Vprop>>& mrp,
                  vector<std::reference_wrapper<Vprop>>& mrm,
                  vector<std::reference_wrapper<Vprop>>& wr,
                  Vprop& b,
                  vector<std::reference_wrapper<Vprop>>& bs, Graph& g,
                  vector<std::reference_wrapper<Graph>>& gs,
                  vector<std::reference_wrapper<boost::any>>& aemat,
                  vector<std::reference_wrapper<boost::any>>& asampler,
                  vector<std::reference_wrapper<boost::any>>& acavity_sampler,
                  vector<std::reference_wrapper<BGraph>>& bg,
                  bool weighted) const
    {
        if (weighted)
        {
            typedef typename property_map_type::apply<DynamicSampler<std::tuple<typename graph_traits<Graph>::edge_descriptor, bool> >,
                                                      GraphInterface::vertex_index_map_t>::type vemap_t;
            vector<std::reference_wrapper<vemap_t>> egroups;
            for (auto& eg : oegroups)
                egroups.push_back(any_cast<vemap_t&>(eg));

            try
            {
                typedef typename get_emat_t::apply<BGraph>::type emat_t;
                vector<std::reference_wrapper<emat_t>> emat;
                for (auto& m : aemat)
                    emat.push_back(any_cast<emat_t&>(m));
                size_t B = num_vertices(bg[0].get());
                size_t max_BE = is_directed::apply<Graph>::type::value ?
                    B * B : (B * (B + 1)) / 2;
                vector<typename Eprop::unchecked_t> umrs;
                for (auto& m : mrs)
                    umrs.push_back(m.get().get_unchecked(max_BE));
                vector<std::reference_wrapper<typename Eprop::unchecked_t>> rumrs;
                for (auto& m : umrs)
                    rumrs.push_back(m);
                dispatch(rumrs, mrp, mrm, wr, b, bs, g, gs, asampler,
                         acavity_sampler, bg, egroups, emat);
            }
            catch (bad_any_cast)
            {
                typedef typename get_ehash_t::apply<BGraph>::type emat_t;
                vector<std::reference_wrapper<emat_t>> emat;
                for (auto& m : aemat)
                    emat.push_back(any_cast<emat_t&>(m));
                dispatch(mrs, mrp, mrm, wr, b, bs, g, gs, asampler,
                         acavity_sampler, bg, egroups, emat);
            }
        }
        else
        {
            typedef typename property_map_type::apply<vector<std::tuple<typename graph_traits<Graph>::edge_descriptor, bool> >,
                                                      GraphInterface::vertex_index_map_t>::type vemap_t;
            vector<std::reference_wrapper<vemap_t>> egroups;
            for (auto& eg : oegroups)
                egroups.push_back(any_cast<vemap_t&>(eg));

            try
            {
                typedef typename get_emat_t::apply<BGraph>::type emat_t;
                vector<std::reference_wrapper<emat_t>> emat;
                for (auto& m : aemat)
                    emat.push_back(any_cast<emat_t&>(m));
                dispatch(mrs, mrp, mrm, wr, b, bs, g, gs, asampler,
                         acavity_sampler, bg, egroups, emat);
            }
            catch (bad_any_cast)
            {
                typedef typename get_ehash_t::apply<BGraph>::type emat_t;
                vector<std::reference_wrapper<emat_t>> emat;
                for (auto& m : aemat)
                    emat.push_back(any_cast<emat_t&>(m));
                dispatch(mrs, mrp, mrm, wr, b, bs, g, gs, asampler,
                         acavity_sampler, bg, egroups, emat);
            }
        }
    }

    template <class Graph, class BGraph, class Egroups, class Emat, class EMprop>
    void dispatch(vector<std::reference_wrapper<EMprop>>& mrs,
                  vector<std::reference_wrapper<Vprop>>& mrp,
                  vector<std::reference_wrapper<Vprop>>& mrm,
                  vector<std::reference_wrapper<Vprop>>& wr,
                  Vprop& b,
                  vector<std::reference_wrapper<Vprop>>& bs, Graph& g,
                  vector<std::reference_wrapper<Graph>>& gs,
                  vector<std::reference_wrapper<boost::any>>& asampler,
                  vector<std::reference_wrapper<boost::any>>& acavity_sampler,
                  vector<std::reference_wrapper<BGraph>>& bgs,
                  vector<std::reference_wrapper<Egroups>>& egroups,
                  vector<std::reference_wrapper<Emat>>& emat) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

        typedef typename property_map<Graph, vertex_index_t>::type vindex_map_t;
        typedef typename property_map_type::apply<Sampler<vertex_t, boost::mpl::false_>,
                                                  vindex_map_t>::type::unchecked_t
            sampler_map_t;

        vector<std::reference_wrapper<sampler_map_t>> sampler;
        for (auto& s : asampler)
            sampler.push_back(any_cast<sampler_map_t&>(s.get()));
        vector<std::reference_wrapper<sampler_map_t>> cavity_sampler;
        for (auto& s : acavity_sampler)
            cavity_sampler.push_back(any_cast<sampler_map_t&>(s.get()));

        bool overlap = overlap_stats[0].get().is_enabled();
        if (!overlap)
        {
            typedef BlockState<Graph, BGraph, EMprop,
                               typename Eprop::unchecked_t,
                               typename Vprop::unchecked_t, Emat,
                               typename Egroups::unchecked_t,
                               typename VEprop::unchecked_t, sampler_map_t,
                               partition_stats_t,
                               overlap_stats_t,
                               typename BMap::value_type,
                               typename Vprop::unchecked_t> state_t;

            vector<state_t> states;
            vector<EntrySet<Graph>> m_entries;
            overlap_stats_t ostats;

            for (size_t i = 0; i < mrs.size(); ++i)
            {
                size_t eidx = random_move ? 1 : max_edge_index[i];

                state_t state = make_block_state(gs[i].get(),
                                                 eweight[i].get().get_unchecked(max_edge_index[i]),
                                                 vweight[i].get().get_unchecked(num_vertices(gs[i].get())),
                                                 bs[i].get().get_unchecked(num_vertices(gs[i].get())),
                                                 bgs[i].get(),
                                                 emat[i].get(),
                                                 mrs[i].get(),
                                                 mrp[i].get().get_unchecked(num_vertices(bgs[i].get())),
                                                 mrm[i].get().get_unchecked(num_vertices(bgs[i].get())),
                                                 wr[i].get().get_unchecked(num_vertices(bgs[i].get())),
                                                 egroups[i].get().get_unchecked(num_vertices(bgs[i].get())),
                                                 esrcpos[i].get().get_unchecked(eidx),
                                                 etgtpos[i].get().get_unchecked(eidx),
                                                 sampler[i].get(), cavity_sampler[i].get(),
                                                 partition_stats[i].get(), ostats,
                                                 block_map[i],
                                                 block_rmap[i].get().get_unchecked(num_vertices(bgs[i].get())),
                                                 free_blocks[i].get(),
                                                 master[i],
                                                 slave[i],
                                                 false);
                states.push_back(state);
                m_entries.emplace_back(num_vertices(bgs[i].get()));
            }

            move_sweep(states, m_entries,
                       wr[0].get().get_unchecked(B),
                       b.get_unchecked(num_vertices(g)),
                       cv.get_unchecked(num_vertices(g)),
                       vmap.get_unchecked(num_vertices(g)),
                       label.get_unchecked(B),
                       vlist, deg_corr,
                       dense, multigraph, beta,
                       eweight[0].get().get_unchecked(max_edge_index[0]),
                       vweight[0].get().get_unchecked(num_vertices(g)),
                       g, sequential, parallel, random_move, c,
                       nmerges,
                       merge_map.get_unchecked(num_vertices(g)),
                       niter, B,
                       verbose, rng, S, nmoves);
        }
        else
        {
            typedef BlockState<Graph, BGraph, EMprop,
                               typename Eprop::unchecked_t,
                               typename Vprop::unchecked_t, Emat,
                               typename Egroups::unchecked_t,
                               typename VEprop::unchecked_t, sampler_map_t,
                               overlap_partition_stats_t,
                               overlap_stats_t,
                               typename BMap::value_type,
                               typename Vprop::unchecked_t> state_t;

            vector<state_t> states;
            for (size_t i = 0; i < mrs.size(); ++i)
            {
                size_t eidx = random_move ? 1 : max_edge_index[i];

                state_t state = make_block_state(gs[i].get(),
                                                 eweight[i].get().get_unchecked(max_edge_index[i]),
                                                 vweight[i].get().get_unchecked(num_vertices(gs[i].get())),
                                                 bs[i].get().get_unchecked(num_vertices(gs[i].get())),
                                                 bgs[i].get(),
                                                 emat[i].get(),
                                                 mrs[i].get(),
                                                 mrp[i].get().get_unchecked(num_vertices(bgs[i].get())),
                                                 mrm[i].get().get_unchecked(num_vertices(bgs[i].get())),
                                                 wr[i].get().get_unchecked(num_vertices(bgs[i].get())),
                                                 egroups[i].get().get_unchecked(num_vertices(bgs[i].get())),
                                                 esrcpos[i].get().get_unchecked(eidx),
                                                 etgtpos[i].get().get_unchecked(eidx),
                                                 sampler[i].get(), cavity_sampler[i].get(),
                                                 overlap_partition_stats[i].get(),
                                                 overlap_stats[i].get(),
                                                 block_map[i],
                                                 block_rmap[i].get().get_unchecked(num_vertices(bgs[i].get())),
                                                 free_blocks[i].get(),
                                                 master[i],
                                                 slave[i],
                                                 false);
                states.push_back(state);
            }
            vector<SingleEntrySet<Graph>> m_entries(states.size());

            if (nmerges == 0)
            {
                if (!node_coherent)
                {
                    move_sweep_overlap(states, m_entries, overlap_stats[0].get(),
                                       wr[0].get().get_unchecked(B),
                                       b.get_unchecked(num_vertices(g)), cv,
                                       vmap, label.get_unchecked(B),
                                       vlist, deg_corr, dense, multigraph, beta,
                                       vweight[0].get().get_unchecked(num_vertices(g)), g,
                                       sequential, parallel, random_move, c, niter,
                                       B, verbose, rng, S, nmoves);
                }
                else
                {
                    vector<EntrySet<Graph>> m_entries;
                    for (auto& state : states)
                        m_entries.emplace_back(num_vertices(state.bg));
                    coherent_move_sweep_overlap(states, m_entries, overlap_stats[0].get(),
                                                wr[0].get().get_unchecked(B),
                                                b.get_unchecked(num_vertices(g)), cv,
                                                vmap, label.get_unchecked(B),
                                                vlist, deg_corr, dense, multigraph, beta,
                                                vweight[0].get().get_unchecked(num_vertices(g)), g,
                                                sequential, random_move, c,
                                                confine_layers, niter,
                                                B, rng, S, nmoves);
                }
            }
            else
            {
                merge_sweep_overlap(states, m_entries, overlap_stats[0].get(),
                                    wr[0].get().get_unchecked(B),
                                    b.get_unchecked(num_vertices(g)), ce, cv,
                                    vmap, label.get_unchecked(B), vlist,
                                    deg_corr, dense, multigraph, g,
                                    random_move, confine_layers, nmerges, niter,
                                    B, rng, S, nmoves);
            }
        }
    }
};

template <class Type>
vector<Type> from_list(boost::python::object list)
{
    vector<Type> v;
    for (int i = 0; i < boost::python::len(list); ++i)
        v.push_back(boost::python::extract<Type>(list[i])());
    return v;
};

template <class Type>
vector<std::reference_wrapper<Type>> from_rlist(boost::python::object list)
{
    vector<std::reference_wrapper<Type>> v;
    for (int i = 0; i < boost::python::len(list); ++i)
        v.emplace_back(boost::python::extract<Type&>(list[i])());
    return v;
};

template <class Type>
vector<std::reference_wrapper<Type>> from_any_list(boost::python::object list)
{
    vector<std::reference_wrapper<Type>> v;
    for (int i = 0; i < boost::python::len(list); ++i)
        v.emplace_back(any_cast<Type&>(boost::python::extract<boost::any&>(list[i])()));
    return v;
};

boost::python::object do_cov_move_sweep(GraphInterface& gi,
                                        boost::any& oce,
                                        boost::any& ocv,
                                        boost::any& ovmap,
                                        boost::python::object& ogis,
                                        boost::python::object& obgi,
                                        boost::python::object& oemat,
                                        boost::python::object& osampler,
                                        boost::python::object& ocavity_sampler,
                                        boost::python::object& omrs,
                                        boost::python::object& omrp,
                                        boost::python::object& omrm,
                                        boost::python::object& owr,
                                        boost::any& ob,
                                        boost::python::object& obs,
                                        bmap_t& bmap,
                                        boost::python::object& obrmap,
                                        boost::python::object& ofree_blocks,
                                        boost::python::object& omaster,
                                        boost::python::object& oslave,
                                        boost::any& olabel, vector<int>& vlist,
                                        bool deg_corr, bool dense,
                                        bool multigraph,
                                        boost::python::object& oeweight,
                                        boost::python::object& ovweight,
                                        boost::python::object& oegroups,
                                        boost::python::object& oesrcpos,
                                        boost::python::object& oetgtpos,
                                        double beta, bool sequential,
                                        bool parallel, bool random_move,
                                        boost::python::object onode_coherent,
                                        double c, bool weighted, size_t nmerges,
                                        boost::any omerge_map,
                                        size_t niter,
                                        size_t B,
                                        boost::python::object& opartition_stats,
                                        boost::python::object& ooverlap_partition_stats,
                                        boost::python::object& ooverlap_stats,
                                        bool verbose, rng_t& rng)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        vmap_t;
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::edge_index_map_t>::type
        emap_t;
    typedef property_map_type::apply<vector<int32_t>,
                                     GraphInterface::vertex_index_map_t>::type
        vvmap_t;
    auto mrs = from_any_list<emap_t>(omrs);
    auto mrp = from_any_list<vmap_t>(omrp);
    auto mrm = from_any_list<vmap_t>(omrm);
    auto wr = from_any_list<vmap_t>(owr);
    vmap_t b = any_cast<vmap_t>(ob);
    auto bs = from_any_list<vmap_t>(obs);

    auto brmap = from_any_list<vmap_t>(obrmap);
    auto free_blocks = from_rlist<vector<size_t>>(ofree_blocks);

    auto eweight = from_any_list<emap_t>(oeweight);
    auto vweight = from_any_list<vmap_t>(ovweight);
    auto egroups = from_rlist<boost::any>(oegroups);
    auto esrcpos = from_any_list<emap_t>(oesrcpos);
    auto etgtpos = from_any_list<emap_t>(oetgtpos);

    vmap_t label = any_cast<vmap_t>(olabel);
    emap_t ce = any_cast<emap_t>(oce);
    vvmap_t cv = any_cast<vvmap_t>(ocv);
    vvmap_t vmap = any_cast<vvmap_t>(ovmap);

    double S = 0;
    size_t nmoves = 0;

    vmap_t merge_map = any_cast<vmap_t>(omerge_map);

    auto gis = from_rlist<GraphInterface>(ogis);

    vector<size_t> eidx;
    for (GraphInterface& g : gis)
        eidx.push_back(g.GetMaxEdgeIndex());

    auto bgi = from_rlist<GraphInterface>(obgi);

    auto partition_stats = from_rlist<partition_stats_t>(opartition_stats);
    auto overlap_partition_stats = from_rlist<overlap_partition_stats_t>(ooverlap_partition_stats);
    auto overlap_stats = from_rlist<overlap_stats_t>(ooverlap_stats);

    auto emat = from_rlist<boost::any>(oemat);

    auto sampler = from_rlist<boost::any>(osampler);
    auto cavity_sampler = from_rlist<boost::any>(ocavity_sampler);

    auto master = from_list<bool>(omaster);
    auto slave = from_list<bool>(oslave);


    bool node_coherent = python::extract<bool>(onode_coherent[0]);
    bool confine_layers = python::extract<bool>(onode_coherent[1]);

    run_action<graph_tool::detail::all_graph_views, boost::mpl::true_>()
        (gi, std::bind(cov_move_sweep_dispatch<emap_t, vmap_t, vvmap_t, emap_t, bmap_t>
                       (ce, cv, vmap, eweight, vweight, egroups, esrcpos, etgtpos,
                        label, vlist, deg_corr, dense, multigraph, beta,
                        sequential, parallel, random_move, node_coherent,
                        confine_layers, c, verbose,
                        gi.GetMaxEdgeIndex(), eidx, nmerges, niter, merge_map,
                        partition_stats, overlap_partition_stats, overlap_stats,
                        master, slave, rng, S, nmoves, bgi, bmap, brmap, free_blocks, B),
                       std::ref(mrs), std::ref(mrp), std::ref(mrm), std::ref(wr),
                       std::ref(b), std::ref(bs), std::ref(bgi), placeholders::_1, std::ref(gis),
                       std::ref(emat), std::ref(sampler), std::ref(cavity_sampler), weighted))();

    return boost::python::make_tuple(S, nmoves);
}

struct covariate_entropy
{
    template <class Graph, class Emap>
    void operator()(Graph& bg, Emap mrs, double& S) const
    {
        for (auto e : edges_range(bg))
            S -= lgamma_fast(mrs[e] + 1);
    }
};

double do_covariate_entropy(GraphInterface& gi, boost::any omrs)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::edge_index_map_t>::type
        emap_t;
    emap_t mrs = any_cast<emap_t>(omrs);

    double S = 0;
    run_action<>()
        (gi, std::bind(covariate_entropy(),
                       placeholders::_1, mrs, std::ref(S)))();
    return S;
}

// void do_create_echash(GraphInterface& gi, size_t l, size_t L,
//                       boost::any& emat_orig, boost::any& emat)
// {
//     run_action<>()(gi, std::bind<void>(create_echash(), placeholders::_1, l, L,
//                                        std::ref(emat_orig), std::ref(emat)))();
// }

void do_ec_hist(GraphInterface& gi, boost::any& aevc, boost::any& aec)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::edge_index_map_t>::type
        emap_t;
    typename emap_t::unchecked_t ec =
        any_cast<emap_t&>(aec).get_unchecked(gi.GetMaxEdgeIndex());
    run_action<>()(gi, std::bind<void>(ec_hist(), placeholders::_1,
                                       placeholders::_2, std::ref(ec)),
                   edge_properties())(aevc);
}


void do_split_graph(GraphInterface& gi, boost::any& aec, boost::any& ab,
                    boost::any& aeweight, boost::any& avweight, boost::any& avc,
                    boost::any& avmap, boost::python::object& ous,
                    boost::python::object& oub,
                    boost::python::object& oueweight,
                    boost::python::object& ouvweight,
                    bmap_t& bmap,
                    boost::python::object& obrmap,
                    boost::python::object& ouvmap)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        vmap_t;
    typedef property_map_type::apply<vector<int32_t>,
                                     GraphInterface::vertex_index_map_t>::type
        vvmap_t;
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::edge_index_map_t>::type
        emap_t;

    emap_t& ec = any_cast<emap_t&>(aec);
    vmap_t& b = any_cast<vmap_t&>(ab);
    vmap_t& vweight = any_cast<vmap_t&>(avweight);
    emap_t& eweight = any_cast<emap_t&>(aeweight);
    vvmap_t& vc = any_cast<vvmap_t&>(avc);
    vvmap_t& vmap = any_cast<vvmap_t&>(avmap);

    auto us = from_rlist<GraphInterface>(ous);
    auto ub = from_any_list<vmap_t>(oub);
    auto uvweight = from_any_list<vmap_t>(ouvweight);
    auto ueweight = from_any_list<emap_t>(oueweight);

    auto rbmap = from_any_list<vmap_t>(obrmap);
    auto uvmap = from_any_list<vmap_t>(ouvmap);

    run_action<>()(gi, std::bind<void>(split_graph(), placeholders::_1,
                                       std::ref(ec), std::ref(b), std::ref(eweight),
                                       std::ref(vweight), std::ref(vc),
                                       std::ref(vmap), std::ref(us),
                                       std::ref(ub),
                                       std::ref(uvweight), std::ref(ueweight),
                                       std::ref(bmap), std::ref(rbmap),
                                       std::ref(uvmap)))();
}

bool bmap_has(const bmap_t& bmap, size_t c, size_t r)
{
    if (c > bmap.size())
        throw GraphException("invalid covariate value:" + lexical_cast<string>(c));
    auto iter = bmap[c].find(r);
    if (iter == bmap[c].end())
        return false;
    return true;
}

size_t bmap_get(const bmap_t& bmap, size_t c, size_t r)
{
    if (c > bmap.size())
        throw GraphException("invalid covariate value:" + lexical_cast<string>(c));
    auto iter = bmap[c].find(r);
    if (iter == bmap[c].end())
        throw GraphException("no mapping for block " + lexical_cast<string>(r)
                             + " in layer " + lexical_cast<string>(c));
    return iter->second;
}

void bmap_set(bmap_t& bmap, size_t c, size_t r, size_t r_u)
{
    if (c > bmap.size())
        throw GraphException("invalid covariate value:" + lexical_cast<string>(c));
    bmap[c][r] = r_u;
}

void bmap_del_c(bmap_t& bmap, size_t c)
{
    if (c > bmap.size())
        throw GraphException("invalid covariate value:" + lexical_cast<string>(c));
    bmap.erase(bmap.begin() + c);
}

bmap_t bmap_copy(const bmap_t& bmap)
{
    return bmap;
}


void export_blockmodel_covariate()
{
    boost::python::class_<bmap_t>("bmap_t")
        .def("has", bmap_has)
        .def("get", bmap_get)
        .def("set", bmap_set)
        .def("del_c", bmap_del_c)
        .def("copy", bmap_copy);

    boost::python::def("cov_move_sweep", do_cov_move_sweep);
    boost::python::def("covariate_entropy", do_covariate_entropy);
    // boost::python::def("create_echash", do_create_echash);
    boost::python::def("ec_hist", do_ec_hist);
    boost::python::def("split_graph", do_split_graph);
}
