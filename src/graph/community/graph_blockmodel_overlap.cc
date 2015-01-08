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
#include "graph_blockmodel.hh"
#include "graph_blockmodel_overlap.hh"


using namespace boost;
using namespace graph_tool;

namespace graph_tool
{
template <class Eprop, class Vprop, class VEprop>
struct move_sweep_overlap_dispatch
{
    move_sweep_overlap_dispatch(Eprop eweight, Vprop vweight,
                                boost::any egroups, VEprop esrcpos,
                                VEprop etgtpos, Vprop label, vector<int>& vlist,
                                bool deg_corr, bool dense, bool multigraph,
                                bool parallel_edges, double beta,
                                bool sequential, bool parallel,
                                bool random_move, double c, bool node_coherent,
                                bool verbose, size_t max_edge_index,
                                size_t nmerges, size_t ntries, Vprop merge_map,
                                overlap_stats_t& overlap_stats,
                                overlap_partition_stats_t& partition_stats,
                                rng_t& rng, double& S, size_t& nmoves,
                                GraphInterface& bgi)

        : eweight(eweight), vweight(vweight), oegroups(egroups), esrcpos(esrcpos),
          etgtpos(etgtpos), label(label), vlist(vlist),
          deg_corr(deg_corr), dense(dense), multigraph(multigraph),
          parallel_edges(parallel_edges), beta(beta), sequential(sequential),
          parallel(parallel), random_move(random_move), c(c),
          node_coherent(node_coherent), verbose(verbose),
          max_edge_index(max_edge_index), nmerges(nmerges), ntries(ntries),
          merge_map(merge_map), overlap_stats(overlap_stats),
          partition_stats(partition_stats), rng(rng), S(S), nmoves(nmoves), bgi(bgi)
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
    bool parallel_edges;
    double beta;
    bool sequential;
    bool parallel;
    bool random_move;
    double c;
    bool node_coherent;
    bool verbose;
    size_t max_edge_index;
    size_t nmerges;
    size_t ntries;
    Vprop merge_map;
    overlap_stats_t& overlap_stats;
    overlap_partition_stats_t& partition_stats;
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
                  boost::any cavity_sampler, BGraph& bg, bool weighted) const
    {
        if (weighted)
        {
            typedef typename property_map_type::apply<DynamicSampler<std::tuple<typename graph_traits<Graph>::edge_descriptor, bool> >,
                                                      GraphInterface::vertex_index_map_t>::type vemap_t;
            vemap_t egroups = any_cast<vemap_t>(oegroups);
            dispatch(mrs, mrp, mrm, wr, b, g, aemat, asampler, cavity_sampler,
                     bg, egroups);
        }
        else
        {
            typedef typename property_map_type::apply<vector<std::tuple<typename graph_traits<Graph>::edge_descriptor, bool> >,
                                                      GraphInterface::vertex_index_map_t>::type vemap_t;
            vemap_t egroups = any_cast<vemap_t>(oegroups);
            dispatch(mrs, mrp, mrm, wr, b, g, aemat, asampler, cavity_sampler,
                     bg, egroups);
        }
    }

    template <class Graph, class BGraph, class Egroups>
    void dispatch(Eprop mrs, Vprop mrp, Vprop mrm, Vprop wr, Vprop b, Graph& g,
                  boost::any& aemat, boost::any asampler, boost::any acavity_sampler,
                  BGraph& bg, Egroups egroups) const
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

            if (nmerges == 0)
            {
                if (!node_coherent)
                {
                    move_sweep_overlap(mrs.get_unchecked(max_BE),
                                       mrp.get_unchecked(num_vertices(bg)),
                                       mrm.get_unchecked(num_vertices(bg)),
                                       wr.get_unchecked(num_vertices(bg)),
                                       b.get_unchecked(num_vertices(g)),
                                       label.get_unchecked(num_vertices(bg)),
                                       vlist, deg_corr, dense, multigraph,
                                       parallel_edges, beta,
                                       eweight.get_unchecked(max_edge_index),
                                       vweight.get_unchecked(num_vertices(g)),
                                       egroups.get_unchecked(num_vertices(bg)),
                                       esrcpos.get_unchecked(eidx),
                                       etgtpos.get_unchecked(eidx), g, bg, emat,
                                       sequential, parallel, random_move, c,
                                       overlap_stats, partition_stats, verbose,
                                       rng, S, nmoves);
                }
                else
                {
                    coherent_move_sweep_overlap(mrs.get_unchecked(max_BE),
                                                mrp.get_unchecked(num_vertices(bg)),
                                                mrm.get_unchecked(num_vertices(bg)),
                                                wr.get_unchecked(num_vertices(bg)),
                                                b.get_unchecked(num_vertices(g)),
                                                label.get_unchecked(num_vertices(bg)), 
                                                vlist, deg_corr, dense, multigraph,
                                                parallel_edges,
                                                eweight.get_unchecked(max_edge_index),
                                                vweight.get_unchecked(num_vertices(g)),
                                                egroups.get_unchecked(num_vertices(bg)),
                                                esrcpos.get_unchecked(eidx),
                                                etgtpos.get_unchecked(eidx), g, bg, emat,
                                                sequential, parallel, random_move, c,
                                                overlap_stats, partition_stats,
                                                verbose, rng, S, nmoves);
                }
            }
            else
            {
                merge_sweep_overlap(mrs.get_unchecked(max_BE),
                                    mrp.get_unchecked(num_vertices(bg)),
                                    mrm.get_unchecked(num_vertices(bg)),
                                    wr.get_unchecked(num_vertices(bg)),
                                    b.get_unchecked(num_vertices(g)),
                                    label.get_unchecked(num_vertices(bg)),
                                    vlist, deg_corr, dense, multigraph,
                                    parallel_edges,
                                    eweight.get_unchecked(max_edge_index),
                                    vweight.get_unchecked(num_vertices(g)),
                                    egroups.get_unchecked(num_vertices(bg)),
                                    esrcpos.get_unchecked(eidx),
                                    etgtpos.get_unchecked(eidx), g, bg, emat,
                                    sampler, cavity_sampler, node_coherent,
                                    sequential, parallel, random_move, c,
                                    nmerges, ntries,
                                    merge_map.get_unchecked(num_vertices(g)),
                                    overlap_stats, partition_stats, verbose,
                                    rng, S, nmoves);
            }
        }
        catch (bad_any_cast&)
        {
            typedef typename get_ehash_t::apply<BGraph>::type emat_t;
            emat_t& emat = any_cast<emat_t&>(aemat);
            if (nmerges == 0)
            {
                move_sweep_overlap(mrs.get_unchecked(num_edges(g)),
                                   mrp.get_unchecked(num_vertices(bg)),
                                   mrm.get_unchecked(num_vertices(bg)),
                                   wr.get_unchecked(num_vertices(bg)),
                                   b.get_unchecked(num_vertices(g)),
                                   label.get_unchecked(num_vertices(bg)), vlist,
                                   deg_corr, dense, multigraph, parallel_edges,
                                   beta, eweight.get_unchecked(max_edge_index),
                                   vweight.get_unchecked(num_vertices(g)),
                                   egroups.get_unchecked(num_vertices(bg)),
                                   esrcpos.get_unchecked(eidx),
                                   etgtpos.get_unchecked(eidx), g, bg, emat,
                                   sequential, parallel, random_move, c,
                                   overlap_stats, partition_stats, verbose, rng,
                                   S, nmoves);
            }
            else
            {
                merge_sweep_overlap(mrs.get_unchecked(num_edges(g)),
                                    mrp.get_unchecked(num_vertices(bg)),
                                    mrm.get_unchecked(num_vertices(bg)),
                                    wr.get_unchecked(num_vertices(bg)),
                                    b.get_unchecked(num_vertices(g)),
                                    label.get_unchecked(num_vertices(bg)),
                                    vlist, deg_corr, dense, multigraph,
                                    parallel_edges,
                                    eweight.get_unchecked(max_edge_index),
                                    vweight.get_unchecked(num_vertices(g)),
                                    egroups.get_unchecked(num_vertices(bg)),
                                    esrcpos.get_unchecked(eidx),
                                    etgtpos.get_unchecked(eidx), g, bg, emat,
                                    sampler, cavity_sampler, node_coherent,
                                    sequential, parallel, random_move, c,
                                    nmerges, ntries,
                                    merge_map.get_unchecked(num_vertices(g)),
                                    overlap_stats, partition_stats, verbose,
                                    rng, S, nmoves);
            }
        }
    }
};


boost::python::object
do_move_sweep_overlap(GraphInterface& gi, GraphInterface& bgi, boost::any& emat,
                      boost::any sampler, boost::any cavity_sampler,
                      boost::any omrs, boost::any omrp, boost::any omrm,
                      boost::any owr, boost::any ob, boost::any olabel,
                      vector<int>& vlist, bool deg_corr, bool dense,
                      bool multigraph, bool parallel_edges, boost::any oeweight,
                      boost::any ovweight, boost::any oegroups,
                      boost::any oesrcpos, boost::any oetgtpos, double beta,
                      bool sequential, bool parallel, bool random_move,
                      double c, bool node_coherent, bool weighted,
                      size_t nmerges, size_t ntries, boost::any omerge_map,
                      overlap_stats_t& overlap_stats,
                      overlap_partition_stats_t& partition_stats, bool verbose,
                      rng_t& rng)
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
        (gi, std::bind(move_sweep_overlap_dispatch<emap_t, vmap_t, vemap_t>
                       (eweight, vweight, oegroups, esrcpos, etgtpos,
                        label, vlist, deg_corr, dense, multigraph, parallel_edges,
                        beta, sequential, parallel, random_move, c, node_coherent,
                        verbose, gi.GetMaxEdgeIndex(), nmerges, ntries, merge_map,
                        overlap_stats, partition_stats, rng, S, nmoves, bgi),
                       mrs, mrp, mrm, wr, b, placeholders::_1,
                       std::ref(emat), sampler, cavity_sampler, weighted))();
    return boost::python::make_tuple(S, nmoves);
}

struct get_overlap_stats
{
    template <class Graph, class Vprop, class VVprop>
    void operator()(Graph& g, Vprop b, VVprop half_edges, Vprop node_index,
                    size_t B, overlap_stats_t& overlap_stats) const
    {
        overlap_stats =  overlap_stats_t(g, b, half_edges, node_index, B);
    }
};

overlap_stats_t
do_get_overlap_stats(GraphInterface& gi, boost::any ob,
                     boost::any ohalf_edges, boost::any onode_index,
                     size_t B)
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

    overlap_stats_t overlap_stats;

    vmap_t b = any_cast<vmap_t>(ob);
    vvmap_t half_edges = any_cast<vvmap_t>(ohalf_edges);
    vmap_t node_index = any_cast<vmap_t>(onode_index);

    run_action<>()(gi, std::bind(get_overlap_stats(), placeholders::_1, b,
                                 half_edges, node_index, B,
                                 std::ref(overlap_stats)))();
    return overlap_stats;
}

struct get_overlap_partition_stats
{
    template <class Graph, class Vprop, class Eprop>
    void operator()(Graph& g, Vprop b, Eprop eweight, size_t N, size_t B,
                    overlap_stats_t& overlap_stats,
                    overlap_partition_stats_t& partition_stats) const
    {
        partition_stats = overlap_partition_stats_t(g, b, overlap_stats,
                                                    eweight, N, B);
    }
};

overlap_partition_stats_t
do_get_overlap_partition_stats(GraphInterface& gi, boost::any ob,
                               boost::any aeweight, size_t N, size_t B,
                               overlap_stats_t& overlap_stats)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        vmap_t;
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::edge_index_map_t>::type
        emap_t;

    overlap_partition_stats_t partition_stats;

    vmap_t b = any_cast<vmap_t>(ob);
    emap_t eweight = any_cast<emap_t>(aeweight);

    run_action<>()(gi, std::bind(get_overlap_partition_stats(),
                                 placeholders::_1, b, eweight, N, B,
                                 std::ref(overlap_stats),
                                 std::ref(partition_stats)))();

    return partition_stats;
}

double do_get_overlap_parallel_entropy(GraphInterface& gi, boost::any ob,
                                       overlap_stats_t& overlap_stats)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        vmap_t;

    vmap_t b = any_cast<vmap_t>(ob);

    double S = 0;
    run_action<>()
        (gi, std::bind(entropy_parallel_edges_overlap(),
                       placeholders::_1, b, std::ref(overlap_stats),
                       std::ref(S)))();
    return S;
}


struct get_eg_overlap
{
    template <class Graph, class EGraph, class EVprop, class VProp,
              class VVProp>
    void operator()(Graph& g, EGraph& eg, EVprop be, VProp b, VProp node_index,
                    VVProp half_edges) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        for (auto e : edges_range(g))
        {
            vertex_t s = get_source(e, g);
            vertex_t t = get_target(e, g);
            vertex_t u = add_vertex(eg);
            vertex_t v = add_vertex(eg);
            add_edge(u, v, eg);
            if (be[e].size() != 2)
                throw GraphException("Edge block property map must have two values per edge");
            b[u] = be[e][0];
            b[v] = be[e][1];
            node_index[u] = s;
            node_index[v] = t;
            half_edges[s].push_back(u);
            half_edges[t].push_back(v);
        }
    }
};

void do_get_eg_overlap(GraphInterface& gi, GraphInterface& egi, boost::any obe,
                       boost::any ob, boost::any onode_index,
                       boost::any ohalf_edges)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        vmap_t;
    typedef property_map_type::apply<vector<int32_t>,
                                     GraphInterface::vertex_index_map_t>::type
        vvmap_t;
    typedef property_map_type::apply<vector<int32_t>,
                                     GraphInterface::edge_index_map_t>::type
        evmap_t;

    vmap_t b = any_cast<vmap_t>(ob);
    evmap_t be = any_cast<evmap_t>(obe);
    vmap_t node_index = any_cast<vmap_t>(onode_index);
    vvmap_t half_edges = any_cast<vvmap_t>(ohalf_edges);

    run_action<>()(gi, std::bind(get_eg_overlap(), placeholders::_1,
                                 std::ref(egi.GetGraph()), be, b, node_index,
                                 half_edges))();
}

struct get_be_overlap
{
    template <class Graph, class EGraph, class EVprop, class VProp>
    void operator()(Graph& g, EGraph& eg, EVprop be, VProp b, VProp node_index)
        const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        for (auto ei : edges_range(eg))
        {
            vertex_t u = source(ei, eg);
            vertex_t v = target(ei, eg);

            size_t s = node_index[u];
            size_t t = node_index[v];

            for (auto e : out_edges_range(s, g))
            {
                if (!be[e].empty() || target(e, g) != t)
                    continue;
                if (is_directed::apply<Graph>::type::value || s < target(e, g))
                    be[e] = {b[u], b[v]};
                else
                    be[e] = {b[v], b[u]};
            }

            for (auto e : in_edges_range(t, g))
            {
                if (!be[e].empty() || source(e, g) != s)
                    continue;
                be[e] = {b[u], b[v]};
            }
        }
    }
};

void do_get_be_overlap(GraphInterface& gi, GraphInterface& egi, boost::any obe,
                       boost::any ob, boost::any onode_index)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        vmap_t;
    typedef property_map_type::apply<vector<int32_t>,
                                     GraphInterface::edge_index_map_t>::type
        evmap_t;
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::edge_index_map_t>::type
        emap_t;

    vmap_t b = any_cast<vmap_t>(ob);
    evmap_t be = any_cast<evmap_t>(obe);
    vmap_t node_index = any_cast<vmap_t>(onode_index);

    run_action<>()(gi, std::bind(get_be_overlap(), placeholders::_1,
                                 std::ref(egi.GetGraph()), be, b,
                                 node_index))();
}


// struct get_overlap_proj
// {
//     template <class Graph, class BGraph, class Eprop, class VProp>
//     void operator()(Graph& g, BGraph& bg, Eprop mrs, VProp b) const
//     {
//         unordered_map<pair<size_t,size_t>, pair<size_t,size_t>,
//                       boost::hash<pair<size_t,size_t> > > edges;

//         size_t pos = 0;
//         typename graph_traits<BGraph>::edge_iterator ei, ei_end;
//         for (tie(ei, ei_end) = boost::edges(bg); ei != ei_end; ++ei)
//         {
//             if (mrs[*ei] == 0)
//                 continue;
//             size_t r = source(*ei, bg);
//             size_t s = target(*ei, bg);
//             if (!is_directed::apply<Graph>::type::value && r > s)
//                 std::swap(r, s);
//             edges[make_pair(r, s)] = make_pair(pos, pos + 1);
//             pos += 2;
//         }

//         typename graph_traits<Graph>::edge_iterator e, e_end;
//         for (tie(e, e_end) = boost::edges(g); e != e_end; ++e)
//         {
//             size_t u = get_source(*e, g);
//             size_t v = get_target(*e, g);
//             size_t r = b[u];
//             size_t s = b[v];
//             if (!is_directed::apply<Graph>::type::value && r > s)
//                 std::swap(r, s);
//             tie(r, s) = edges[make_pair(r, s)];
//             if (!is_directed::apply<Graph>::type::value && b[u] > b[v])
//                 std::swap(r, s);
//             b[u] = r;
//             b[v] = s;
//         }
//     }
// };

// void do_get_overlap_proj(GraphInterface& gi, GraphInterface& bgi,
//                          boost::any omrs, boost::any ob)
// {
//     typedef property_map_type::apply<int32_t,
//                                      GraphInterface::vertex_index_map_t>::type
//         vmap_t;
//     typedef property_map_type::apply<int32_t,
//                                      GraphInterface::edge_index_map_t>::type
//         emap_t;

//     emap_t mrs = any_cast<emap_t>(omrs);
//     vmap_t b = any_cast<vmap_t>(ob);

//     run_action<>()(gi, std::bind(get_overlap_proj(), placeholders::_1,
//                                  std::ref(bgi.GetGraph()), mrs, b))();
// }

struct get_be_from_b_overlap
{
    template <class Graph, class EVprop, class VProp>
    void operator()(Graph& g, EVprop be, VProp b)
        const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        typename graph_traits<Graph>::edge_iterator ei, ei_end;
        for (auto e : edges_range(g))
        {
            vertex_t s = get_source(e, g);
            vertex_t t = get_target(e, g);
            be[e] = {b[s], b[t]};
        }
    }
};

void do_get_be_from_b_overlap(GraphInterface& gi, boost::any obe, boost::any ob)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        vmap_t;
    typedef property_map_type::apply<vector<int32_t>,
                                     GraphInterface::edge_index_map_t>::type
        evmap_t;

    vmap_t b = any_cast<vmap_t>(ob);
    evmap_t be = any_cast<evmap_t>(obe);

    run_action<>()(gi, std::bind(get_be_from_b_overlap(), placeholders::_1,
                                 be, b))();
}


struct get_bv_overlap
{
    template <class Graph, class VProp, class VVProp>
    void operator()(Graph& g, VProp b, VProp node_index, VVProp bv, VVProp bc_in,
                    VVProp bc_out, VVProp bc_total) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

        typedef unordered_map<int, int> map_t;
        vector<map_t> hist_in;
        vector<map_t> hist_out;

        for (auto v : vertices_range(g))
        {
            if (out_degree(v, g) > 0)
            {
                vertex_t s = node_index[v];
                if (s >= hist_out.size())
                    hist_out.resize(s + 1);
                hist_out[s][b[v]]++;
            }

            if (in_degreeS()(v, g) > 0)
            {
                vertex_t t = node_index[v];
                if (t >= hist_in.size())
                    hist_in.resize(t + 1);
                hist_in[t][b[v]]++;
            }
        }

        size_t N = max(hist_in.size(), hist_out.size());
        hist_in.resize(N);
        hist_out.resize(N);

        set<size_t> rs;
        for (size_t i = 0; i < N; ++i)
        {
            rs.clear();
            for (auto iter = hist_out[i].begin(); iter != hist_out[i].end(); ++iter)
                rs.insert(iter->first);
            for (auto iter = hist_in[i].begin(); iter != hist_in[i].end(); ++iter)
                rs.insert(iter->first);
            // if (rs.empty())
            //     throw GraphException("Cannot have empty overlapping block membership!");
            for (size_t r : rs)
            {
                bv[i].push_back(r);

                auto iter_in = hist_in[i].find(r);
                if (iter_in != hist_in[i].end())
                    bc_in[i].push_back(iter_in->second);
                else
                    bc_in[i].push_back(0);

                auto iter_out = hist_out[i].find(r);
                if (iter_out != hist_out[i].end())
                    bc_out[i].push_back(iter_out->second);
                else
                    bc_out[i].push_back(0);

                bc_total[i].push_back(bc_in[i].back() +
                                      bc_out[i].back());
            }
        }
    }
};

void do_get_bv_overlap(GraphInterface& gi, boost::any ob,  boost::any onode_index,
                       boost::any obv, boost::any obc_in, boost::any obc_out,
                       boost::any obc_total)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        vmap_t;
    typedef property_map_type::apply<vector<int32_t>,
                                     GraphInterface::vertex_index_map_t>::type
        vvmap_t;
    typedef property_map_type::apply<vector<int32_t>,
                                     GraphInterface::edge_index_map_t>::type
        evmap_t;

    vmap_t b = any_cast<vmap_t>(ob);
    vmap_t node_index = any_cast<vmap_t>(onode_index);
    vvmap_t bv = any_cast<vvmap_t>(obv);
    vvmap_t bc_in = any_cast<vvmap_t>(obc_in);
    vvmap_t bc_out = any_cast<vvmap_t>(obc_out);
    vvmap_t bc_total = any_cast<vvmap_t>(obc_total);

    run_action<>()(gi, std::bind(get_bv_overlap(), placeholders::_1, b,
                                 node_index, bv, bc_in, bc_out, bc_total))();
 }

struct get_wr_overlap
{
    template <class Graph, class VProp, class VVProp>
    void operator()(Graph& g, VVProp bv, VProp wr) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

        for (auto v : vertices_range(g))
        {
            for (size_t i = 0; i < bv[v].size(); ++i)
                wr[bv[v][i]]++;
        }
    }
};

void do_get_wr_overlap(GraphInterface& gi, boost::any obv,
                       boost::any owr)
{
    typedef property_map_type::apply<vector<int32_t>,
                                     GraphInterface::vertex_index_map_t>::type
        vvmap_t;
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        vmap_t;

    vvmap_t bv = any_cast<vvmap_t>(obv);
    vmap_t wr = any_cast<vmap_t>(owr);

    run_action<>()(gi, std::bind(get_wr_overlap(), placeholders::_1, bv, wr))();
}

struct get_nodeset_overlap
{
    template <class Graph, class VProp, class VVProp>
    void operator()(Graph& g, VProp node_index, VVProp half_edges)
        const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        for (auto e : edges_range(g))
        {
            vertex_t s = get_source(e, g);
            vertex_t t = get_target(e, g);
            half_edges[node_index[s]].push_back(s);
            half_edges[node_index[t]].push_back(t);
        }
    }
};

void do_get_nodeset_overlap(GraphInterface& gi, boost::any onode_index,
                            boost::any ohalf_edges)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        vmap_t;
    typedef property_map_type::apply<vector<int32_t>,
                                     GraphInterface::vertex_index_map_t>::type
        vvmap_t;

    vmap_t node_index = any_cast<vmap_t>(onode_index);
    vvmap_t half_edges = any_cast<vvmap_t>(ohalf_edges);

    run_action<>()(gi, std::bind(get_nodeset_overlap(), placeholders::_1,
                                 node_index, half_edges))();
}

struct get_augmented_overlap
{
    template <class Graph, class VProp>
    void operator()(Graph& g, VProp b, VProp node_index, VProp br_map,
                    vector<int32_t>& br_b, vector<int32_t>& br_ni) const
    {

        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

        unordered_map<std::tuple<int, int>, size_t> idx_map;
        vector<std::tuple<int, int>> idx_rmap;
        size_t pos = 0;

        for (auto v : vertices_range(g))
        {
            size_t vi = node_index[v];
            auto br = std::make_tuple(b[v], vi);
            size_t idx;
            auto iter = idx_map.find(br);
            if (iter != idx_map.end())
            {
                idx = iter->second;
            }
            else
            {
                idx = pos;
                idx_map[br] = pos++;
                idx_rmap.push_back(br);
            }
            br_map[v] = idx;
        }

        for (size_t i = 0; i < idx_rmap.size(); ++i)
        {
            auto& br = idx_rmap[i];
            br_b.push_back(get<0>(br));
            br_ni.push_back(get<1>(br));
        }
    }
};

void do_get_augmented_overlap(GraphInterface& gi, boost::any ob,
                              boost::any onode_index, boost::any obr_map,
                              vector<int32_t>& br_b, vector<int32_t>& br_ni)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        vmap_t;

    vmap_t b = any_cast<vmap_t>(ob);
    vmap_t node_index = any_cast<vmap_t>(onode_index);
    vmap_t br_map = any_cast<vmap_t>(obr_map);

    run_action<>()(gi, std::bind(get_augmented_overlap(), placeholders::_1,
                                 b, node_index, br_map, std::ref(br_b),
                                 std::ref(br_ni)))();
}


struct get_overlap_split
{
    template <class Graph, class VVProp, class VProp>
    void operator()(Graph& g, VVProp bv, VProp b) const
    {
        unordered_map<vector<int>, size_t> bvset;

        for (auto v : vertices_range(g))
        {
            auto r = bv[v];
            auto iter = bvset.find(r);
            if (iter == bvset.end())
                iter = bvset.insert(make_pair(r, bvset.size())).first;
            b[v] = iter->second;
        }
    }
};

void do_get_overlap_split(GraphInterface& gi, boost::any obv, boost::any ob)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        vmap_t;
    typedef property_map_type::apply<vector<int32_t>,
                                     GraphInterface::vertex_index_map_t>::type
        vvmap_t;

    vvmap_t bv = any_cast<vvmap_t>(obv);
    vmap_t b = any_cast<vmap_t>(ob);

    run_action<>()(gi, std::bind(get_overlap_split(),
                                 placeholders::_1, bv, b))();
}

} // namespace graph_tool


void export_blockmodel_overlap()
{
    using namespace boost::python;

    class_<overlap_stats_t>("overlap_stats")
        .def("is_enabled", &overlap_stats_t::is_enabled);
    class_<overlap_partition_stats_t>("overlap_partition_stats")
        .def("is_enabled", &overlap_partition_stats_t::is_enabled)
        .def("get_partition_dl", &overlap_partition_stats_t::get_partition_dl)
        .def("get_deg_dl", &overlap_partition_stats_t::get_deg_dl);

    def("move_sweep_overlap", do_move_sweep_overlap);
    def("init_overlap_stats", do_get_overlap_stats);
    def("init_overlap_partition_stats", do_get_overlap_partition_stats);

    def("overlap_parallel_entropy", do_get_overlap_parallel_entropy);

    // def("get_overlap_proj", do_get_overlap_proj);

    def("get_eg_overlap", do_get_eg_overlap);
    def("get_be_overlap", do_get_be_overlap);
    def("get_be_from_b_overlap", do_get_be_from_b_overlap);
    def("get_bv_overlap", do_get_bv_overlap);
    def("get_wr_overlap", do_get_wr_overlap);
    def("get_nodeset_overlap", do_get_nodeset_overlap);
    def("get_augmented_overlap", do_get_augmented_overlap);
    def("get_overlap_split", do_get_overlap_split);
}
