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


#define BOOST_PYTHON_MAX_ARITY 30
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

#include "random.hh"

#include "config.h"
#include "graph_blockmodel.hh"

using namespace boost;

// ====================
// Entropy calculation
// ====================

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
        (gi, boost::bind<void>(entropy(), mrs, mrp, mrm, wr, deg_corr, _1,
                               ref(S)))();
    return S;
}

double lbinom(double N, double k)
{
    return lgamma(N + 1) - lgamma(N - k + 1) - lgamma(k + 1);
}

double SB(double B, size_t N, size_t E, bool directed)
{
    double x;
    if (directed)
        x = (B * B);
    else
        x = (B * (B + 1)) / 2;
    return lbinom(x + E - 1, E) / E + N * log(B) / double(E);
}


// ===============================
// Block merge and merge distance
// ===============================

double do_dist(GraphInterface& gi, size_t r, size_t s, boost::any omrs,
               boost::any omrp, boost::any omrm, boost::any owr, bool deg_corr)
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

    double d;

    run_action<graph_tool::detail::all_graph_views, mpl::true_>()
        (gi, boost::bind<void>(dist(), r, s, mrs, mrp, mrm, wr, deg_corr, _1, ref(d)))();

    return d;
}

python::object do_min_dist(GraphInterface& gi, int n, boost::any omrs,
                           boost::any omrp, boost::any omrm, boost::any owr,
                           bool deg_corr, vector<int>& vlist, rng_t& rng)
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

    boost::tuple<double, size_t, size_t> min_d;
    run_action<graph_tool::detail::all_graph_views, mpl::true_>()
        (gi, boost::bind<void>(min_dist<rng_t>(rng), n, mrs, mrp, mrm, wr, deg_corr, _1,
                               ref(vlist), ref(min_d)))();

    return python::make_tuple(get<0>(min_d), get<1>(min_d), get<2>(min_d));
}


void do_b_join(GraphInterface& gi, size_t r, size_t s, boost::any omrs,
             boost::any omrp, boost::any omrm, boost::any owr, bool deg_corr,
             vector<int>& vlist)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        vmap_t;
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::edge_index_map_t>::type
        emap_t;
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        vimap_t;
    emap_t mrs = any_cast<emap_t>(omrs);
    vmap_t mrp = any_cast<vmap_t>(omrp);
    vmap_t mrm = any_cast<vmap_t>(omrm);
    vmap_t wr = any_cast<vmap_t>(owr);

    run_action<graph_tool::detail::all_graph_views, mpl::true_>()
        (gi, boost::bind<void>(b_join(r, s), mrs, mrp, mrm, wr, deg_corr,
                               ref(vlist), _1))();
}

// ===============================
// Block moves
// ===============================

boost::any do_create_emat(GraphInterface& gi, size_t B)
{
    boost::any emat;
    run_action<>()(gi, boost::bind<void>(create_emat(), _1, boost::ref(emat), B))();
    return emat;
}


template <class Eprop, class Vprop>
struct remove_vertex_dispatch
{
    remove_vertex_dispatch(Eprop eweight, Vprop vweight,
                        GraphInterface& bgi)
        : eweight(eweight), vweight(vweight), bgi(bgi) {}

    Eprop eweight;
    Vprop vweight;
    GraphInterface& bgi;

    template <class Graph>
    void operator()(size_t v, Eprop mrs, Vprop mrp, Vprop mrm,
                    Vprop wr, Vprop b, Graph& g, boost::any& emat) const
    {
        if (is_directed::apply<Graph>::type::value)
        {
            dispatch(v, mrs, mrp, mrm, wr, b, g, emat, bgi.GetGraph());
        }
        else
        {
            UndirectedAdaptor<GraphInterface::multigraph_t> ug(bgi.GetGraph());
            dispatch(v, mrs, mrp, mrm, wr, b, g, emat, ug);
        }
    }

    template <class Graph, class BGraph>
    void dispatch(size_t v, Eprop mrs, Vprop mrp, Vprop mrm, Vprop wr, Vprop b,
                  Graph& g, boost::any& aemat, BGraph& bg) const
    {
        typedef typename get_emat_t::apply<BGraph>::type emat_t;
        emat_t& emat = any_cast<emat_t&>(aemat);
        remove_vertex(v, mrs, mrp, mrm, wr, b, eweight, vweight, g, bg, emat);
    }
};



void do_remove_vertex(GraphInterface& gi, GraphInterface& gbi, size_t v,
                      boost::any omrs, boost::any omrp, boost::any omrm,
                      boost::any owr, boost::any ob, boost::any oeweight,
                      boost::any ovweight, boost::any& emat)
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
    vmap_t b = any_cast<vmap_t>(ob);
    emap_t eweight = any_cast<emap_t>(oeweight);
    vmap_t vweight = any_cast<vmap_t>(ovweight);

    run_action<graph_tool::detail::all_graph_views, mpl::true_>()
        (gi, boost::bind<void>(remove_vertex_dispatch<emap_t, vmap_t>(eweight, vweight, gbi),
                               v, mrs, mrp, mrm, wr, b, _1, ref(emat)))();
}

template <class Eprop, class Vprop>
struct add_vertex_dispatch
{
    add_vertex_dispatch(Eprop eweight, Vprop vweight,
                        GraphInterface& bgi)
        : eweight(eweight), vweight(vweight), bgi(bgi) {}

    Eprop eweight;
    Vprop vweight;
    GraphInterface& bgi;

    template <class Graph>
    void operator()(pair<size_t, size_t> vr, Eprop mrs, Vprop mrp, Vprop mrm,
                    Vprop wr, Vprop b, Graph& g, boost::any& emat) const
    {
        if (is_directed::apply<Graph>::type::value)
        {
            dispatch(vr.first, vr.second, mrs, mrp, mrm, wr, b, g, emat,
                     bgi.GetGraph());
        }
        else
        {
            UndirectedAdaptor<GraphInterface::multigraph_t> ug(bgi.GetGraph());
            dispatch(vr.first, vr.second, mrs, mrp, mrm, wr, b, g, emat, ug);
        }
    }

    template <class Graph, class BGraph>
    void dispatch(size_t v, size_t r, Eprop mrs, Vprop mrp, Vprop mrm,
                  Vprop wr, Vprop b, Graph& g, boost::any& aemat, BGraph& bg) const
    {
        typedef typename get_emat_t::apply<BGraph>::type emat_t;
        emat_t& emat = any_cast<emat_t&>(aemat);
        add_vertex(v, r, mrs, mrp, mrm, wr, b, eweight, vweight, g, bg, emat);
    }
};


void do_add_vertex(GraphInterface& gi, GraphInterface& gbi, size_t v, size_t r,
                   boost::any omrs, boost::any omrp, boost::any omrm,
                   boost::any owr, boost::any ob, boost::any ovweight,
                   boost::any oeweight, boost::any& emat)
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
    vmap_t b = any_cast<vmap_t>(ob);
    emap_t eweight = any_cast<emap_t>(oeweight);
    vmap_t vweight = any_cast<vmap_t>(ovweight);

    run_action<graph_tool::detail::all_graph_views, mpl::true_>()
        (gi, boost::bind<void>(add_vertex_dispatch<emap_t, vmap_t>(eweight, vweight, gbi),
                               make_pair(v, r), mrs, mrp, mrm, wr,
                               b, _1, ref(emat)))();
}


template <class Eprop, class Vprop>
struct move_vertex_dispatch
{
    move_vertex_dispatch(Eprop eweight, Vprop vweight, bool deg_corr, double& S,
                         GraphInterface& bgi)
        : eweight(eweight), vweight(vweight), deg_corr(deg_corr), S(S), bgi(bgi) {}

    Eprop eweight;
    Vprop vweight;
    bool deg_corr;
    double& S;
    GraphInterface& bgi;

    template <class Graph>
    void operator()(pair<size_t, size_t> vnr, Eprop mrs, Vprop mrp, Vprop mrm,
                    Vprop wr, Vprop b, Graph& g, boost::any& emat) const
    {
        if (is_directed::apply<Graph>::type::value)
        {
            dispatch(vnr.first, vnr.second, mrs, mrp, mrm, wr, b, g, emat,
                     bgi.GetGraph());
        }
        else
        {
            UndirectedAdaptor<GraphInterface::multigraph_t> ug(bgi.GetGraph());
            dispatch(vnr.first, vnr.second, mrs, mrp, mrm, wr, b, g, emat,
                     ug);
        }
    }

    template <class Graph, class BGraph>
    void dispatch(size_t v, size_t nr, Eprop mrs, Vprop mrp, Vprop mrm,
                  Vprop wr, Vprop b, Graph& g, boost::any& aemat, BGraph& bg) const
    {
        typedef typename get_emat_t::apply<BGraph>::type emat_t;
        emat_t& emat = any_cast<emat_t&>(aemat);
        S = move_vertex(v, nr, mrs, mrp, mrm, wr, b, deg_corr, eweight, vweight,
                        g, bg, emat, true);
    }
};

double do_move_vertex(GraphInterface& gi, GraphInterface& bgi, boost::any& emat,
                      size_t v, size_t nr, boost::any omrs, boost::any omrp,
                      boost::any omrm, boost::any owr, boost::any ob,
                      bool deg_corr, boost::any oeweight, boost::any ovweight)
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
    vmap_t b = any_cast<vmap_t>(ob);
    emap_t eweight = any_cast<emap_t>(oeweight);
    vmap_t vweight = any_cast<vmap_t>(ovweight);

    double S;
    run_action<graph_tool::detail::all_graph_views, mpl::true_>()
        (gi, boost::bind<void>(move_vertex_dispatch<emap_t, vmap_t>(eweight, vweight, deg_corr, S, bgi),
                               make_pair(v, nr), mrs, mrp, mrm, wr,
                               b, _1, ref(emat)))();
    return S;
}


//============
// Main loops
//============


template <class Eprop, class Vprop, class VEprop>
struct move_sweep_dispatch
{
    move_sweep_dispatch(Eprop eweight, Vprop vweight, boost::any egroups,
                        VEprop esrcpos, VEprop etgtpos, Vprop label, size_t L,
                        vector<int>& vlist, bool deg_corr, double beta,
                        bool sequential, bool verbose, size_t max_edge_index,
                        rng_t& rng, double& S, size_t& nmoves,
                        GraphInterface& bgi)

        : eweight(eweight), vweight(vweight), oegroups(egroups), esrcpos(esrcpos),
          etgtpos(etgtpos), label(label), L(L), vlist(vlist),
          deg_corr(deg_corr), beta(beta), verbose(verbose),
          max_edge_index(max_edge_index),
          rng(rng), S(S), nmoves(nmoves), bgi(bgi)
    {}

    Eprop eweight;
    Vprop vweight;
    boost::any oegroups;
    VEprop esrcpos;
    VEprop etgtpos;
    Vprop label;
    size_t L;
    size_t n;
    vector<int>& vlist;
    bool deg_corr;
    double beta;
    bool sequential;
    bool verbose;
    size_t max_edge_index;
    rng_t& rng;
    double& S;
    size_t& nmoves;
    GraphInterface& bgi;

    template <class Graph>
    void operator()(Eprop mrs, Vprop mrp, Vprop mrm, Vprop wr, Vprop b,
                    Graph& g, boost::any& emat) const
    {
        if (is_directed::apply<Graph>::type::value)
        {
            dispatch(mrs, mrp, mrm, wr, b, g, emat, bgi.GetGraph());
        }
        else
        {
            UndirectedAdaptor<GraphInterface::multigraph_t> ug(bgi.GetGraph());
            dispatch(mrs, mrp, mrm, wr, b, g, emat, ug);
        }
    }


    template <class Graph, class BGraph>
    void dispatch(Eprop mrs, Vprop mrp, Vprop mrm, Vprop wr, Vprop b, Graph& g,
                  boost::any& aemat, BGraph& bg) const
    {
        const GraphWrap<BGraph> wbg(bg, bgi);
        typedef typename get_emat_t::apply<BGraph>::type emat_t;
        emat_t& emat = any_cast<emat_t&>(aemat);

        typedef typename property_map_type::apply<vector<boost::tuple<typename graph_traits<Graph>::edge_descriptor, bool, size_t> >,
                                                  GraphInterface::vertex_index_map_t>::type vemap_t;
        vemap_t egroups = any_cast<vemap_t>(oegroups);

        //make sure the properties are _unchecked_, since otherwise it affects
        //performance

        move_sweep(mrs.get_unchecked(emat.num_elements()),
                   mrp.get_unchecked(num_vertices(bg)),
                   mrm.get_unchecked(num_vertices(bg)),
                   wr.get_unchecked(num_vertices(bg)),
                   b.get_unchecked(num_vertices(g)),
                   label.get_unchecked(num_vertices(g)),
                   L, vlist, deg_corr, beta,
                   eweight.get_unchecked(max_edge_index + 1),
                   vweight.get_unchecked(num_vertices(g)),
                   egroups.get_unchecked(num_vertices(bg)),
                   esrcpos.get_unchecked(max_edge_index + 1),
                   etgtpos.get_unchecked(max_edge_index + 1),
                   g, wbg, emat, sequential, verbose, rng, S, nmoves);
    }

};


python::object do_move_sweep(GraphInterface& gi, GraphInterface& bgi,
                             boost::any& emat, boost::any omrs, boost::any omrp,
                             boost::any omrm, boost::any owr, boost::any ob,
                             boost::any olabel, size_t L, vector<int>& vlist,
                             bool deg_corr, boost::any oeweight,
                             boost::any ovweight, boost::any oegroups,
                             boost::any oesrcpos, boost::any oetgtpos,
                             double beta, bool sequential, bool verbose,
                             rng_t& rng)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        vmap_t;
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::edge_index_map_t>::type
        emap_t;
    typedef property_map_type::apply<vector<int32_t>,
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
    run_action<graph_tool::detail::all_graph_views, mpl::true_>()
        (gi, boost::bind<void>(move_sweep_dispatch<emap_t, vmap_t, vemap_t>
                               (eweight, vweight, oegroups, esrcpos, etgtpos,
                                label, L, vlist, deg_corr, beta,
                                sequential, verbose, gi.GetMaxEdgeIndex(),
                                rng, S, nmoves, bgi),
                               mrs, mrp, mrm, wr, b, _1, ref(emat)))();
    return python::make_tuple(S, nmoves);
}

boost::any do_build_egroups(GraphInterface& gi, GraphInterface& bgi,
                            boost::any ob, boost::any oeweights,
                            boost::any oesrcpos, boost::any oetgtpos)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        vmap_t;
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::edge_index_map_t>::type
        emap_t;
    typedef property_map_type::apply<vector<int32_t>,
                                     GraphInterface::edge_index_map_t>::type
        vemap_t;
    vmap_t b = any_cast<vmap_t>(ob);

    vemap_t esrcpos = any_cast<vemap_t>(oesrcpos);
    vemap_t etgtpos = any_cast<vemap_t>(oetgtpos);
    emap_t eweights = any_cast<emap_t>(oeweights);

    boost::any oegroups;
    run_action<graph_tool::detail::all_graph_views, mpl::true_>()
        (gi, boost::bind<void>(build_egroups(), b, ref(oegroups), esrcpos, etgtpos,
                               eweights, _1, bgi.GetVertexIndex()))();
    return oegroups;
}

struct collect_edge_marginals_dispatch
{
    template <class Graph, class Vprop, class MEprop>
    void operator()(Graph& g, size_t B, Vprop cb, MEprop p,
                    tr1::tuple<boost::any, GraphInterface&> abg) const
    {
        typedef typename Graph::orig_graph_t graph_t;
        const GraphWrap<graph_t> bg(*any_cast<graph_t*>(get<0>(abg)), get<1>(abg));
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

    run_action<graph_tool::detail::all_graph_views, mpl::true_>()
        (gi, boost::bind<void>(collect_edge_marginals_dispatch(),
                               _1, B, b, p,
                               tr1::tuple<boost::any, GraphInterface&>(gbi.GetGraphView(), gbi)))();
}

python::tuple do_bethe_entropy(GraphInterface& gi, size_t B, boost::any op,
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
    run_action<graph_tool::detail::all_graph_views, mpl::true_>()
        (gi, boost::bind<void>(bethe_entropy(),
                               _1, B, p, pv, ref(H), ref(sH),
                               ref(Hmf), ref(sHmf)))();
    return python::make_tuple(H, sH, Hmf, sHmf);
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

    run_action<graph_tool::detail::all_graph_views, mpl::true_>()
        (gi, boost::bind<void>(collect_vertex_marginals_dispatch(),
                               _1, b, _2),
         vertex_scalar_vector_properties())(op);
}

vector<int32_t> get_vector(size_t n)
{
    return vector<int32_t>(n);
}

void export_blockmodel()
{
    using namespace boost::python;

    def("get_vector", get_vector);

    def("add_vertex", do_add_vertex);
    def("remove_vertex", do_remove_vertex);
    def("move_vertex", do_move_vertex);

    def("create_emat", do_create_emat);
    def("build_egroups", do_build_egroups);

    def("move_sweep", do_move_sweep);

    def("join", do_b_join);
    def("dist", do_dist);
    def("min_dist", do_min_dist);
    def("entropy", do_get_ent);

    def("edge_marginals", do_collect_edge_marginals);
    def("bethe_entropy", do_bethe_entropy);
    def("vertex_marginals", do_collect_vertex_marginals);

    def("SB", SB);
}
