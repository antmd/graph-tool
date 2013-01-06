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

#include "graph_filtering.hh"
#include "graph.hh"
#include "graph_properties.hh"
#include "graph_exceptions.hh"

#if (GCC_VERSION >= 40400)
#   include <tr1/random>
#else
#   include <boost/tr1/random.hpp>
#endif

#include <boost/lambda/bind.hpp>

#include "graph_sfdp.hh"
#include "random.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;


void sfdp_layout(GraphInterface& g, boost::any pos, boost::any vweight,
                 boost::any eweight, boost::any pin, python::object spring_parms,
                 double theta, double init_step, double step_schedule,
                 size_t max_level, double epsilon, size_t max_iter,
                 bool adaptive, bool verbose)
{
    typedef ConstantPropertyMap<int32_t,GraphInterface::vertex_t> vweight_map_t;
    typedef ConstantPropertyMap<int32_t,GraphInterface::edge_t> eweight_map_t;
    typedef mpl::push_back<vertex_scalar_properties, vweight_map_t>::type
        vertex_props_t;
    typedef mpl::push_back<edge_scalar_properties, eweight_map_t>::type
        edge_props_t;

    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        group_map_t;

    double C = python::extract<double>(spring_parms[0]);
    double K = python::extract<double>(spring_parms[1]);
    double p = python::extract<double>(spring_parms[2]);
    double gamma = python::extract<double>(spring_parms[3]);
    double mu = python::extract<double>(spring_parms[4]);
    double mu_p = python::extract<double>(spring_parms[5]);
    group_map_t groups =
        any_cast<group_map_t>(python::extract<any>(spring_parms[6]));

    if(vweight.empty())
        vweight = vweight_map_t(1);
    if(eweight.empty())
        eweight = eweight_map_t(1);

    typedef property_map_type::apply<uint8_t,
                                     GraphInterface::vertex_index_map_t>::type
        pin_map_t;
    pin_map_t pin_map = any_cast<pin_map_t>(pin);

    run_action<graph_tool::detail::never_directed>()
        (g,
         bind<void>(get_sfdp_layout(C, K, p, theta, gamma, mu, mu_p, init_step,
                                    step_schedule, max_level, epsilon,
                                    max_iter, adaptive),
                    _1, g.GetVertexIndex(), _2, _3, _4,
                    pin_map.get_unchecked(num_vertices(g.GetGraph())),
                    groups.get_unchecked(num_vertices(g.GetGraph())), verbose),
         vertex_floating_vector_properties(), vertex_props_t(), edge_props_t())
        (pos, vweight, eweight);
}

struct do_propagate_pos
{
    template <class Graph, class CoarseGraph, class VertexMap, class PosMap,
              class RNG>
    void operator()(Graph& g, CoarseGraph* cg, VertexMap& vmap,
                    boost::any acvmap, PosMap pos, boost::any acpos,
                    double delta, RNG& rng) const
    {
        typename PosMap::checked_t cpos =
            any_cast<typename PosMap::checked_t>(acpos);
        typename VertexMap::checked_t cvmap =
            any_cast<typename VertexMap::checked_t>(acvmap);
        typedef typename property_traits<VertexMap>::value_type c_t;
        typedef typename property_traits<PosMap>::value_type pos_t;
        typedef typename pos_t::value_type val_t;

        tr1::variate_generator<RNG&, tr1::uniform_real<val_t> >
            noise(rng, tr1::uniform_real<val_t>(-delta, delta));

        tr1::unordered_map<c_t, pos_t, boost::hash<c_t> >
            cmap(num_vertices(*cg));
        int i, N = num_vertices(*cg);
        #pragma omp parallel for default(shared) private(i)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<CoarseGraph>::vertex_descriptor v =
                vertex(i, *cg);
            if (v == graph_traits<CoarseGraph>::null_vertex())
                continue;
            cmap[cvmap[v]] = cpos[v];
        }

        N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v =
                vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            pos[v] = cmap[vmap[v]];
            {
                if (delta > 0)
                {
                    #pragma omp critical
                    for (size_t j = 0; j < pos[v].size(); ++j)
                        pos[v][j] += noise();
                }
            }
        }
    }
};

struct get_pointers
{
    template <class List>
    struct apply
    {
        typedef typename mpl::transform<List,
                                        mpl::quote1<add_pointer> >::type type;
    };
};

void propagate_pos(GraphInterface& gi, GraphInterface& cgi, boost::any vmap,
                   boost::any cvmap, boost::any pos, boost::any cpos,
                   double delta, rng_t& rng)
{
    typedef mpl::vector<property_map_type::apply
                            <int32_t,
                             GraphInterface::vertex_index_map_t>::type>::type
        vmaps_t;


    run_action<>()
        (gi, bind<void>(do_propagate_pos(),
                        _1, _2, _3, cvmap, _4, cpos, delta, ref(rng)),
         get_pointers::apply<graph_tool::detail::all_graph_views>::type(),
         vmaps_t(), vertex_floating_vector_properties())
        (cgi.GetGraphView(), vmap, pos);
}

struct do_propagate_pos_mivs
{
    template <class Graph, class MIVSMap, class PosMap,
              class RNG>
    void operator()(Graph& g, MIVSMap mivs, PosMap pos, double delta, RNG& rng) const
    {
        typedef typename property_traits<PosMap>::value_type pos_t;
        typedef typename pos_t::value_type val_t;

        tr1::variate_generator<RNG&, tr1::uniform_real<val_t> >
            noise(rng, tr1::uniform_real<val_t>(-delta, delta));

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v =
                vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            if (mivs[v])
                continue;
            pos[v].resize(2);
            pos[v][0] = pos[v][1] = 0;
            size_t count = 0;
            typename graph_traits<Graph>::adjacency_iterator a, a_end;
            for(tie(a, a_end) = adjacent_vertices(v, g); a != a_end; ++a)
            {
                if (!mivs[*a])
                    continue;
                for (size_t j = 0; j < pos[v].size(); ++j)
                    pos[v][j] += pos[*a][j];
                ++count;
            }

            if (count == 0)
                throw ValueException("invalid MIVS! Vertex has no neighbours "
                                     "belonging to the set!");

            if (count == 1)
            {
                if (delta > 0)
                {
                    #pragma omp critical
                    for (size_t j = 0; j < pos[v].size(); ++j)
                        pos[v][j] += noise();
                }
            }
            else
            {
                for (size_t j = 0; j < pos[v].size(); ++j)
                    pos[v][j] /= count;
            }
        }
    }
};


void propagate_pos_mivs(GraphInterface& gi, boost::any mivs, boost::any pos,
                        double delta, rng_t& rng)
{
    run_action<>()
        (gi, bind<void>(do_propagate_pos_mivs(),
                        _1, _2, _3, delta, ref(rng)),
         vertex_scalar_properties(), vertex_floating_vector_properties())
        (mivs, pos);
}


struct do_avg_dist
{
    template <class Graph, class PosMap>
    void operator()(Graph& g, PosMap pos, double& ad) const
    {
        int i, N = num_vertices(g);
        size_t count = 0;
        double d = 0;
        #pragma omp parallel for default(shared) private(i) \
            reduction(+: d, count)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v =
                vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            typename graph_traits<Graph>::adjacency_iterator a, a_end;
            for(tie(a, a_end) = adjacent_vertices(v, g); a != a_end; ++a)
            {
                d += dist(pos[v], pos[*a]);
                count++;
            }
        }
        d /= count;
        ad = d;
    }
};


double avg_dist(GraphInterface& gi, boost::any pos)
{
    double d;
    run_action<>()
        (gi, bind<void>(do_avg_dist(), _1, _2, ref(d)),
         vertex_scalar_vector_properties()) (pos);
    return d;
}

#include <boost/python.hpp>

void export_sfdp()
{
    python::def("sfdp_layout", &sfdp_layout);
    python::def("propagate_pos", &propagate_pos);
    python::def("propagate_pos_mivs", &propagate_pos_mivs);
    python::def("avg_dist", &avg_dist);
}
