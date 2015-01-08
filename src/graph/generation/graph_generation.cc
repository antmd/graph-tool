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

#include "graph.hh"
#include "graph_util.hh"
#include "graph_filtering.hh"
#include "graph_generation.hh"
#include "sampler.hh"
#include "dynamic_sampler.hh"
#include <boost/python.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

class PythonFuncWrap
{
public:
    PythonFuncWrap(boost::python::object o): _o(o) {}

    pair<size_t, size_t> operator()(size_t i) const
    {
        boost::python::object ret = _o(i);
        return boost::python::extract<pair<size_t,size_t> >(ret);
    }

    size_t operator()(size_t i, bool) const
    {
        boost::python::object ret = _o(i);
        return boost::python::extract<size_t>(ret);
    }

private:
    boost::python::object _o;
};

void generate_graph(GraphInterface& gi, size_t N, boost::python::object deg_sample,
                    bool no_parallel, bool no_self_loops, bool undirected,
                    rng_t& rng, bool verbose, bool verify)
{
    typedef graph_tool::detail::get_all_graph_views::apply<
    graph_tool::detail::filt_scalar_type, boost::mpl::bool_<false>,
        boost::mpl::bool_<false>, boost::mpl::bool_<false>,
        boost::mpl::bool_<true>, boost::mpl::bool_<true> >::type graph_views;

    if (undirected)
        gi.SetDirected(false);

    run_action<graph_views>()
        (gi, std::bind(gen_graph(), placeholders::_1, N,
                       PythonFuncWrap(deg_sample),
                       no_parallel, no_self_loops,
                       std::ref(rng), verbose, verify))();
}

size_t random_rewire(GraphInterface& gi, string strat, size_t niter,
                     bool no_sweep, bool self_loops, bool parallel_edges,
                     bool alias, bool traditional, bool persist,
                     boost::python::object corr_prob, boost::any apin,
                     boost::any block, bool cache, rng_t& rng, bool verbose);
void predecessor_graph(GraphInterface& gi, GraphInterface& gpi,
                       boost::any pred_map);
void line_graph(GraphInterface& gi, GraphInterface& lgi,
                boost::any edge_index);
boost::python::tuple graph_union(GraphInterface& ugi, GraphInterface& gi,
                          boost::any avprop);
void vertex_property_union(GraphInterface& ugi, GraphInterface& gi,
                           boost::any p_vprop, boost::any p_eprop,
                           boost::any uprop, boost::any prop);
void edge_property_union(GraphInterface& ugi, GraphInterface& gi,
                         boost::any p_vprop, boost::any p_eprop,
                         boost::any uprop, boost::any prop);
void triangulation(GraphInterface& gi, boost::python::object points, boost::any pos,
                   string type, bool periodic);
void lattice(GraphInterface& gi, boost::python::object oshape, bool periodic);
void geometric(GraphInterface& gi, boost::python::object opoints, double r,
               boost::python::object orange, bool periodic, boost::any pos);
void price(GraphInterface& gi, size_t N, double gamma, double c, size_t m,
           rng_t& rng);
void complete(GraphInterface& gi, size_t N, bool directed, bool self_loops);
void circular(GraphInterface& gi, size_t N, size_t k, bool directed, bool self_loops);

using namespace boost::python;

BOOST_PYTHON_MODULE(libgraph_tool_generation)
{
    def("gen_graph", &generate_graph);
    def("random_rewire", &random_rewire);
    def("predecessor_graph", &predecessor_graph);
    def("line_graph", &line_graph);
    def("graph_union", &graph_union);
    def("vertex_property_union", &vertex_property_union);
    def("edge_property_union", &edge_property_union);
    def("triangulation", &triangulation);
    def("lattice", &lattice);
    def("geometric", &geometric);
    def("price", &price);
    def("complete", &complete);
    def("circular", &circular);

    class_<Sampler<int, boost::mpl::false_>>("Sampler",
                                             init<const vector<int>&, const vector<double>&>())
        .def("sample", &Sampler<int, boost::mpl::false_>::sample<rng_t>,
             return_value_policy<copy_const_reference>());

    class_<DynamicSampler<int>>("DynamicSampler",
                                init<const vector<int>&,
                                     const vector<double>&>())
        .def("sample", &DynamicSampler<int>::sample<rng_t>,
             return_value_policy<copy_const_reference>())
        .def("insert", &DynamicSampler<int>::insert)
        .def("remove", &DynamicSampler<int>::remove)
        .def("reset", &DynamicSampler<int>::reset)
        .def("rebuild", &DynamicSampler<int>::rebuild);
}
