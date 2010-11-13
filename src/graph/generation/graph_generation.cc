// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007-2010 Tiago de Paula Peixoto <tiago@skewed.de>
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
#include <boost/python.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

typedef tr1::mt19937 rng_t;

class PythonFuncWrap
{
public:
    PythonFuncWrap(python::object o): _o(o) {}

    pair<size_t, size_t> operator()() const
    {
        python::object ret = _o();
        return python::extract<pair<size_t,size_t> >(ret);
    }

    size_t operator()(bool) const
    {
        python::object ret = _o();
        return python::extract<size_t>(ret);
    }

private:
    python::object _o;
};

void generate_random_graph(GraphInterface& gi, size_t N,
                           python::object deg_sample,
                           bool uncorrelated, bool no_parallel,
                           bool no_self_loops, bool undirected,
                           size_t seed, bool verbose, bool verify)
{
    typedef graph_tool::detail::get_all_graph_views::apply<
    graph_tool::detail::scalar_pairs, mpl::bool_<false>,
        mpl::bool_<false>, mpl::bool_<false>,
        mpl::bool_<true>, mpl::bool_<true> >::type graph_views;

    if (undirected)
        gi.SetDirected(false);

    if (uncorrelated)
    {
        run_action<graph_views>()
            (gi, bind<void>(gen_random_graph(), _1, N,
                            PythonFuncWrap(deg_sample),
                            no_parallel, no_self_loops,
                            seed, verbose, verify))();
    }
    else
    {
        run_action<graph_views>()
            (gi, bind<void>(gen_random_graph(), _1, N,
                            PythonFuncWrap(deg_sample),
                            no_parallel, no_self_loops,
                            seed, verbose, verify))();
    }
    gi.ReIndexEdges();
}

void random_rewire(GraphInterface& gi, string strat, bool self_loops,
                   bool parallel_edges, python::object corr_prob, size_t seed,
                   bool verbose);
void predecessor_graph(GraphInterface& gi, GraphInterface& gpi,
                       boost::any pred_map);
void line_graph(GraphInterface& gi, GraphInterface& lgi,
                boost::any edge_index);
python::tuple graph_union(GraphInterface& ugi, GraphInterface& gi);
void vertex_property_union(GraphInterface& ugi, GraphInterface& gi,
                           boost::any p_vprop, boost::any p_eprop,
                           boost::any uprop, boost::any prop);
void edge_property_union(GraphInterface& ugi, GraphInterface& gi,
                         boost::any p_vprop, boost::any p_eprop,
                         boost::any uprop, boost::any prop);
void triangulation(GraphInterface& gi, python::object points, boost::any pos,
                   string type, bool periodic);
void lattice(GraphInterface& gi, python::object oshape, bool periodic);
void geometric(GraphInterface& gi, python::object opoints, double r,
               python::object orange, bool periodic, boost::any pos);
void price(GraphInterface& gi, size_t N, double gamma, double c, size_t m,
           size_t seed);

using namespace boost::python;

BOOST_PYTHON_MODULE(libgraph_tool_generation)
{
    def("gen_random_graph", &generate_random_graph);
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
}
