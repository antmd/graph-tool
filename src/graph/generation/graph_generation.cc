// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007  Tiago de Paula Peixoto <tiago@forked.de>
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

    double operator()(pair<size_t, size_t> deg, pair<size_t, size_t> degl) const
    {
        python::object ret = _o(python::make_tuple(deg.first, deg.second),
                                python::make_tuple(degl.first, degl.second));
        return python::extract<double>(ret);
    }

private:
    python::object _o;
};

void generate_random_graph(GraphInterface& gi, size_t N,
                           python::object deg_sample,
                           python::object corr_deg_sample,
                           bool uncorrelated, bool no_parallel,
                           bool no_self_loops, bool undirected,
                           size_t seed, bool verbose)
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
            (gi, bind<void>(gen_random_graph<mpl::bool_<false> >(N), _1,
                            PythonFuncWrap(deg_sample),
                            PythonFuncWrap(corr_deg_sample),
                            no_parallel, no_self_loops,
                            undirected, seed, verbose))();
    }
    else
    {
        run_action<graph_views>()
            (gi, bind<void>(gen_random_graph<mpl::bool_<true> >(N), _1,
                            PythonFuncWrap(deg_sample),
                            PythonFuncWrap(corr_deg_sample),
                            no_parallel, no_self_loops,
                            undirected, seed, verbose))();
    }
    gi.ReIndexEdges();
}

void random_rewire(GraphInterface& gi, string strat, bool self_loops,
                   bool parallel_edges, size_t seed);
void predecessor_graph(GraphInterface& gi, GraphInterface& gpi,
                       boost::any pred_map);

using namespace boost::python;

BOOST_PYTHON_MODULE(libgraph_tool_generation)
{
    def("gen_random_graph", &generate_random_graph);
    def("random_rewire", &random_rewire);
    def("predecessor_graph", &predecessor_graph);
}
