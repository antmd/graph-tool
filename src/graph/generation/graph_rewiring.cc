// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007-2010 Tiago de Paula Peixoto <tiago@forked.de>
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
#include "graph_filtering.hh"

#if (GCC_VERSION >= 40400)
#   include <tr1/random>
#else
#   include <boost/tr1/random.hpp>
#endif

#include "graph_rewiring.hh"

#include <boost/bind.hpp>
#include <boost/python.hpp>

using namespace graph_tool;
using namespace boost;

class PythonFuncWrap
{
public:
    PythonFuncWrap(python::object o): _o(o) {}

    double operator()(pair<size_t, size_t> deg, pair<size_t, size_t> degl) const
    {
        python::object ret = _o(python::make_tuple(deg.first, deg.second),
                                python::make_tuple(degl.first, degl.second));
        return python::extract<double>(ret);
    }

private:
    python::object _o;
};

void random_rewire(GraphInterface& gi, string strat, bool self_loops,
                   bool parallel_edges, python::object corr_prob, size_t seed,
                   bool verbose)
{
    rng_t rng(static_cast<rng_t::result_type>(seed));
    PythonFuncWrap corr(corr_prob);

    if (strat == "erdos")
        run_action<graph_tool::detail::never_reversed>()
            (gi, boost::bind<void>(graph_rewire<ErdosRewireStrategy>(),
                                   _1, gi.GetEdgeIndex(), boost::ref(corr),
                                   boost::ref(rng), self_loops, parallel_edges,
                                   verbose))();
    else if (strat == "uncorrelated")
        run_action<graph_tool::detail::never_reversed>()
            (gi, boost::bind<void>(graph_rewire<RandomRewireStrategy>(),
                                   _1, gi.GetEdgeIndex(), boost::ref(corr),
                                   boost::ref(rng), self_loops, parallel_edges,
                                   verbose))();
    else if (strat == "correlated")
        run_action<graph_tool::detail::never_reversed>()
            (gi, boost::bind<void>(graph_rewire<CorrelatedRewireStrategy>(),
                                   _1, gi.GetEdgeIndex(), boost::ref(corr),
                                   boost::ref(rng), self_loops, parallel_edges,
                                   verbose))();
    else if (strat == "probabilistic")
        run_action<>()
            (gi, boost::bind<void>(graph_rewire<ProbabilisticRewireStrategy>(),
                                   _1, gi.GetEdgeIndex(), boost::ref(corr),
                                   boost::ref(rng), self_loops, parallel_edges,
                                   verbose))();
    else
        throw ValueException("invalid random rewire strategy: " + strat);
}
