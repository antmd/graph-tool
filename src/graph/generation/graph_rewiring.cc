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

#include "graph_python_interface.hh"
#include "graph.hh"
#include "graph_filtering.hh"

#include <boost/bind.hpp>
#include <boost/python.hpp>

#include "graph_rewiring.hh"

using namespace graph_tool;
using namespace boost;


class PythonFuncWrap
{
public:
    PythonFuncWrap(python::object o): _o(o) {}

    double operator()(pair<size_t, size_t> deg1, pair<size_t, size_t> deg2)
        const
    {
        python::object ret = _o(python::make_tuple(deg1.first, deg1.second),
                                python::make_tuple(deg2.first, deg2.second));
        return python::extract<double>(ret);
    }

    template <class Type>
    double operator()(const Type& deg1, const Type& deg2) const
    {
        python::object ret = _o(python::object(deg1), python::object(deg2));
        return python::extract<double>(ret);
    }

private:
    python::object _o;
};


struct graph_rewire_block
{
    template <class Graph, class EdgeIndexMap, class CorrProb, class BlockProp>
    void operator()(Graph& g, EdgeIndexMap edge_index, CorrProb corr_prob,
                    pair<bool, bool> rest, BlockProp block_prop,
                    pair<size_t, bool> iter_sweep,
                    pair<bool, bool> cache_verbose, size_t& pcount, rng_t& rng)
        const
    {
        graph_rewire<ProbabilisticRewireStrategy>()
            (g, edge_index, corr_prob, rest.first, rest.second, iter_sweep,
             cache_verbose, pcount, rng, PropertyBlock<BlockProp>(block_prop));
    }
};


size_t random_rewire(GraphInterface& gi, string strat, size_t niter,
                     bool no_sweep, bool self_loops, bool parallel_edges,
                     python::object corr_prob, boost::any block,
                     bool cache, rng_t& rng, bool verbose)
{
    PythonFuncWrap corr(corr_prob);
    size_t pcount = 0;

    if (strat == "erdos")
        run_action<graph_tool::detail::never_reversed>()
            (gi, boost::bind<void>(graph_rewire<ErdosRewireStrategy>(),
                                   _1, gi.GetEdgeIndex(), boost::ref(corr),
                                   self_loops, parallel_edges,
                                   make_pair(niter, no_sweep),
                                   make_pair(cache, verbose),
                                   boost::ref(pcount), boost::ref(rng)))();
    else if (strat == "uncorrelated")
        run_action<graph_tool::detail::never_reversed>()
            (gi, boost::bind<void>(graph_rewire<RandomRewireStrategy>(),
                                   _1, gi.GetEdgeIndex(), boost::ref(corr),
                                   self_loops, parallel_edges,
                                   make_pair(niter, no_sweep),
                                   make_pair(cache, verbose),
                                   boost::ref(pcount), boost::ref(rng)))();
    else if (strat == "correlated")
        run_action<graph_tool::detail::never_reversed>()
            (gi, boost::bind<void>(graph_rewire<CorrelatedRewireStrategy>(),
                                   _1, gi.GetEdgeIndex(), boost::ref(corr),
                                   self_loops, parallel_edges,
                                   make_pair(niter, no_sweep),
                                   make_pair(cache, verbose),
                                   boost::ref(pcount), boost::ref(rng)))();
    else if (strat == "probabilistic")
        run_action<>()
            (gi, boost::bind<void>(graph_rewire<ProbabilisticRewireStrategy>(),
                                   _1, gi.GetEdgeIndex(), boost::ref(corr),
                                   self_loops, parallel_edges,
                                   make_pair(niter, no_sweep),
                                   make_pair(cache, verbose),
                                   boost::ref(pcount), boost::ref(rng)))();
    else if (strat == "blockmodel")
        run_action<>()
            (gi, boost::bind<void>(graph_rewire_block(),
                                   _1, gi.GetEdgeIndex(), boost::ref(corr),
                                   make_pair(self_loops, parallel_edges), _2,
                                   make_pair(niter, no_sweep),
                                   make_pair(cache, verbose),
                                   boost::ref(pcount), boost::ref(rng)),
             vertex_properties())(block);
    else
        throw ValueException("invalid random rewire strategy: " + strat);
    return pcount;
}
