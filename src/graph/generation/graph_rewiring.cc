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
    PythonFuncWrap(boost::python::object o): _o(o) {}

    double operator()(pair<size_t, size_t> deg1, pair<size_t, size_t> deg2)
        const
    {
        boost::python::object ret = _o(boost::python::make_tuple(deg1.first, deg1.second),
                                       boost::python::make_tuple(deg2.first, deg2.second));
        return boost::python::extract<double>(ret);
    }

    template <class Type>
    double operator()(const Type& deg1, const Type& deg2) const
    {
        boost::python::object ret = _o(boost::python::object(deg1), boost::python::object(deg2));
        return boost::python::extract<double>(ret);
    }

private:
    boost::python::object _o;
};

struct graph_rewire_block
{
    graph_rewire_block(bool alias, bool traditional) : alias(alias), traditional(traditional) {}
    bool alias;
    bool traditional;

    template <class Graph, class EdgeIndexMap, class CorrProb, class BlockProp>
    void operator()(Graph& g, EdgeIndexMap edge_index, CorrProb corr_prob,
                    pair<bool, bool> rest, BlockProp block_prop,
                    pair<size_t, bool> iter_sweep,
                    std::tuple<bool, bool, bool> cache_verbose, size_t& pcount, rng_t& rng)
        const
    {
        if (traditional)
        {
            graph_rewire<TradBlockRewireStrategy>()
                (g, edge_index, corr_prob, rest.first, rest.second, iter_sweep,
                 cache_verbose, pcount, rng, PropertyBlock<BlockProp>(block_prop));
        }
        else
        {
            if (alias)
                graph_rewire<AliasProbabilisticRewireStrategy>()
                    (g, edge_index, corr_prob, rest.first, rest.second, iter_sweep,
                     cache_verbose, pcount, rng, PropertyBlock<BlockProp>(block_prop));
            else
                graph_rewire<ProbabilisticRewireStrategy>()
                    (g, edge_index, corr_prob, rest.first, rest.second, iter_sweep,
                     cache_verbose, pcount, rng, PropertyBlock<BlockProp>(block_prop));
        }
    }
};


size_t random_rewire(GraphInterface& gi, string strat, size_t niter,
                     bool no_sweep, bool self_loops, bool parallel_edges,
                     bool alias, bool traditional, bool persist,
                     boost::python::object corr_prob, boost::any block, bool cache,
                     rng_t& rng, bool verbose)
{
    PythonFuncWrap corr(corr_prob);
    size_t pcount = 0;

    if (strat == "erdos")
        run_action<graph_tool::detail::never_reversed>()
            (gi, std::bind(graph_rewire<ErdosRewireStrategy>(),
                           placeholders::_1, gi.GetEdgeIndex(),
                           std::ref(corr),
                           self_loops, parallel_edges,
                           make_pair(niter, no_sweep),
                           std::make_tuple(persist, cache, verbose),
                           std::ref(pcount), std::ref(rng)))();
    else if (strat == "uncorrelated")
        run_action<graph_tool::detail::never_reversed>()
            (gi, std::bind(graph_rewire<RandomRewireStrategy>(),
                           placeholders::_1, gi.GetEdgeIndex(), std::ref(corr),
                           self_loops, parallel_edges,
                           make_pair(niter, no_sweep),
                           std::make_tuple(persist, cache, verbose),
                           std::ref(pcount), std::ref(rng)))();
    else if (strat == "correlated")
        run_action<graph_tool::detail::never_reversed>()
            (gi, std::bind(graph_rewire<CorrelatedRewireStrategy>(),
                           placeholders::_1, gi.GetEdgeIndex(), std::ref(corr),
                           self_loops, parallel_edges,
                           make_pair(niter, no_sweep),
                           std::make_tuple(persist, cache, verbose),
                           std::ref(pcount), std::ref(rng)))();
    else if (strat == "probabilistic")
        run_action<>()
            (gi, std::bind(graph_rewire<ProbabilisticRewireStrategy>(),
                           placeholders::_1, gi.GetEdgeIndex(), std::ref(corr),
                           self_loops, parallel_edges,
                           make_pair(niter, no_sweep),
                           std::make_tuple(persist, cache, verbose),
                           std::ref(pcount), std::ref(rng)))();
    else if (strat == "blockmodel")
        run_action<>()
            (gi, std::bind(graph_rewire_block(alias, traditional),
                           placeholders::_1, gi.GetEdgeIndex(),
                           std::ref(corr),
                           make_pair(self_loops, parallel_edges),
                           placeholders::_2,
                           make_pair(niter, no_sweep),
                           std::make_tuple(persist, cache, verbose),
                           std::ref(pcount), std::ref(rng)),
             vertex_properties())(block);
    else
        throw ValueException("invalid random rewire strategy: " + strat);
    return pcount;
}
