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

    typedef pair<size_t, size_t> deg_t;

    double operator()(deg_t deg1, deg_t deg2)
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

    template <class ProbMap>
    void get_probs(ProbMap& probs) const
    {
        typedef typename ProbMap::key_type::first_type block_t;
        if (PyObject_HasAttrString(_o.ptr(), "__getitem__"))
        {
            int N = boost::python::len(_o);
            for (int i = 0; i < N; ++i)
            {
                block_t ks = boost::python::extract<block_t>(_o[i][0])();
                block_t kt = boost::python::extract<block_t>(_o[i][1])();
                double p = boost::python::extract<double>(_o[i][2])();
                if (std::isnan(p) || std::isinf(p) || p <= 0)
                    continue;
                probs[make_pair(ks, kt)] += p;
            }
        }
    }

private:
    boost::python::object _o;
};

struct graph_rewire_block
{
    graph_rewire_block(bool alias, bool traditional) : alias(alias), traditional(traditional) {}
    bool alias;
    bool traditional;

    template <class Graph, class EdgeIndexMap, class CorrProb, class PinMap,
              class BlockProp>
    void operator()(Graph& g, EdgeIndexMap edge_index, CorrProb corr_prob,
                    PinMap pin, pair<bool, bool> rest, BlockProp block_prop,
                    pair<size_t, bool> iter_sweep,
                    std::tuple<bool, bool, bool> cache_verbose, size_t& pcount, rng_t& rng)
        const
    {
        if (traditional)
        {
            graph_rewire<TradBlockRewireStrategy>()
                (g, edge_index, corr_prob, pin, rest.first, rest.second, iter_sweep,
                 cache_verbose, pcount, rng, PropertyBlock<BlockProp>(block_prop));
        }
        else
        {
            if (alias)
                graph_rewire<AliasProbabilisticRewireStrategy>()
                    (g, edge_index, corr_prob, pin, rest.first, rest.second, iter_sweep,
                     cache_verbose, pcount, rng, PropertyBlock<BlockProp>(block_prop));
            else
                graph_rewire<ProbabilisticRewireStrategy>()
                    (g, edge_index, corr_prob, pin,  rest.first, rest.second, iter_sweep,
                     cache_verbose, pcount, rng, PropertyBlock<BlockProp>(block_prop));
        }
    }
};


struct graph_rewire_correlated
{
    template <class Graph, class EdgeIndexMap, class CorrProb, class PinMap,
              class BlockProp>
    void operator()(Graph& g, EdgeIndexMap edge_index, CorrProb corr_prob,
                    PinMap pin, bool self_loops, bool parallel_edges,
                    pair<size_t, bool> iter_sweep,
                    std::tuple<bool, bool, bool> cache_verbose,
                    size_t& pcount, rng_t& rng, BlockProp block_prop) const
    {
        graph_rewire<CorrelatedRewireStrategy>()
            (g, edge_index, corr_prob, pin, self_loops, parallel_edges, iter_sweep,
             cache_verbose, pcount, rng, PropertyBlock<BlockProp>(block_prop));
    }
};


size_t random_rewire(GraphInterface& gi, string strat, size_t niter,
                     bool no_sweep, bool self_loops, bool parallel_edges,
                     bool alias, bool traditional, bool persist,
                     boost::python::object corr_prob, boost::any apin,
                     boost::any block, bool cache, rng_t& rng, bool verbose)
{
    PythonFuncWrap corr(corr_prob);
    size_t pcount = 0;

    typedef property_map_type::apply<uint8_t,
                                     GraphInterface::edge_index_map_t>::type
        emap_t;

    emap_t::unchecked_t pin =
        any_cast<emap_t>(apin).get_unchecked(gi.GetMaxEdgeIndex());


    if (strat == "erdos")
    {
        run_action<graph_tool::detail::never_reversed>()
            (gi, std::bind(graph_rewire<ErdosRewireStrategy>(),
                           placeholders::_1, gi.GetEdgeIndex(),
                           std::ref(corr), pin,
                           self_loops, parallel_edges,
                           make_pair(niter, no_sweep),
                           std::make_tuple(persist, cache, verbose),
                           std::ref(pcount), std::ref(rng)))();
    }
    else if (strat == "uncorrelated")
    {
        run_action<graph_tool::detail::never_reversed>()
            (gi, std::bind(graph_rewire<RandomRewireStrategy>(),
                           placeholders::_1, gi.GetEdgeIndex(), std::ref(corr),
                           pin, self_loops, parallel_edges,
                           make_pair(niter, no_sweep),
                           std::make_tuple(persist, cache, verbose),
                           std::ref(pcount), std::ref(rng)))();
    }
    else if (strat == "correlated")
    {
        if (block.empty())
        {
            run_action<graph_tool::detail::never_reversed>()
                (gi, std::bind(graph_rewire<CorrelatedRewireStrategy>(),
                               placeholders::_1, gi.GetEdgeIndex(), std::ref(corr),
                               pin, self_loops, parallel_edges,
                               make_pair(niter, no_sweep),
                               std::make_tuple(persist, cache, verbose),
                               std::ref(pcount), std::ref(rng)))();
        }
        else
        {
            run_action<graph_tool::detail::never_reversed>()
                (gi, std::bind(graph_rewire_correlated(),
                               placeholders::_1, gi.GetEdgeIndex(), std::ref(corr),
                               pin, self_loops, parallel_edges,
                               make_pair(niter, no_sweep),
                               std::make_tuple(persist, cache, verbose),
                               std::ref(pcount), std::ref(rng),
                               placeholders::_2),
                 vertex_properties())(block);
        }
    }
    else if (strat == "probabilistic")
    {
        run_action<>()
            (gi, std::bind(graph_rewire<ProbabilisticRewireStrategy>(),
                           placeholders::_1, gi.GetEdgeIndex(), std::ref(corr),
                           pin, self_loops, parallel_edges,
                           make_pair(niter, no_sweep),
                           std::make_tuple(persist, cache, verbose),
                           std::ref(pcount), std::ref(rng)))();
    }
    else if (strat == "blockmodel")
    {
        run_action<>()
            (gi, std::bind(graph_rewire_block(alias, traditional),
                           placeholders::_1, gi.GetEdgeIndex(),
                           std::ref(corr), pin,
                           make_pair(self_loops, parallel_edges),
                           placeholders::_2,
                           make_pair(niter, no_sweep),
                           std::make_tuple(persist, cache, verbose),
                           std::ref(pcount), std::ref(rng)),
             vertex_properties())(block);
    }
    else
    {
        throw ValueException("invalid random rewire strategy: " + strat);
    }
    return pcount;
}
