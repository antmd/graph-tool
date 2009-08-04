#! /usr/bin/env python
# graph_tool.py -- a general graph manipulation python module
#
# Copyright (C) 2007 Tiago de Paula Peixoto <tiago@forked.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
``graph_tool.generation`` - Random Graph Generation
---------------------------------------------------
"""

from .. dl_import import dl_import
dl_import("import libgraph_tool_generation")

from .. core import Graph
import sys, numpy

__all__ = ["random_graph", "random_rewire"]

def _corr_wrap(i, j, corr):
    return corr(i[1], j[1])

def random_graph(N, deg_sampler, deg_corr=None, directed=True,
                 parallel=False, self_loops=False,
                 seed=0, verbose=False):
    r"""
    Generate a random graph, with a given degree distribution and correlation.

    Parameters
    ----------
    N : int
        Number of vertices in the graph.
    deg_sampler : function
        A degree sampler function which is called without arguments, and returns
        a tuple of ints representing the in and out-degree of a given vertex (or
        a single int for undirected graphs, representing the out-degree). This
        function is called once per vertex, but may be called more times, if the
        degree sequence cannot be used to build a graph.
    deg_corr : function (optional, default: None)
        A function which give the degree correlation of the graph. It should be
        callable with two parameters: the in,out-degree pair of the source
        vertex an edge, and the in,out-degree pair of the target of the same
        edge (for undirected graphs, both parameters are single values). The
        function should return a number proportional to the probability of such
        an edge existing in the generated graph.
    directed : bool (optional, default: True)
        Whether the generated graph should be directed.
    parallel : bool (optional, default: False)
        If True, parallel edges are allowed.
    self_loops : bool (optional, default: False)
        If True, self-loops are allowed.
    seed : int (optional, default: 0)
        Seed for the random number generator. If seed=0, a random value is
        chosen.

    Returns
    -------
    random_graph : Graph
        The generated graph.

    See Also
    --------
    random_rewire: in place graph shuffling

    Notes
    -----
    The algorithm maintains a list of all available source and target degree
    pairs, such that the deg_corr function is called only once with the same
    parameters.

    The uncorrelated case, the complexity is :math:`O(V+E)`. For the correlated
    case the worst-case complexity is :math:`O(V^2)`, but the typical case has
    complexity :math:`O(V + E\log N_k + N_k^2)`, where :math:`N_k < V` is the
    number of different degrees sampled (or in,out-degree pairs).

    Examples
    --------

    >>> from numpy.random import randint, random, seed, poisson
    >>> from pylab import *
    >>> seed(42)

    This is a degree sampler which uses rejection sampling to sample from the
    distribution :math:`P(k)\propto 1/k`, up to a maximum.

    >>> def sample_k(max):
    ...     accept = False
    ...     while not accept:
    ...         k = randint(1,max+1)
    ...         accept = random() < 1.0/k
    ...     return k
    ...

    The following generates a random undirected graph with degree distribution
    :math:`P(k)\propto 1/k` (with k_max=40) and an *assortative* degree
    correlation of the form:

    .. math::

        P(i,k) \propto \frac{1}{1+|i-k|}

    >>> g = gt.random_graph(1000, lambda: sample_k(40),
    ...                     lambda i,k: 1.0/(1+abs(i-k)), directed=False)
    >>> gt.scalar_assortativity(g, "out")
    (0.59472179721535989, 0.011919463022240388)

    The following samples an in,out-degree pair from the joint distribution:

    .. math::

        p(j,k) = \frac{1}{2}\frac{e^{-m_1}m_1^j}{j!}\frac{e^{-m_1}m_1^k}{k!} +
                 \frac{1}{2}\frac{e^{-m_2}m_2^j}{j!}\frac{e^{-m_2}m_2^k}{k!}

    with :math:`m_1 = 4` and :math:`m_2 = 20`.

    >>> def deg_sample():
    ...    if random() > 0.5:
    ...        return poisson(4), poisson(4)
    ...    else:
    ...        return poisson(20), poisson(20)
    ...

    The following generates a random directed graph with this distribution, and
    plots the combined degree correlation.

    >>> g = gt.random_graph(20000, deg_sample)
    >>>
    >>> hist = gt.combined_corr_hist(g, "in", "out")
    >>> imshow(hist[0], interpolation="nearest")
    <...>
    >>> colorbar()
    <...>
    >>> xlabel("in degree")
    <...>
    >>> ylabel("out degree")
    <...>
    >>> savefig("combined-deg-hist.png")

    .. figure:: combined-deg-hist.png
        :align: center

        Combined degree histogram.

    A correlated directed graph can be build as follows. Consider the following
    degree correlation:

    .. math::

         P(j',k'|j,k)=\frac{e^{-k}k^{j'}}{j'!}
         \frac{e^{-(20-j)}(20-j)^{k'}}{k'!}

    i.e., the in->out correlation is "disassortative", the out->in correlation
    is "assortative", and everything else is uncorrelated.
    We will use a flat degree distribution in the range [1,20).

    >>> p = scipy.stats.poisson
    >>> g = gt.random_graph(20000, lambda: (sample_k(19), sample_k(19)),
    ...                                     lambda a,b: (p.pmf(a[0],b[1])*
    ...                                                  p.pmf(a[1],20-b[0])))

    Lets plot the average degree correlations to check.

    >>> clf()
    >>> corr = gt.avg_neighbour_corr(g, "in", "in")
    >>> errorbar(corr[2], corr[0], yerr=corr[1], fmt="o-", label="<in> vs in")
    (...)
    >>> corr = gt.avg_neighbour_corr(g, "in", "out")
    >>> errorbar(corr[2], corr[0], yerr=corr[1], fmt="o-", label="<out> vs in")
    (...)
    >>> corr = gt.avg_neighbour_corr(g, "out", "in")
    >>> errorbar(corr[2], corr[0], yerr=corr[1], fmt="o-", label="<in> vs out")
    (...)
    >>> corr = gt.avg_neighbour_corr(g, "out", "out")
    >>> errorbar(corr[2], corr[0], yerr=corr[1], fmt="o-", label="<out> vs out")
    (...)
    >>> legend(loc="best")
    <...>
    >>> xlabel("source degree")
    <...>
    >>> ylabel("average target degree")
    <...>
    >>> savefig("deg-corr-dir.png")

    .. figure:: deg-corr-dir.png
        :align: center

        Average nearest neighbour correlations.

    """
    if seed == 0:
        seed = numpy.random.randint(0, sys.maxint)
    g = Graph()
    if deg_corr == None:
        uncorrelated = True
    else:
        uncorrelated = False
    if not directed and deg_corr != None:
        corr = lambda i,j: _corr_wrap(i, j, deg_corr)
    else:
        corr = deg_corr
    libgraph_tool_generation.gen_random_graph(g._Graph__graph, N,
                                              deg_sampler, corr,
                                              uncorrelated, not parallel,
                                              not self_loops, not directed,
                                              seed, verbose)
    g.set_directed(directed)
    return g

def random_rewire(g, strat="uncorrelated", self_loops = False,
                  parallel_edges = False, seed = 0):
    if seed != 0:
        seed = random.randint(0, sys.maxint)
    if g.is_reversed():
        was_reversed = True
    else:
        was_reversed = False
    g.set_reversed(False)
    libgraph_tool_generation.random_rewire(g._Graph__graph, strat, self_loops,
                                           parallel_edges, seed)
    if was_reversed:
        g.set_reversed(True)

