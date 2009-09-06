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

from .. core import Graph, _check_prop_scalar, _prop
import sys, numpy

__all__ = ["random_graph", "random_rewire", "predecessor_tree", "line_graph",
           "graph_union"]

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

def random_rewire(g, strat= "uncorrelated", parallel_edges = False,
                  self_loops = False, seed = 0):
    r"""
    Shuffled the graph in-place. The degrees (either in or out) of each vertex
    are always the same, but otherwise the edges are randomly placed. If
    strat == "correlated", the degree correlations are also maintained: The new
    source and target of each edge both have the same in and out-degree.

    Parameters
    ----------
    g : Graph
        Graph to be shuffled. The graph will be modified.
    strat : string (optional, default: "uncorrelated")
        If strat == "uncorrelated" only the degrees of the vertices will be
        maintained, nothing else. If strat == "correlated", additionally the new
        source and target of each edge both have the same in and out-degree.
    parallel : bool (optional, default: False)
        If True, parallel edges are allowed.
    self_loops : bool (optional, default: False)
        If True, self-loops are allowed.
    seed : int (optional, default: 0)
        Seed for the random number generator. If seed == 0, a random value is
        chosen.

    Returns
    -------
    None

    See Also
    --------
    random_graph: random graph generation

    Notes
    -----

    Each edge gets swapped at least once, so the overall complexity is
    :math:`O(E)`.

    Examples
    --------

    Some small graphs for visualization.

    >>> from numpy.random import zipf, seed
    >>> from pylab import *
    >>> seed(42)
    >>> g = gt.random_graph(1000, lambda: sample_k(10),
    ...                     lambda i,j: exp(abs(i-j)), directed=False)
    >>> gt.graph_draw(g, layout="arf", output="rewire_orig.png")
    <...>
    >>> gt.random_rewire(g, "correlated")
    >>> gt.graph_draw(g, layout="arf", output="rewire_corr.png")
    <...>
    >>> gt.random_rewire(g)
    >>> gt.graph_draw(g, layout="arf", output="rewire_uncorr.png")
    <...>


    .. figure:: rewire_orig.png
        :align: center

        Original graph. (It is a `ridiculogram <http://www.youtube.com/watch?v=YS-asmU3p_4>`_).

    .. figure:: rewire_corr.png
        :align: center

        Shuffled graph, with degree correlations.

    .. figure:: rewire_uncorr.png
        :align: center

        Shuffled graph, without degree correlations.

    We can try some larger graphs to get better statistics.

    >>> clf()
    >>> g = gt.random_graph(20000, lambda: sample_k(20),
    ...                     lambda i,j: exp(abs(i-j)), directed=False)
    >>> corr = gt.avg_neighbour_corr(g, "out", "out")
    >>> errorbar(corr[2], corr[0], yerr=corr[1], fmt="o", label="original")
    (...)
    >>> gt.random_rewire(g, "correlated")
    >>> corr = gt.avg_neighbour_corr(g, "out", "out")
    >>> errorbar(corr[2], corr[0], yerr=corr[1], fmt="o", label="correlated")
    (...)
    >>> gt.random_rewire(g)
    >>> corr = gt.avg_neighbour_corr(g, "out", "out")
    >>> errorbar(corr[2], corr[0], yerr=corr[1], fmt="o", label="uncorrelated")
    (...)
    >>> xlabel("$k$")
    <...>
    >>> ylabel(r"$\left<k_{nn}\right>$")
    <...>
    >>> legend(loc="best")
    <...>
    >>> savefig("shuffled-stats.png")

    .. figure:: shuffled-stats.png
        :align: center

        Average degree correlations for the different shuffled and non-shuffled
        graphs. The shuffled graph with correlations displays exactly the same
        correlation as the original graph.

    Now let's do it for a directed graph. See
    :func:`~graph_tool.generation.random_graph` for more details.

    >>> p = scipy.stats.poisson
    >>> g = gt.random_graph(20000, lambda: (sample_k(19), sample_k(19)),
    ...                                     lambda a,b: (p.pmf(a[0],b[1])*
    ...                                                  p.pmf(a[1],20-b[0])))
    >>> clf()
    >>> corr = gt.avg_neighbour_corr(g, "in", "out")
    >>> errorbar(corr[2], corr[0], yerr=corr[1], fmt="o-", label="<out> vs in")
    (...)
    >>> corr = gt.avg_neighbour_corr(g, "out", "in")
    >>> errorbar(corr[2], corr[0], yerr=corr[1], fmt="o-", label="<in> vs out")
    (...)
    >>> gt.random_rewire(g, "correlated")
    >>> corr = gt.avg_neighbour_corr(g, "in", "out")
    >>> errorbar(corr[2], corr[0], yerr=corr[1], fmt="o-",
    ...          label="<out> vs in, correlated")
    (...)
    >>> corr = gt.avg_neighbour_corr(g, "out", "in")
    >>> errorbar(corr[2], corr[0], yerr=corr[1], fmt="o-",
    ...          label="<in> vs out, correlated")
    (...)
    >>> gt.random_rewire(g, "uncorrelated")
    >>> corr = gt.avg_neighbour_corr(g, "in", "out")
    >>> errorbar(corr[2], corr[0], yerr=corr[1], fmt="o-",
    ...          label="<out> vs in, uncorrelated")
    (...)
    >>> corr = gt.avg_neighbour_corr(g, "out", "in")
    >>> errorbar(corr[2], corr[0], yerr=corr[1], fmt="o-",
    ...          label="<in> vs out, uncorrelated")
    (...)
    >>> legend(loc="best")
    <...>
    >>> xlabel("source degree")
    <...>
    >>> ylabel("average target degree")
    <...>
    >>> savefig("shuffled-deg-corr-dir.png")

    .. figure:: shuffled-deg-corr-dir.png
        :align: center

        Average degree correlations for the different shuffled and non-shuffled
        directed graphs. The shuffled graph with correlations displays exactly
        the same correlation as the original graph.
    """

    if seed != 0:
        seed = random.randint(0, sys.maxint)

    g.stash_filter(reversed=True)
    libgraph_tool_generation.random_rewire(g._Graph__graph, strat, self_loops,
                                           parallel_edges, seed)
    g.pop_filter(reversed=True)

def predecessor_tree(g, pred_map):
    """Return a graph from a list of predecessors given by
    the 'pred_map' vertex property."""

    _check_prop_scalar(pred_map, "pred_map")
    pg = Graph()
    libgraph_tool_generation.predecessor_graph(g._Graph__graph,
                                               pg._Graph__graph,
                                               _prop("v", g, pred_map))
    return pg

def line_graph(g):
    lg = Graph(directed=g.is_directed())

    vertex_map = lg.new_vertex_property("int64_t")

    libgraph_tool_generation.line_graph(g._Graph__graph,
                                        lg._Graph__graph,
                                        _prop("v", lg, vertex_map))
    return lg, vertex_map

def graph_union(g1, g2, props=[], include=False):
    if not include:
        g1 = Graph(g1)
    g1.stash_filter(directed=True)
    g1.set_directed(True)
    g2.stash_filter(directed=True)
    g2.set_directed(True)
    n_props = []

    try:
        vmap, emap = libgraph_tool_generation.graph_union(g1._Graph__graph,
                                                          g2._Graph__graph)
        for p in props:
            p1, p2 = p
            if not include:
                p1 = g1.copy_property(p1)
            if p2.value_type() != p1.value_type():
                p2 = g2.copy_property(p2, value_type=p1.value_type())
            if p1.key_type() == 'v':
                libgraph_tool_generation.\
                      vertex_property_union(g1._Graph__graph, g2._Graph__graph,
                                            vmap, emap,
                                            _prop(p1.key_type(), g1, p1),
                                            _prop(p2.key_type(), g2, p2))
            else:
                libgraph_tool_generation.\
                      edge_property_union(g1._Graph__graph, g2._Graph__graph,
                                          vmap, emap,
                                          _prop(p1.key_type(), g1, p1),
                                          _prop(p2.key_type(), g2, p2))
            n_props.append(p1)
    finally:
        g1.pop_filter(directed=True)
        g2.pop_filter(directed=True)

    if len(n_props) > 0:
        return g1, n_props
    else:
        return g1
