#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# graph_tool -- a general graph manipulation python module
#
# Copyright (C) 2006-2013 Tiago de Paula Peixoto <tiago@skewed.de>
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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.s

"""
``graph_tool.generation`` - Random graph generation
---------------------------------------------------

Summary
+++++++

.. autosummary::
   :nosignatures:

   random_graph
   random_rewire
   predecessor_tree
   line_graph
   graph_union
   triangulation
   lattice
   geometric_graph
   price_network

Contents
++++++++
"""

from __future__ import division, absolute_import, print_function

from .. dl_import import dl_import
dl_import("from . import libgraph_tool_generation")

from .. import Graph, GraphView, _check_prop_scalar, _prop, _limit_args, _gt_type, _get_rng
from .. stats import label_parallel_edges, label_self_loops
import inspect
import types
import sys, numpy, numpy.random

__all__ = ["random_graph", "random_rewire", "predecessor_tree", "line_graph",
           "graph_union", "triangulation", "lattice", "geometric_graph",
           "price_network"]


def random_graph(N, deg_sampler, deg_corr=None, cache_probs=True, directed=True,
                 parallel_edges=False, self_loops=False, blockmodel=None,
                 block_type="int", degree_block=False,
                 random=True, mix_time=10, verbose=False):
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

        Optionally, you can also pass a function which receives one or two
        arguments: If ``blockmodel == None``, the single argument passed will
        be the index of the vertex which will receive the degree.
        If ``blockmodel != None``, the first value passed will be the vertex
        index, and the second will be the block value of the vertex.
    deg_corr : function (optional, default: ``None``)
        A function which gives the degree correlation of the graph. It should be
        callable with two parameters: the in,out-degree pair of the source
        vertex an edge, and the in,out-degree pair of the target of the same
        edge (for undirected graphs, both parameters are single values). The
        function should return a number proportional to the probability of such
        an edge existing in the generated graph.

        If ``blockmodel != None``, the value passed to the function will be the
        block value of the respective vertices, not the in/out-degree pairs.
    cache_probs : bool (optional, default: ``True``)
        If ``True``, the probabilities returned by the ``deg_corr`` parameter
        will be cached internally. This is crucial for good performance, since
        in this case the supplied python function is called only a few times,
        and not at every attempted edge rewire move. However, in the case were
        the different parameter combinations to the probability function is very
        large, the memory requirements to keep the cache may be very large.
    directed : bool (optional, default: ``True``)
        Whether the generated graph should be directed.
    parallel_edges : bool (optional, default: ``False``)
        If ``True``, parallel edges are allowed.
    self_loops : bool (optional, default: ``False``)
        If ``True``, self-loops are allowed.
    blockmodel : list or :class:`~numpy.ndarray` or function (optional, default: ``None``)
        If supplied, the graph will be sampled from a blockmodel ensemble. If
        the value is a list or a :class:`~numpy.ndarray`, it must have
        ``len(block_model) == N``, and the values will define to which block
        each vertex belongs.

        If this value is a function, it will be used to sample the block
        types. It must be callable either with no arguments or with a single
        argument which will be the vertex index. In either case it must return
        a type compatible with the ``block_type`` parameter.
    block_type : string (optional, default: ``"int"``)
        Value type of block labels. Valid only if ``blockmodel != None``.
    degree_block : bool (optional, default: ``False``)
        If ``True``, the degree of each vertex will be appended to block labels
        when constructing the blockmodel, such that the resulting block type
        will be a pair :math:`(r, k)`, where :math:`r` is the original block
        label.
    random : bool (optional, default: ``True``)
        If ``True``, the returned graph is randomized. Otherwise a deterministic
        placement of the edges will be used.
    mix_time : int (optional, default: ``10``)
        Number of edge sweeps to perform in order to mix the graph. This value
        is ignored if ``parallel_edges == self_loops == True`` and
        ``strat != "probabilistic"``.
    verbose : bool (optional, default: ``False``)
        If ``True``, verbose information is displayed.

    Returns
    -------
    random_graph : :class:`~graph_tool.Graph`
        The generated graph.
    blocks : :class:`~graph_tool.PropertyMap`
        A vertex property map with the block values. This is only returned if
        ``blockmodel != None``.

    See Also
    --------
    random_rewire: in place graph shuffling

    Notes
    -----
    The algorithm makes sure the degree sequence is graphical (i.e. realizable)
    and keeps re-sampling the degrees if is not. With a valid degree sequence,
    the edges are placed deterministically, and later the graph is shuffled with
    the :func:`~graph_tool.generation.random_rewire` function, with the
    ``mix_time`` parameter passed as ``n_iter``.

    The complexity is :math:`O(V + E)` if parallel edges are allowed, and
    :math:`O(V + E \times\text{mix-time})` if parallel edges are not allowed.


    .. note ::

        If ``parallel_edges == False`` this algorithm only guarantees that the
        returned graph will be a random sample from the desired ensemble if
        ``mix_time`` is sufficiently large. The algorithm implements an
        efficient Markov chain based on edge swaps, with a mixing time which
        depends on the degree distribution and correlations desired. If degree
        correlations are provided, the mixing time tends to be larger.

        If ``strat == "probabilistic"``, the Markov chain still needs to be
        mixed, even if parallel edges and self-loops are allowed. In this case
        the Markov chain is implemented using the Metropolis-Hastings
        [metropolis-equations-1953]_ [hastings-monte-carlo-1970]_
        acceptance/rejection algorithm.

    Examples
    --------

    .. testcode::
       :hide:

       from numpy.random import randint, random, seed, poisson
       from pylab import *
       seed(43)
       gt.seed_rng(42)

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
    ...                     lambda i, k: 1.0 / (1 + abs(i - k)), directed=False,
    ...                     mix_time=100)
    >>> gt.scalar_assortativity(g, "out")
    (0.606974415793802, 0.010937932958096811)

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
    >>>
    >>> clf()
    >>> imshow(hist[0].T, interpolation="nearest", origin="lower")
    <...>
    >>> colorbar()
    <...>
    >>> xlabel("in-degree")
    <...>
    >>> ylabel("out-degree")
    <...>
    >>> savefig("combined-deg-hist.pdf")

    .. testcode::
       :hide:

       savefig("combined-deg-hist.png")

    .. figure:: combined-deg-hist.*
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
    ...                     lambda a,b: (p.pmf(a[0], b[1]) *
    ...                                  p.pmf(a[1], 20 - b[0])),
    ...                     mix_time=100)

    Lets plot the average degree correlations to check.

    >>> clf()
    >>> axes([0.1,0.15,0.63,0.8])
    <...>
    >>> corr = gt.avg_neighbour_corr(g, "in", "in")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="o-",
    ...         label=r"$\left<\text{in}\right>$ vs in")
    <...>
    >>> corr = gt.avg_neighbour_corr(g, "in", "out")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="o-",
    ...         label=r"$\left<\text{out}\right>$ vs in")
    <...>
    >>> corr = gt.avg_neighbour_corr(g, "out", "in")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="o-",
    ...          label=r"$\left<\text{in}\right>$ vs out")
    <...>
    >>> corr = gt.avg_neighbour_corr(g, "out", "out")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="o-",
    ...          label=r"$\left<\text{out}\right>$ vs out")
    <...>
    >>> legend(bbox_to_anchor=(1.01, 0.5), loc="center left", borderaxespad=0.)
    <...>
    >>> xlabel("Source degree")
    <...>
    >>> ylabel("Average target degree")
    <...>
    >>> savefig("deg-corr-dir.pdf")

    .. testcode::
       :hide:

       savefig("deg-corr-dir.png")

    .. figure:: deg-corr-dir.*
        :align: center

        Average nearest neighbour correlations.


    **Blockmodels**


    The following example shows how a stochastic blockmodel
    [holland-stochastic-1983]_ [karrer-stochastic-2011]_ can be generated. We
    will consider a system of 10 blocks, which form communities. The connection
    probability will be given by

    >>> def corr(a, b):
    ...    if a == b:
    ...        return 0.999
    ...    else:
    ...        return 0.001

    The blockmodel can be generated as follows.

    >>> g, bm = gt.random_graph(5000, lambda: poisson(10), directed=False,
    ...                         blockmodel=lambda: randint(10), deg_corr=corr,
    ...                         mix_time=500)
    >>> gt.graph_draw(g, vertex_fill_color=bm, output="blockmodel.pdf")
    <...>

    .. testcode::
       :hide:

       gt.graph_draw(g, vertex_fill_color=bm, output="blockmodel.png")

    .. figure:: blockmodel.*
        :align: center

        Simple blockmodel with 10 blocks.


    References
    ----------
    .. [metropolis-equations-1953]  Metropolis, N.; Rosenbluth, A.W.;
       Rosenbluth, M.N.; Teller, A.H.; Teller, E. "Equations of State
       Calculations by Fast Computing Machines". Journal of Chemical Physics 21
       (6): 1087-1092 (1953). :doi:`10.1063/1.1699114`
    .. [hastings-monte-carlo-1970] Hastings, W.K. "Monte Carlo Sampling Methods
       Using Markov Chains and Their Applications". Biometrika 57 (1): 97-109 (1970).
       :doi:`10.1093/biomet/57.1.97`
    .. [holland-stochastic-1983] Paul W. Holland, Kathryn Blackmond Laskey, and
       Samuel Leinhardt, "Stochastic blockmodels: First steps," Social Networks
       5, no. 2: 109-13 (1983) :doi:`10.1016/0378-8733(83)90021-7`
    .. [karrer-stochastic-2011] Brian Karrer and M. E. J. Newman, "Stochastic
       blockmodels and community structure in networks," Physical Review E 83,
       no. 1: 016107 (2011) :doi:`10.1103/PhysRevE.83.016107` :arxiv:`1008.3926`
    """

    g = Graph()
    if deg_corr == None:
        uncorrelated = True
    else:
        uncorrelated = False

    if (type(blockmodel) is types.FunctionType or
        type(blockmodel) is types.LambdaType):
        btype = block_type
        bm = []
        if len(inspect.getargspec(blockmodel)[0]) == 0:
            for i in range(N):
                bm.append(blockmodel())
        else:
            for i in range(N):
                bm.append(blockmodel(i))
        blockmodel = bm
    elif blockmodel is not None:
        btype = _gt_type(blockmodel[0])

    if len(inspect.getargspec(deg_sampler)[0]) > 0:
        if blockmodel is not None:
            sampler = lambda i: deg_sampler(i, blockmodel[i])
        else:
            sampler = deg_sampler
    else:
        sampler = lambda i: deg_sampler()

    libgraph_tool_generation.gen_graph(g._Graph__graph, N, sampler,
                                       uncorrelated, not parallel_edges,
                                       not self_loops, not directed,
                                       _get_rng(), verbose, True)
    g.set_directed(directed)

    if degree_block:
        if btype in ["object", "string"] or "vector" in btype:
            btype = "object"
        elif btype in ["int", "int32_t", "bool"]:
            btype = "vector<int32_t>"
        elif btype in ["long", "int64_t"]:
            btype = "vector<int64_t>"
        elif btype in ["double"]:
            btype = "vector<double>"
        elif btype in ["long double"]:
            btype = "vector<long double>"

    if blockmodel is not None:
        bm = g.new_vertex_property(btype)
        if btype in ["object", "string"] or "vector" in btype:
            for v in g.vertices():
                if not degree_block:
                    bm[v] = blockmodel[int(v)]
                else:
                    if g.is_directed():
                        bm[v] = (blockmodel[int(v)], v.in_degree(),
                                 v.out_degree())
                    else:
                        bm[v] = (blockmodel[int(v)], v.out_degree())
        else:
            try:
                bm.a = blockmodel
            except ValueError:
                bm = g.new_vertex_property("object")
                for v in g.vertices():
                    bm[v] = blockmodel[int(v)]
    else:
        bm = None

    if parallel_edges and self_loops and deg_corr is None:
        mix_time = 1
    if random:
        if deg_corr is not None:
            random_rewire(g, strat="probabilistic", n_iter=mix_time,
                          parallel_edges=parallel_edges, deg_corr=deg_corr,
                          cache_probs=cache_probs, self_loops=self_loops,
                          blockmodel=bm, verbose=verbose)
        else:
            random_rewire(g, parallel_edges=parallel_edges, n_iter=mix_time,
                          self_loops=self_loops, verbose=verbose)

    if bm is None:
        return g
    else:
        return g, bm


@_limit_args({"strat": ["erdos", "correlated", "uncorrelated",
                        "probabilistic"]})
def random_rewire(g, strat="uncorrelated", n_iter=1, edge_sweep=True,
                  parallel_edges=False, self_loops=False, deg_corr=None,
                  cache_probs=True, blockmodel=None, ret_fail=False,
                  verbose=False):
    r"""
    Shuffle the graph in-place.

    If ``strat != "erdos"``, the degrees (either in or out) of each vertex are
    always the same, but otherwise the edges are randomly placed. If
    ``strat == "correlated"``, the degree correlations are also maintained: The
    new source and target of each edge both have the same in and out-degree. If
    ``strat == "probabilistic"``, then edges are rewired according to the degree
    correlation given by the parameter ``deg_corr``.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be shuffled. The graph will be modified.
    strat : string (optional, default: ``"uncorrelated"``)
        If ``strat == "erdos"``, the resulting graph will be entirely random. If
        ``strat == "uncorrelated"`` only the degrees of the vertices will be
        maintained, nothing else. If ``strat == "correlated"``, additionally the
        new source and target of each edge both have the same in and out-degree.
        If ``strat == "probabilistic"``, than edges are rewired according to the
        degree correlation given by the parameter ``deg_corr``.
    n_iter : int (optional, default: ``1``)
        Number of iterations. If ``edge_sweep == True``, each iteration
        corresponds to an entire "sweep" over all edges. Otherwise this
        corresponds to the total number of edges which are randomly chosen for a
        swap attempt (which may repeat).
    edge_sweep : bool (optional, default: ``True``)
        If ``True``, each iteration will perform an entire "sweep" over the
        edges, where each edge is visited once in random order, and a edge swap
        is attempted.
    parallel : bool (optional, default: ``False``)
        If ``True``, parallel edges are allowed.
    self_loops : bool (optional, default: ``False``)
        If ``True``, self-loops are allowed.
    deg_corr : function (optional, default: ``None``)
        A function which gives the degree correlation of the graph. It should be
        callable with two parameters: the in,out-degree pair of the source
        vertex an edge, and the in,out-degree pair of the target of the same
        edge (for undirected graphs, both parameters are single values). The
        function should return a number proportional to the probability of such
        an edge existing in the generated graph. This parameter is ignored,
        unless ``strat == "probabilistic"``.

        If ``blockmodel != None``, the value passed to the function will be the
        block value of the respective vertices, not the in/out-degree pairs.
    cache_probs : bool (optional, default: ``True``)
        If ``True``, the probabilities returned by the ``deg_corr`` parameter
        will be cached internally. This is crucial for good performance, since
        in this case the supplied python function is called only a few times,
        and not at every attempted edge rewire move. However, in the case were
        the different parameter combinations to the probability function is very
        large, the memory requirements to keep the cache may be very large.
    blockmodel : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        If supplied, the graph will be rewired to conform to a blockmodel
        ensemble. The value must be a vertex property map which defines the
        block of each vertex.
    ret_fail : bool (optional, default: ``False``)
        If ``True``, the number of failed edge moves (due to parallel edges or
        self-loops) is returned.
    verbose : bool (optional, default: ``False``)
        If ``True``, verbose information is displayed.


    Returns
    -------
    fail_count : int
        Number of failed edge moves (due to parallel edges or self-loops).

    See Also
    --------
    random_graph: random graph generation

    Notes
    -----
    This algorithm iterates through all the edges in the network and tries to
    swap its target or source with the target or source of another edge.

    .. note::

        If ``parallel_edges = False``, parallel edges are not placed during
        rewiring. In this case, the returned graph will be a uncorrelated sample
        from the desired ensemble only if ``n_iter`` is sufficiently large. The
        algorithm implements an efficient Markov chain based on edge swaps, with
        a mixing time which depends on the degree distribution and correlations
        desired. If degree probabilistic correlations are provided, the mixing
        time tends to be larger.

        If ``strat == "probabilistic"``, the Markov chain still needs to be
        mixed, even if parallel edges and self-loops are allowed. In this case
        the Markov chain is implemented using the Metropolis-Hastings
        [metropolis-equations-1953]_ [hastings-monte-carlo-1970]_
        acceptance/rejection algorithm.


    Each edge is tentatively swapped once per iteration, so the overall
    complexity is :math:`O(V + E \times \text{n-iter})`. If ``edge_sweep ==
    False``, the complexity becomes :math:`O(V + E + \text{n-iter})`.

    Examples
    --------

    Some small graphs for visualization.

    .. testcode::
       :hide:

       from numpy.random import random, seed
       from pylab import *
       seed(43)
       gt.seed_rng(42)

    >>> g, pos = gt.triangulation(random((1000,2)))
    >>> pos = gt.arf_layout(g)
    >>> gt.graph_draw(g, pos=pos, output="rewire_orig.pdf", output_size=(300, 300))
    <...>

    .. testcode::
       :hide:

       gt.graph_draw(g, pos=pos, output="rewire_orig.png", output_size=(300, 300))

    >>> gt.random_rewire(g, "correlated")
    >>> pos = gt.arf_layout(g)
    >>> gt.graph_draw(g, pos=pos, output="rewire_corr.pdf", output_size=(300, 300))
    <...>

    .. testcode::
       :hide:

       gt.graph_draw(g, pos=pos, output="rewire_corr.png", output_size=(300, 300))

    >>> gt.random_rewire(g)
    >>> pos = gt.arf_layout(g)
    >>> gt.graph_draw(g, pos=pos, output="rewire_uncorr.pdf", output_size=(300, 300))
    <...>

    .. testcode::
       :hide:

       gt.graph_draw(g, pos=pos, output="rewire_uncorr.png", output_size=(300, 300))

    >>> gt.random_rewire(g, "erdos")
    >>> pos = gt.arf_layout(g)
    >>> gt.graph_draw(g, pos=pos, output="rewire_erdos.pdf", output_size=(300, 300))
    <...>

    .. testcode::
       :hide:

       gt.graph_draw(g, pos=pos, output="rewire_erdos.png", output_size=(300, 300))

    Some `ridiculograms <http://www.youtube.com/watch?v=YS-asmU3p_4>`_ :

    .. image:: rewire_orig.*
    .. image:: rewire_corr.*
    .. image:: rewire_uncorr.*
    .. image:: rewire_erdos.*

    **From left to right**: Original graph; Shuffled graph, with degree correlations;
    Shuffled graph, without degree correlations; Shuffled graph, with random degrees.

    We can try with larger graphs to get better statistics, as follows.

    >>> figure()
    <...>
    >>> g = gt.random_graph(30000, lambda: sample_k(20),
    ...                     lambda i, j: exp(abs(i-j)), directed=False,
    ...                     mix_time=100)
    >>> corr = gt.avg_neighbour_corr(g, "out", "out")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="o-", label="Original")
    <...>
    >>> gt.random_rewire(g, "correlated")
    >>> corr = gt.avg_neighbour_corr(g, "out", "out")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="*", label="Correlated")
    <...>
    >>> gt.random_rewire(g)
    >>> corr = gt.avg_neighbour_corr(g, "out", "out")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="o-", label="Uncorrelated")
    <...>
    >>> gt.random_rewire(g, "erdos")
    >>> corr = gt.avg_neighbour_corr(g, "out", "out")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="o-", label=r"Erd\H{o}s")
    <...>
    >>> xlabel("$k$")
    <...>
    >>> ylabel(r"$\left<k_{nn}\right>$")
    <...>
    >>> legend(loc="best")
    <...>
    >>> savefig("shuffled-stats.pdf")

    .. testcode::
       :hide:

       savefig("shuffled-stats.png")


    .. figure:: shuffled-stats.*
        :align: center

        Average degree correlations for the different shuffled and non-shuffled
        graphs. The shuffled graph with correlations displays exactly the same
        correlation as the original graph.

    Now let's do it for a directed graph. See
    :func:`~graph_tool.generation.random_graph` for more details.

    >>> p = scipy.stats.poisson
    >>> g = gt.random_graph(20000, lambda: (sample_k(19), sample_k(19)),
    ...                     lambda a, b: (p.pmf(a[0], b[1]) * p.pmf(a[1], 20 - b[0])),
    ...                     mix_time=100)
    >>> figure()
    <...>
    >>> axes([0.1,0.15,0.6,0.8])
    <...>
    >>> corr = gt.avg_neighbour_corr(g, "in", "out")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="o-",
    ...          label=r"$\left<\text{o}\right>$ vs i")
    <...>
    >>> corr = gt.avg_neighbour_corr(g, "out", "in")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="o-",
    ...          label=r"$\left<\text{i}\right>$ vs o")
    <...>
    >>> gt.random_rewire(g, "correlated")
    >>> corr = gt.avg_neighbour_corr(g, "in", "out")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="o-",
    ...          label=r"$\left<\text{o}\right>$ vs i, corr.")
    <...>
    >>> corr = gt.avg_neighbour_corr(g, "out", "in")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="o-",
    ...          label=r"$\left<\text{i}\right>$ vs o, corr.")
    <...>
    >>> gt.random_rewire(g, "uncorrelated")
    >>> corr = gt.avg_neighbour_corr(g, "in", "out")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="o-",
    ...          label=r"$\left<\text{o}\right>$ vs i, uncorr.")
    <...>
    >>> corr = gt.avg_neighbour_corr(g, "out", "in")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="o-",
    ...          label=r"$\left<\text{i}\right>$ vs o, uncorr.")
    <...>
    >>> legend(bbox_to_anchor=(1.01, 0.5), loc="center left", borderaxespad=0.)
    <...>
    >>> xlabel("Source degree")
    <...>
    >>> ylabel("Average target degree")
    <...>
    >>> savefig("shuffled-deg-corr-dir.pdf")

    .. testcode::
       :hide:

       savefig("shuffled-deg-corr-dir.png")

    .. figure:: shuffled-deg-corr-dir.*
        :align: center

        Average degree correlations for the different shuffled and non-shuffled
        directed graphs. The shuffled graph with correlations displays exactly
        the same correlation as the original graph.

    References
    ----------
    .. [metropolis-equations-1953]  Metropolis, N.; Rosenbluth, A.W.;
       Rosenbluth, M.N.; Teller, A.H.; Teller, E. "Equations of State
       Calculations by Fast Computing Machines". Journal of Chemical Physics 21
       (6): 1087-1092 (1953). :doi:`10.1063/1.1699114`
    .. [hastings-monte-carlo-1970] Hastings, W.K. "Monte Carlo Sampling Methods
       Using Markov Chains and Their Applications". Biometrika 57 (1): 97-109 (1970).
       :doi:`10.1093/biomet/57.1.97`
    .. [holland-stochastic-1983] Paul W. Holland, Kathryn Blackmond Laskey, and
       Samuel Leinhardt, "Stochastic blockmodels: First steps," Social Networks
       5, no. 2: 109-13 (1983) :doi:`10.1016/0378-8733(83)90021-7`
    .. [karrer-stochastic-2011] Brian Karrer and M. E. J. Newman, "Stochastic
       blockmodels and community structure in networks," Physical Review E 83,
       no. 1: 016107 (2011) :doi:`10.1103/PhysRevE.83.016107` :arxiv:`1008.3926`

    """
    if not parallel_edges:
        p = label_parallel_edges(g)
        if p.a.max() != 0:
            raise ValueError("Parallel edge detected. Can't rewire " +
                             "graph without parallel edges if it " +
                             "already contains parallel edges!")
    if not self_loops:
        l = label_self_loops(g)
        if l.a.max() != 0:
            raise ValueError("Self-loop detected. Can't rewire graph " +
                             "without self-loops if it already contains" +
                             " self-loops!")

    if (deg_corr is not None and not g.is_directed()) and blockmodel is None:
        corr = lambda i, j: deg_corr(i[1], j[1])
    else:
        corr = deg_corr

    if strat != "probabilistic":
        g = GraphView(g, reversed=False)
    elif blockmodel is not None:
        strat = "blockmodel"
    pcount = libgraph_tool_generation.random_rewire(g._Graph__graph, strat,
                                                    n_iter, not edge_sweep,
                                                    self_loops, parallel_edges,
                                                    corr, _prop("v", g, blockmodel),
                                                    cache_probs,
                                                    _get_rng(), verbose)
    if ret_fail:
        return pcount


def predecessor_tree(g, pred_map):
    """Return a graph from a list of predecessors given by the ``pred_map`` vertex property."""

    _check_prop_scalar(pred_map, "pred_map")
    pg = Graph()
    libgraph_tool_generation.predecessor_graph(g._Graph__graph,
                                               pg._Graph__graph,
                                               _prop("v", g, pred_map))
    return pg


def line_graph(g):
    """Return the line graph of the given graph `g`.

    Notes
    -----
    Given an undirected graph G, its line graph L(G) is a graph such that

        * each vertex of L(G) represents an edge of G; and
        * two vertices of L(G) are adjacent if and only if their corresponding
          edges share a common endpoint ("are adjacent") in G.

    For a directed graph, the second criterion becomes:

       * Two vertices representing directed edges from u to v and from w to x in
         G are connected by an edge from uv to wx in the line digraph when v =
         w.


    Examples
    --------

    >>> g = gt.collection.data["lesmis"]
    >>> lg, vmap = gt.line_graph(g)
    >>> gt.graph_draw(g, pos=g.vp["pos"], output="lesmis.pdf")
    <...>
    >>> pos = gt.graph_draw(lg, output="lesmis-lg.pdf")

    .. testcode::
       :hide:

       gt.graph_draw(g, pos=g.vp["pos"], output="lesmis.png")
       pos = gt.graph_draw(lg, pos=pos, output="lesmis-lg.png")


    .. figure:: lesmis.png
       :align: left

       Coappearances of characters in Victor Hugo's novel "Les Miserables".

    .. figure:: lesmis-lg.png
       :align: right

       Line graph of the coappearance network on the left.

    References
    ----------
    .. [line-wiki] http://en.wikipedia.org/wiki/Line_graph
    """
    lg = Graph(directed=g.is_directed())

    vertex_map = lg.new_vertex_property("int64_t")

    libgraph_tool_generation.line_graph(g._Graph__graph,
                                        lg._Graph__graph,
                                        _prop("v", lg, vertex_map))
    return lg, vertex_map


def graph_union(g1, g2, intersection=None, props=None, include=False):
    """Return the union of graphs g1 and g2, composed of all edges and vertices
    of g1 and g2, without overlap.

    Parameters
    ----------
    g1 : :class:`~graph_tool.Graph`
       First graph in the union.
    g2 : :class:`~graph_tool.Graph`
       Second graph in the union.
    intersection : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
       Vertex property map owned by `g1` which maps each of each of its vertices
       to vertex indexes belonging to `g2`. Negative values mean no mapping
       exists, and thus both vertices in `g1` and `g2` will be present in the
       union graph.
    props : list of tuples of :class:`~graph_tool.PropertyMap` (optional, default: ``[]``)
       Each element in this list must be a tuple of two PropertyMap objects. The
       first element must be a property of `g1`, and the second of `g2`. The
       values of the property maps are propagated into the union graph, and
       returned.
    include : bool (optional, default: ``False``)
       If true, graph `g2` is inserted into `g1` which is modified. If false, a
       new graph is created, and both graphs remain unmodified.

    Returns
    -------
    ug : :class:`~graph_tool.Graph`
        The union graph
    props : list of :class:`~graph_tool.PropertyMap` objects
        List of propagated properties.  This is only returned if `props` is not
        empty.

    Examples
    --------

    .. testcode::
       :hide:

       from numpy.random import random, seed
       from pylab import *
       seed(42)
       gt.seed_rng(42)

    >>> g = gt.triangulation(random((300,2)))[0]
    >>> ug = gt.graph_union(g, g)
    >>> uug = gt.graph_union(g, ug)
    >>> pos = gt.sfdp_layout(g)
    >>> gt.graph_draw(g, pos=pos, output_size=(300,300), output="graph_original.pdf")
    <...>

    .. testcode::
       :hide:

       gt.graph_draw(g, pos=pos, output_size=(300,300), output="graph_original.png")

    >>> pos = gt.sfdp_layout(ug)
    >>> gt.graph_draw(ug, pos=pos, output_size=(300,300), output="graph_union.pdf")
    <...>

    .. testcode::
       :hide:

       gt.graph_draw(ug, pos=pos, output_size=(300,300), output="graph_union.png")

    >>> pos = gt.sfdp_layout(uug)
    >>> gt.graph_draw(uug, pos=pos, output_size=(300,300), output="graph_union2.pdf")
    <...>

    .. testcode::
       :hide:

       gt.graph_draw(uug, pos=pos, output_size=(300,300), output="graph_union2.png")


    .. image:: graph_original.*
    .. image:: graph_union.*
    .. image:: graph_union2.*

    """
    if props == None:
        props = []
    if not include:
        g1 = Graph(g1)
    if intersection is None:
        intersection = g1.new_vertex_property("int32_t")
        intersection.a = 0
    else:
        intersection = intersection.copy("int32_t")
        intersection.a[intersection.a >= 0] += 1
        intersection.a[intersection.a < 0] = 0

    g1.stash_filter(directed=True)
    g1.set_directed(True)
    g2.stash_filter(directed=True)
    g2.set_directed(True)
    n_props = []

    try:
        vmap, emap = libgraph_tool_generation.graph_union(g1._Graph__graph,
                                                          g2._Graph__graph,
                                                          _prop("v", g1,
                                                                intersection))
        for p1, p2 in props:
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


@_limit_args({"type": ["simple", "delaunay"]})
def triangulation(points, type="simple", periodic=False):
    r"""
    Generate a 2D or 3D triangulation graph from a given point set.

    Parameters
    ----------
    points : :class:`~numpy.ndarray`
        Point set for the triangulation. It may be either a N x d array, where N
        is the number of points, and d is the space dimension (either 2 or 3).
    type : string (optional, default: ``'simple'``)
        Type of triangulation. May be either 'simple' or 'delaunay'.
    periodic : bool (optional, default: ``False``)
        If ``True``, periodic boundary conditions will be used. This is
        parameter is valid only for type="delaunay", and is otherwise ignored.

    Returns
    -------
    triangulation_graph : :class:`~graph_tool.Graph`
        The generated graph.
    pos : :class:`~graph_tool.PropertyMap`
        Vertex property map with the Cartesian coordinates.

    See Also
    --------
    random_graph: random graph generation

    Notes
    -----

    A triangulation [cgal-triang]_ is a division of the convex hull of a point
    set into triangles, using only that set as triangle vertices.

    In simple triangulations (`type="simple"`), the insertion of a point is done
    by locating a face that contains the point, and splitting this face into
    three new faces (the order of insertion is therefore important). If the
    point falls outside the convex hull, the triangulation is restored by
    flips. Apart from the location, insertion takes a time O(1). This bound is
    only an amortized bound for points located outside the convex hull.

    Delaunay triangulations (`type="delaunay"`) have the specific empty sphere
    property, that is, the circumscribing sphere of each cell of such a
    triangulation does not contain any other vertex of the triangulation in its
    interior. These triangulations are uniquely defined except in degenerate
    cases where five points are co-spherical. Note however that the CGAL
    implementation computes a unique triangulation even in these cases.

    Examples
    --------
    .. testcode::
       :hide:

       from numpy.random import random, seed
       from pylab import *
       seed(42)
       gt.seed_rng(42)
    >>> points = random((500, 2)) * 4
    >>> g, pos = gt.triangulation(points)
    >>> weight = g.new_edge_property("double") # Edge weights corresponding to
    ...                                        # Euclidean distances
    >>> for e in g.edges():
    ...    weight[e] = sqrt(sum((array(pos[e.source()]) -
    ...                          array(pos[e.target()]))**2))
    >>> b = gt.betweenness(g, weight=weight)
    >>> b[1].a *= 100
    >>> gt.graph_draw(g, pos=pos, output_size=(300,300), vertex_fill_color=b[0],
    ...               edge_pen_width=b[1], output="triang.pdf")
    <...>

    .. testcode::
       :hide:

       gt.graph_draw(g, pos=pos, output_size=(300,300), vertex_fill_color=b[0],
                     edge_pen_width=b[1], output="triang.png")

    >>> g, pos = gt.triangulation(points, type="delaunay")
    >>> weight = g.new_edge_property("double")
    >>> for e in g.edges():
    ...    weight[e] = sqrt(sum((array(pos[e.source()]) -
    ...                          array(pos[e.target()]))**2))
    >>> b = gt.betweenness(g, weight=weight)
    >>> b[1].a *= 120
    >>> gt.graph_draw(g, pos=pos, output_size=(300,300), vertex_fill_color=b[0],
    ...               edge_pen_width=b[1], output="triang-delaunay.pdf")
    <...>

    .. testcode::
       :hide:

       gt.graph_draw(g, pos=pos, output_size=(300,300), vertex_fill_color=b[0],
                     edge_pen_width=b[1], output="triang-delaunay.png")


    2D triangulation of random points:

    .. image:: triang.*
    .. image:: triang-delaunay.*

    *Left:* Simple triangulation. *Right:* Delaunay triangulation. The vertex
    colors and the edge thickness correspond to the weighted betweenness
    centrality.

    References
    ----------
    .. [cgal-triang] http://www.cgal.org/Manual/last/doc_html/cgal_manual/Triangulation_3/Chapter_main.html

    """

    if points.shape[1] not in [2, 3]:
        raise ValueError("points array must have shape N x d, with d either 2 or 3.")
    # copy points to ensure continuity and correct data type
    points = numpy.array(points, dtype='float64')
    if points.shape[1] == 2:
        npoints = numpy.zeros((points.shape[0], 3))
        npoints[:,:2] = points
        points = npoints
    g = Graph(directed=False)
    pos = g.new_vertex_property("vector<double>")
    libgraph_tool_generation.triangulation(g._Graph__graph, points,
                                           _prop("v", g, pos), type, periodic)
    return g, pos


def lattice(shape, periodic=False):
    r"""
    Generate a N-dimensional square lattice.

    Parameters
    ----------
    shape : list or :class:`~numpy.ndarray`
        List of sizes in each dimension.
    periodic : bool (optional, default: ``False``)
        If ``True``, periodic boundary conditions will be used.

    Returns
    -------
    lattice_graph : :class:`~graph_tool.Graph`
        The generated graph.

    See Also
    --------
    triangulation: 2D or 3D triangulation
    random_graph: random graph generation

    Examples
    --------
    .. testcode::
       :hide:

       gt.seed_rng(42)

    >>> g = gt.lattice([10,10])
    >>> pos = gt.sfdp_layout(g, cooling_step=0.95, epsilon=1e-2)
    >>> gt.graph_draw(g, pos=pos, output_size=(300,300), output="lattice.pdf")
    <...>

    .. testcode::
       :hide:

       gt.graph_draw(g, pos=pos, output_size=(300,300), output="lattice.png")

    >>> g = gt.lattice([10,20], periodic=True)
    >>> pos = gt.sfdp_layout(g, cooling_step=0.95, epsilon=1e-2)
    >>> gt.graph_draw(g, pos=pos, output_size=(300,300), output="lattice_periodic.pdf")
    <...>

    .. testcode::
       :hide:

       gt.graph_draw(g, pos=pos, output_size=(300,300), output="lattice_periodic.png")

    >>> g = gt.lattice([10,10,10])
    >>> pos = gt.sfdp_layout(g, cooling_step=0.95, epsilon=1e-2)
    >>> gt.graph_draw(g, pos=pos, output_size=(300,300), output="lattice_3d.pdf")
    <...>

    .. testcode::
       :hide:

       gt.graph_draw(g, pos=pos, output_size=(300,300), output="lattice_3d.png")


    .. image:: lattice.*
    .. image:: lattice_periodic.*
    .. image:: lattice_3d.*

    *Left:* 10x10 2D lattice. *Middle:* 10x20 2D periodic lattice (torus).
    *Right:* 10x10x10 3D lattice.

    References
    ----------
    .. [lattice] http://en.wikipedia.org/wiki/Square_lattice

    """

    g = Graph(directed=False)
    libgraph_tool_generation.lattice(g._Graph__graph, shape, periodic)
    return g


def geometric_graph(points, radius, ranges=None):
    r"""
    Generate a geometric network form a set of N-dimensional points.

    Parameters
    ----------
    points : list or :class:`~numpy.ndarray`
        List of points. This must be a two-dimensional array, where the rows are
        coordinates in a N-dimensional space.
    radius : float
        Pairs of points with an euclidean distance lower than this parameters
        will be connected.
    ranges : list or :class:`~numpy.ndarray` (optional, default: ``None``)
        If provided, periodic boundary conditions will be assumed, and the
        values of this parameter it will be used as the ranges in all
        dimensions. It must be a two-dimensional array, where each row will
        cointain the lower and upper bound of each dimension.

    Returns
    -------
    geometric_graph : :class:`~graph_tool.Graph`
        The generated graph.
    pos : :class:`~graph_tool.PropertyMap`
        A vertex property map with the position of each vertex.

    Notes
    -----
    A geometric graph [geometric-graph]_ is generated by connecting points
    embedded in a N-dimensional euclidean space which are at a distance equal to
    or smaller than a given radius.

    See Also
    --------
    triangulation: 2D or 3D triangulation
    random_graph: random graph generation
    lattice : N-dimensional square lattice

    Examples
    --------
    .. testcode::
       :hide:

       from numpy.random import random, seed
       from pylab import *
       seed(42)
       gt.seed_rng(42)

    >>> points = random((500, 2)) * 4
    >>> g, pos = gt.geometric_graph(points, 0.3)
    >>> gt.graph_draw(g, pos=pos, output_size=(300,300), output="geometric.pdf")
    <...>

    .. testcode::
       :hide:

       gt.graph_draw(g, pos=pos, output_size=(300,300), output="geometric.png")

    >>> g, pos = gt.geometric_graph(points, 0.3, [(0,4), (0,4)])
    >>> pos = gt.graph_draw(g, output_size=(300,300), output="geometric_periodic.pdf")

    .. testcode::
       :hide:

       gt.graph_draw(g, pos=pos, output_size=(300,300), output="geometric_periodic.png")


    .. image:: geometric.*
    .. image:: geometric_periodic.*

    *Left:* Geometric network with random points. *Right:* Same network, but
     with periodic boundary conditions.

    References
    ----------
    .. [geometric-graph] Jesper Dall and Michael Christensen, "Random geometric
       graphs", Phys. Rev. E 66, 016121 (2002), :doi:`10.1103/PhysRevE.66.016121`

    """

    g = Graph(directed=False)
    pos = g.new_vertex_property("vector<double>")
    if type(points) != numpy.ndarray:
        points = numpy.array(points)
    if len(points.shape) < 2:
        raise ValueError("points list must be a two-dimensional array!")
    if ranges is not None:
        periodic = True
        if type(ranges) != numpy.ndarray:
            ranges = numpy.array(ranges, dtype="float")
        else:
            ranges = array(ranges, dtype="float")
    else:
        periodic = False
        ranges = ()

    libgraph_tool_generation.geometric(g._Graph__graph, points, float(radius),
                                       ranges, periodic,
                                       _prop("v", g, pos))
    return g, pos


def price_network(N, m=1, c=None, gamma=1, directed=True, seed_graph=None):
    r"""A generalized version of Price's -- or Barabási-Albert if undirected -- preferential attachment network model.

    Parameters
    ----------
    N : int
        Size of the network.
    m : int (optional, default: ``1``)
        Out-degree of newly added vertices.
    c : float (optional, default: ``1 if directed == True else 0``)
        Constant factor added to the probability of a vertex receiving an edge
        (see notes below).
    gamma : float (optional, default: ``1``)
        Preferential attachment power (see notes below).
    directed : bool (optional, default: ``True``)
        If ``True``, a Price network is generated. If ``False``, a
        Barabási-Albert network is generated.
    seed_graph : :class:`~graph_tool.Graph` (optional, default: ``None``)
        If provided, this graph will be used as the starting point of the
        algorithm.

    Returns
    -------
    price_graph : :class:`~graph_tool.Graph`
        The generated graph.

    Notes
    -----

    The (generalized) [price]_ network is either a directed or undirected graph
    (the latter is called a Barabási-Albert network), generated dynamically by
    at each step adding a new vertex, and connecting it to :math:`m` other
    vertices, chosen with probability :math:`\pi` defined as:

    .. math::

        \pi \propto k^\gamma + c

    where :math:`k` is the in-degree of the vertex (or simply the degree in the
    undirected case). If :math:`\gamma=1`, the tail of resulting in-degree
    distribution of the directed case is given by

    .. math::

        P_{k_\text{in}} \sim k_\text{in}^{-(2 + c/m)},

    or for the undirected case

    .. math::

        P_{k} \sim k^{-(3 + c/m)}.

    However, if :math:`\gamma \ne 1`, the in-degree distribution is not
    scale-free (see [dorogovtsev-evolution]_ for details).

    Note that if `seed_graph` is not given, the algorithm will *always* start
    with one node if :math:`c > 0`, or with two nodes with a link between them
    otherwise. If :math:`m > 1`, the degree of the newly added vertices will be
    vary dynamically as :math:`m'(t) = \min(m, N(t))`, where :math:`N(t)` is the
    number of vertices added so far. If this behaviour is undesired, a proper
    seed graph with :math:`N \ge m` vertices must be provided.

    This algorithm runs in :math:`O(N\log N)` time.

    See Also
    --------
    triangulation: 2D or 3D triangulation
    random_graph: random graph generation
    lattice : N-dimensional square lattice
    geometric_graph : N-dimensional geometric network

    Examples
    --------
    .. testcode::
       :hide:

       gt.seed_rng(42)

    >>> g = gt.price_network(100000)
    >>> gt.graph_draw(g, pos=gt.sfdp_layout(g, epsilon=1e-2, cooling_step=0.95),
    ...               vertex_fill_color=g.vertex_index, vertex_size=2,
    ...               edge_pen_width=1, output="price-network.png")
    <...>
    >>> g = gt.price_network(100000, c=0.1)
    >>> gt.graph_draw(g, pos=gt.sfdp_layout(g, epsilon=1e-2, cooling_step=0.95),
    ...               vertex_fill_color=g.vertex_index, vertex_size=2,
    ...               edge_pen_width=1, output="price-network-broader.png")
    <...>

    .. figure:: price-network.png
        :align: center

        Price network with :math:`N=10^5` nodes and :math:`c=1`.  The colors
        represent the order in which vertices were added.

    .. figure:: price-network-broader.png
        :align: center

        Price network with :math:`N=10^5` nodes and :math:`c=0.1`.  The colors
        represent the order in which vertices were added.


    References
    ----------

    .. [yule] Yule, G. U. "A Mathematical Theory of Evolution, based on the
       Conclusions of Dr. J. C. Willis, F.R.S.". Philosophical Transactions of
       the Royal Society of London, Ser. B 213: 21-87, 1925,
       :doi:`10.1098/rstb.1925.0002`
    .. [price] Derek De Solla Price, "A general theory of bibliometric and other
       cumulative advantage processes", Journal of the American Society for
       Information Science, Volume 27, Issue 5, pages 292-306, September 1976,
       :doi:`10.1002/asi.4630270505`
    .. [barabasi-albert] Barabási, A.-L., and Albert, R., "Emergence of
       scaling in random networks", Science, 286, 509, 1999,
       :doi:`10.1126/science.286.5439.509`
    .. [dorogovtsev-evolution] S. N. Dorogovtsev and J. F. F. Mendes, "Evolution
       of networks", Advances in Physics, 2002, Vol. 51, No. 4, 1079-1187,
       :doi:`10.1080/00018730110112519`
    """

    if c is None:
        c = 1 if directed else 0

    if seed_graph is None:
        g = Graph(directed=directed)
        if c > 0:
            g.add_vertex()
        else:
            g.add_vertex(2)
            g.add_edge(g.vertex(1), g.vertex(0))
        N -= g.num_vertices()
    else:
        g = seed_graph
    libgraph_tool_generation.price(g._Graph__graph, N, gamma, c, m, _get_rng())
    return g
