#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# graph_tool -- a general graph manipulation python module
#
# Copyright (C) 2007-2010 Tiago de Paula Peixoto <tiago@forked.de>
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

from .. dl_import import dl_import
dl_import("import libgraph_tool_generation")

from .. core import Graph, _check_prop_scalar, _prop, _limit_args
from .. stats import label_parallel_edges, label_self_loops
import sys, numpy, numpy.random

__all__ = ["random_graph", "random_rewire", "predecessor_tree", "line_graph",
           "graph_union", "triangulation", "lattice", "geometric_graph",
           "price_network"]


def random_graph(N, deg_sampler, deg_corr=None, directed=True,
                 parallel_edges=False, self_loops=False, random=True,
                 verbose=False):
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
        A function which gives the degree correlation of the graph. It should be
        callable with two parameters: the in,out-degree pair of the source
        vertex an edge, and the in,out-degree pair of the target of the same
        edge (for undirected graphs, both parameters are single values). The
        function should return a number proportional to the probability of such
        an edge existing in the generated graph.
    directed : bool (optional, default: True)
        Whether the generated graph should be directed.
    parallel_edges : bool (optional, default: False)
        If True, parallel edges are allowed.
    self_loops : bool (optional, default: False)
        If True, self-loops are allowed.
    random : bool (optional, default: True)
        If True, the returned graph is randomized.
    verbose : bool (optional, default: False)
        If True, verbose information is displayed.

    Returns
    -------
    random_graph : :class:`~graph_tool.Graph`
        The generated graph.

    See Also
    --------
    random_rewire: in place graph shuffling

    Notes
    -----
    The algorithm makes sure the degree sequence is graphical (i.e. realizable)
    and keeps re-sampling the degrees if is not. With a valid degree sequence,
    the edges are placed deterministically, and later the graph is shuffled with
    the :func:`~graph_tool.generation.random_rewire` function.

    The complexity is :math:`O(V+E)` if parallel edges are allowed, and
    :math:`O(V+E\log N_k)` if parallel edges are not allowed, where :math:`N_k <
    V` is the number of different degrees sampled (or in,out-degree pairs).

    References
    ----------
    [deg-sequence] http://en.wikipedia.org/wiki/Degree_%28graph_theory%29#Degree_sequence

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
    (0.62986894481988553, 0.011101504846821255)

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
    >>> xlabel("in-degree")
    <...>
    >>> ylabel("out-degree")
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

    >>> figure(figsize=(6,3))
    <...>
    >>> axes([0.1,0.15,0.63,0.8])
    <...>
    >>> corr = gt.avg_neighbour_corr(g, "in", "in")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="o-",
    ...         label=r"$\left<\text{in}\right>$ vs in")
    (...)
    >>> corr = gt.avg_neighbour_corr(g, "in", "out")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="o-",
    ...         label=r"$\left<\text{out}\right>$ vs in")
    (...)
    >>> corr = gt.avg_neighbour_corr(g, "out", "in")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="o-",
    ...          label=r"$\left<\text{in}\right>$ vs out")
    (...)
    >>> corr = gt.avg_neighbour_corr(g, "out", "out")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="o-",
    ...          label=r"$\left<\text{out}\right>$ vs out")
    (...)
    >>> legend(loc=(1.05,0.5))
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
    seed = numpy.random.randint(0, sys.maxint)
    g = Graph()
    if deg_corr == None:
        uncorrelated = True
    else:
        uncorrelated = False
    libgraph_tool_generation.gen_random_graph(g._Graph__graph, N, deg_sampler,
                                              uncorrelated, not parallel_edges,
                                              not self_loops, not directed,
                                              seed, verbose, True)
    g.set_directed(directed)
    if random:
        random_rewire(g, parallel_edges=parallel_edges,
                      self_loops=self_loops, verbose=verbose)
        if deg_corr != None:
            random_rewire(g, strat="probabilistic",
                          parallel_edges=parallel_edges, deg_corr=deg_corr,
                          self_loops=self_loops, verbose=verbose)
    return g


@_limit_args({"strat": ["erdos", "correlated", "uncorrelated", "probabilistic"]})
def random_rewire(g, strat="uncorrelated", parallel_edges=False,
                  self_loops=False, deg_corr=None, verbose=False):
    r"""
    Shuffle the graph in-place. If `strat` != "erdos", the degrees (either in or
    out) of each vertex are always the same, but otherwise the edges are
    randomly placed. If `strat` = "correlated", the degree correlations are
    also maintained: The new source and target of each edge both have the same
    in and out-degree. If `strat` = "probabilistic", than edges are rewired
    according to the degree correlation given by the parameter `deg_corr`.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be shuffled. The graph will be modified.
    strat : string (optional, default: "uncorrelated")
        If `strat` = "erdos", the resulting graph will be entirely random. If
        `strat` = "uncorrelated" only the degrees of the vertices will be
        maintained, nothing else. If `strat` = "correlated", additionally the
        new source and target of each edge both have the same in and out-degree.
        If `strat` = "probabilistic", than edges are rewired according to the
        degree correlation given by the parameter `deg_corr`.
    parallel : bool (optional, default: False)
        If True, parallel edges are allowed.
    self_loops : bool (optional, default: False)
        If True, self-loops are allowed.
    deg_corr : function (optional, default: None)
        A function which gives the degree correlation of the graph. It should be
        callable with two parameters: the in,out-degree pair of the source
        vertex an edge, and the in,out-degree pair of the target of the same
        edge (for undirected graphs, both parameters are single values). The
        function should return a number proportional to the probability of such
        an edge existing in the generated graph. This parameter is ignored,
        unless `strat` = "probabilistic".
    verbose : bool (optional, default: False)
        If True, verbose information is displayed.

    See Also
    --------
    random_graph: random graph generation

    Notes
    -----
    This algorithm iterates through all the edges in the network and tries to
    swap its target or source with the target or source of another edge.

    .. note::
        If `parallel_edges` = False, parallel edges are not placed during
        rewiring. In this case, for some special graphs it may be necessary to
        call the function more than once to obtain a graph which corresponds to
        a uniform sample from the ensemble. But typically, if the graph is
        sufficiently large, a single call should be enough.

    Each edge gets swapped at least once, so the overall complexity is
    :math:`O(E)`. If `strat` = "probabilistic" the complexity is
    :math:`O(E\log N_k)`,  where :math:`N_k < V` is the number of different
    degrees (or in,out-degree pairs).


    Examples
    --------

    Some small graphs for visualization.

    >>> from numpy.random import random, seed
    >>> from pylab import *
    >>> seed(42)
    >>> g, pos = gt.triangulation(random((1000,2)))
    >>> gt.graph_draw(g, layout="arf", output="rewire_orig.png", size=(6,6))
    <...>
    >>> gt.random_rewire(g, "correlated")
    >>> gt.graph_draw(g, layout="arf", output="rewire_corr.png", size=(6,6))
    <...>
    >>> gt.random_rewire(g)
    >>> gt.graph_draw(g, layout="arf", output="rewire_uncorr.png", size=(6,6))
    <...>
    >>> gt.random_rewire(g, "erdos")
    >>> gt.graph_draw(g, layout="arf", output="rewire_erdos.png", size=(6,6))
    <...>

    Some `ridiculograms <http://www.youtube.com/watch?v=YS-asmU3p_4>`_ :

    .. image:: rewire_orig.png
    .. image:: rewire_corr.png
    .. image:: rewire_uncorr.png
    .. image:: rewire_erdos.png

    *From left to right:* Original graph; Shuffled graph, with degree
    correlations; Shuffled graph, without degree correlations; Shuffled graph,
    with random degrees.

    We can try some larger graphs to get better statistics.

    >>> figure()
    <...>
    >>> g = gt.random_graph(30000, lambda: sample_k(20),
    ...                     lambda i,j: exp(abs(i-j)), directed=False)
    >>> corr = gt.avg_neighbour_corr(g, "out", "out")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="o-", label="original")
    (...)
    >>> gt.random_rewire(g, "correlated")
    >>> corr = gt.avg_neighbour_corr(g, "out", "out")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="*", label="correlated")
    (...)
    >>> gt.random_rewire(g)
    >>> corr = gt.avg_neighbour_corr(g, "out", "out")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="o-", label="uncorrelated")
    (...)
    >>> gt.random_rewire(g, "erdos")
    >>> corr = gt.avg_neighbour_corr(g, "out", "out")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="o-", label="Erdos")
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
    ...                     lambda a,b: (p.pmf(a[0],b[1])*p.pmf(a[1],20-b[0])))
    >>> figure(figsize=(6,3))
    <...>
    >>> axes([0.1,0.15,0.6,0.8])
    <...>
    >>> corr = gt.avg_neighbour_corr(g, "in", "out")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="o-",
    ...          label=r"$\left<\text{o}\right>$ vs i")
    (...)
    >>> corr = gt.avg_neighbour_corr(g, "out", "in")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="o-",
    ...          label=r"$\left<\text{i}\right>$ vs o")
    (...)
    >>> gt.random_rewire(g, "correlated")
    >>> corr = gt.avg_neighbour_corr(g, "in", "out")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="o-",
    ...          label=r"$\left<\text{o}\right>$ vs i, corr.")
    (...)
    >>> corr = gt.avg_neighbour_corr(g, "out", "in")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="o-",
    ...          label=r"$\left<\text{i}\right>$ vs o, corr.")
    (...)
    >>> gt.random_rewire(g, "uncorrelated")
    >>> corr = gt.avg_neighbour_corr(g, "in", "out")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="o-",
    ...          label=r"$\left<\text{o}\right>$ vs i, uncorr.")
    (...)
    >>> corr = gt.avg_neighbour_corr(g, "out", "in")
    >>> errorbar(corr[2][:-1], corr[0], yerr=corr[1], fmt="o-",
    ...          label=r"$\left<\text{i}\right>$ vs o, uncorr.")
    (...)
    >>> legend(loc=(1.05,0.45))
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

    seed = numpy.random.randint(0, sys.maxint)

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

    if deg_corr != None and  not g.is_directed():
        corr = lambda i, j: deg_corr(i[1], j[1])
    else:
        corr = deg_corr

    if corr == None:
        g.stash_filter(reversed=True)
    try:
        libgraph_tool_generation.random_rewire(g._Graph__graph, strat,
                                               self_loops, parallel_edges,
                                               corr, seed, verbose)
    finally:
        if corr == None:
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


def graph_union(g1, g2, props=None, include=False):
    """Return the union of graphs g1 and g2, composed of all edges and vertices
    of g1 and g2, without overlap.

    Parameters
    ----------
    g1 : :class:`~graph_tool.Graph`
       First graph in the union.
    g2 : :class:`~graph_tool.Graph`
       Second graph in the union.
    props : list of tuples of :class:`~graph_tool.PropertyMap` (optional, default: [])
       Each element in this list must be a tuple of two PropertyMap objects. The
       first element must be a property of `g1`, and the second of `g2`. The
       values of the property maps are propagated into the union graph, and
       returned.
    include : bool (optional, default: False)
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

    >>> from numpy.random import random, seed
    >>> seed(42)
    >>> g = gt.triangulation(random((300,2)))[0]
    >>> ug = gt.graph_union(g, g)
    >>> uug = gt.graph_union(g, ug)
    >>> gt.graph_draw(g, layout="arf", size=(8,8), output="graph_original.png")
    <...>
    >>> gt.graph_draw(ug, layout="arf", size=(8,8), output="graph_union.png")
    <...>
    >>> gt.graph_draw(uug, layout="arf", size=(8,8), output="graph_union2.png")
    <...>

    .. image:: graph_original.png
    .. image:: graph_union.png
    .. image:: graph_union2.png

    """
    if props == None:
        props = []
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


@_limit_args({"type": ["simple", "delaunay"]})
def triangulation(points, type="simple", periodic=False):
    r"""
    Generate a 2D or 3D triangulation graph from a given point set.

    Parameters
    ----------
    points : :class:`~numpy.ndarray`
        Point set for the triangulation. It may be either a N x d array, where N
        is the number of points, and d is the space dimension (either 2 or 3).
    type : string (optional, default: 'simple')
        Type of triangulation. May be either 'simple' or 'delaunay'.
    periodic : bool (optional, default: False)
        If True, periodic boundary conditions will be used. This is parameter is
        valid only for type="delaunay", and is otherwise ignored.

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
    >>> from numpy.random import seed, random
    >>> seed(42)
    >>> points = random((500, 2)) * 4
    >>> g, pos = gt.triangulation(points)
    >>> weight = g.new_edge_property("double") # Edge weights corresponding to
    ...                                        # Euclidean distances
    >>> for e in g.edges():
    ...    weight[e] = sqrt(sum((array(pos[e.source()]) -
    ...                          array(pos[e.target()]))**2))
    >>> b = gt.betweenness(g, weight=weight)
    >>> b[1].a *= 100
    >>> gt.graph_draw(g, pos=pos, pin=True, size=(8,8), vsize=0.07, vcolor=b[0],
    ...               eprops={"penwidth":b[1]}, output="triang.png")
    <...>
    >>> g, pos = gt.triangulation(points, type="delaunay")
    >>> weight = g.new_edge_property("double")
    >>> for e in g.edges():
    ...    weight[e] = sqrt(sum((array(pos[e.source()]) -
    ...                          array(pos[e.target()]))**2))
    >>> b = gt.betweenness(g, weight=weight)
    >>> b[1].a *= 120
    >>> gt.graph_draw(g, pos=pos, pin=True, size=(8,8), vsize=0.07, vcolor=b[0],
    ...               eprops={"penwidth":b[1]}, output="triang-delaunay.png")
    <...>

    2D triangulation of random points:

    .. image:: triang.png
    .. image:: triang-delaunay.png

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
    periodic : bool (optional, default: False)
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
    >>> g = gt.lattice([10,10])
    >>> gt.graph_draw(g, size=(8,8), output="lattice.png")
    <...>
    >>> g = gt.lattice([10,20], periodic=True)
    >>> gt.graph_draw(g, size=(8,8), output="lattice_periodic.png")
    <...>
    >>> g = gt.lattice([10,10,10])
    >>> gt.graph_draw(g, size=(8,8), output="lattice_3d.png")
    <...>

    .. image:: lattice.png
    .. image:: lattice_periodic.png
    .. image:: lattice_3d.png

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
    ranges : list or :class:`~numpy.ndarray` (optional, default: None)
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
    >>> from numpy.random import seed, random
    >>> seed(42)
    >>> points = random((500, 2)) * 4
    >>> g, pos = gt.geometric_graph(points, 0.3)
    >>> gt.graph_draw(g, pos=pos, pin=True, size=(8,8), output="geometric.png")
    <...>
    >>> g, pos = gt.geometric_graph(points, 0.3, [(0,4), (0,4)])
    >>> gt.graph_draw(g, size=(8,8), output="geometric_periodic.png")
    <...>

    .. image:: geometric.png
    .. image:: geometric_periodic.png

    *Left:* Geometric network with random points. *Right:* Same network, but
     with periodic boundary conditions.

    References
    ----------
    .. [geometric-graph] Jesper Dall and Michael Christensen, "Random geometric
       graphs", Phys. Rev. E 66, 016121 (2002), DOI: 10.1103/PhysRevE.66.016121

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
    m : int (optional, default: 1)
        Out-degree of newly added vertices.
    c : float (optional, default: 1 if ``directed==True`` else 0)
        Constant factor added to the probability of a vertex receiving an edge
        (see notes below).
    gamma : float (optional, default: 1)
        Preferential attachment power (see notes below).
    directed : bool (optional, default: True)
        If ``True``, a Price network is generated. If ``False``, a
        Barabási-Albert network is generated.
    seed_graph : :class:`~graph_tool.Graph` (optional, default: None)
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
    vertices, chosen with probability:

    .. math::

        P \propto k^\gamma + c

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

    This algorithm runs in :math:`O(N\log N)` time.

    See Also
    --------
    triangulation: 2D or 3D triangulation
    random_graph: random graph generation
    lattice : N-dimensional square lattice
    geometric_graph : N-dimensional geometric network

    Examples
    --------
    >>> from numpy.random import seed, random
    >>> seed(42)
    >>> g = gt.price_network(100000)
    >>> gt.graph_draw(g, layout="sfdp", size=(12,12), vcolor=g.vertex_index,
    ...               output="price-network.png")
    <...>
    >>> g = gt.price_network(100000, c=0.1)
    >>> gt.graph_draw(g, layout="sfdp", size=(12,12), vcolor=g.vertex_index,
    ...               output="price-network-broader.png")
    <...>

    .. image:: price-network.png
    .. image:: price-network-broader.png

    Price networks with :math:`N=10^5` nodes. **Left:** :math:`c=1`, **Right:**
    :math:`c=0.1`. The colors represent the order in which vertices were
    added.

    References
    ----------

    .. [yule] Yule, G. U. "A Mathematical Theory of Evolution, based on the
       Conclusions of Dr. J. C. Willis, F.R.S.". Philosophical Transactions of
       the Royal Society of London, Ser. B 213: 21–87, 1925,
       doi:10.1098/rstb.1925.0002
    .. [price] Derek De Solla Price, "A general theory of bibliometric and other
       cumulative advantage processes", Journal of the American Society for
       Information Science, Volume 27, Issue 5, pages 292–306, September 1976,
       DOI: 10.1002/asi.4630270505
    .. [barabasi-albert] Barabási, A.-L., and Albert, R., "Emergence of
       scaling in random networks", Science, 286, 509, 1999;
       DOI: 10.1126/science.286.5439.509
    .. [dorogovtsev-evolution] S. N. Dorogovtsev and J. F. F. Mendes, "Evolution
       of networks", Advances in Physics, 2002, Vol. 51, No. 4, 1079-1187,
       DOI: 10.1080/00018730110112519
    """

    if c is None:
        c = 1 if directed else 0

    if seed_graph is None:
        if directed:
            g = Graph()
            g.add_vertex(m)
        else:
            N_s = m + 1 if m % 2 != 0 else m + 2
            g = random_graph(N_s, lambda: 1, directed=False)
        N -= g.num_vertices()
    else:
        g = seed_graph
        if g.num_vertices() < m:
            raise ValueError("seed_graph has number of vertices < m!")
    seed = numpy.random.randint(0, sys.maxint)
    libgraph_tool_generation.price(g._Graph__graph, N, gamma, c, m, seed)
    return g
