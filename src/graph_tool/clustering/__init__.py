#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# graph_tool -- a general graph manipulation python module
#
# Copyright (C) 2006-2015 Tiago de Paula Peixoto <tiago@skewed.de>
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
``graph_tool.clustering`` - Clustering coefficients
---------------------------------------------------

Provides algorithms for calculation of clustering coefficients,
aka. transitivity.

Summary
+++++++

.. autosummary::
   :nosignatures:

   local_clustering
   global_clustering
   extended_clustering
   motifs
   motif_significance

Contents
++++++++
"""

from __future__ import division, absolute_import, print_function

from .. dl_import import dl_import
dl_import("from . import libgraph_tool_clustering as _gt")

from .. import _degree, _prop, Graph, GraphView, PropertyMap, _get_rng
from .. topology import isomorphism
from .. generation import random_rewire
from .. stats import vertex_hist

from collections import defaultdict
from numpy import *
from numpy import random
import sys

__all__ = ["local_clustering", "global_clustering", "extended_clustering",
           "motifs", "motif_significance"]


def local_clustering(g, prop=None, undirected=True):
    r"""
    Return the local clustering coefficients for all vertices.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    prop : :class:`~graph_tool.PropertyMap` or string, optional
        Vertex property map where results will be stored. If specified, this
        parameter will also be the return value.
    undirected : bool (default: True)
        Calculate the *undirected* clustering coefficient, if graph is directed
        (this option has no effect if the graph is undirected).

    Returns
    -------
    prop : :class:`~graph_tool.PropertyMap`
        Vertex property containing the clustering coefficients.

    See Also
    --------
    global_clustering: global clustering coefficient
    extended_clustering: extended (generalized) clustering coefficient
    motifs: motif counting

    Notes
    -----
    The local clustering coefficient [watts-collective-1998]_ :math:`c_i` is
    defined as

    .. math::
       c_i = \frac{|\{e_{jk}\}|}{k_i(k_i-1)} :\, v_j,v_k \in N_i,\, e_{jk} \in E

    where :math:`k_i` is the out-degree of vertex :math:`i`, and

    .. math::
       N_i = \{v_j : e_{ij} \in E\}

    is the set of out-neighbours of vertex :math:`i`. For undirected graphs the
    value of :math:`c_i` is normalized as

    .. math::
       c'_i = 2c_i.

    The implemented algorithm runs in :math:`O(|V|\left< k\right>^3)` time,
    where :math:`\left< k\right>` is the average out-degree.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    .. testcode::
       :hide:

       np.random.seed(42)
       gt.seed_rng(42)

    >>> g = gt.random_graph(1000, lambda: (5,5))
    >>> clust = gt.local_clustering(g)
    >>> print(gt.vertex_average(g, clust))
    (0.006177777777777778, 0.0003696618407994121)

    References
    ----------
    .. [watts-collective-1998] D. J. Watts and Steven Strogatz, "Collective
       dynamics of 'small-world' networks", Nature, vol. 393, pp 440-442, 1998.
       :doi:`10.1038/30918`
    """

    if prop == None:
        prop = g.new_vertex_property("double")
    if g.is_directed() and undirected:
        g = GraphView(g, directed=False)
    _gt.local_clustering(g._Graph__graph, _prop("v", g, prop))
    return prop


def global_clustering(g):
    r"""
    Return the global clustering coefficient.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.

    Returns
    -------
    c : tuple of floats
        Global clustering coefficient and standard deviation (jacknife method)

    See Also
    --------
    local_clustering: local clustering coefficient
    extended_clustering: extended (generalized) clustering coefficient
    motifs: motif counting

    Notes
    -----
    The global clustering coefficient [newman-structure-2003]_ :math:`c` is
    defined as

    .. math::
       c = 3 \times \frac{\text{number of triangles}}
                          {\text{number of connected triples}}

    The implemented algorithm runs in :math:`O(|V|\left< k\right>^3)` time,
    where :math:`\left< k\right>` is the average (total) degree.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    .. testcode::
       :hide:

       np.random.seed(42)
       gt.seed_rng(42)

    >>> g = gt.random_graph(1000, lambda: (5,5))
    >>> print(gt.global_clustering(g))
    (0.006177777777777778, 0.0003700318726720911)

    References
    ----------
    .. [newman-structure-2003] M. E. J. Newman, "The structure and function of
       complex networks", SIAM Review, vol. 45, pp. 167-256, 2003,
       :doi:`10.1137/S003614450342480`
    """

    c = _gt.global_clustering(g._Graph__graph)
    return c


def extended_clustering(g, props=None, max_depth=3, undirected=False):
    r"""
    Return the extended clustering coefficients for all vertices.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    props : list of :class:`~graph_tool.PropertyMap` objects, optional
        list of vertex property maps where results will be stored. If specified,
        this parameter will also be the return value.
    max_depth : int, optional
        Maximum clustering order (default: 3).
    undirected : bool, optional
        Calculate the *undirected* clustering coefficients, if graph is directed
        (this option has no effect if the graph is undirected).

    Returns
    -------
    prop : list of :class:`~graph_tool.PropertyMap` objects
        List of vertex properties containing the clustering coefficients.

    See Also
    --------
    local_clustering: local clustering coefficient
    global_clustering: global clustering coefficient
    motifs: motif counting

    Notes
    -----
    The extended clustering coefficient :math:`c^d_i` of order :math:`d` is
    defined as

    .. math::
       c^d_i = \frac{\left|\right\{ \{u,v\}; u,v \in N_i | d_{G(V\setminus
       \{i\})}(u,v) = d \left\}\right|}{{\left|N_i\right| \choose 2}},

    where :math:`d_G(u,v)` is the shortest distance from vertex :math:`u` to
    :math:`v` in graph :math:`G`, and

    .. math::
       N_i = \{v_j : e_{ij} \in E\}

    is the set of out-neighbours of :math:`i`. According to the above
    definition, we have that the traditional local clustering coefficient is
    recovered for :math:`d=1`, i.e., :math:`c^1_i = c_i`.

    The implemented algorithm runs in
    :math:`O(|V|\left<k\right>^{2+\text{max-depth}})` worst time, where
    :math:`\left< k\right>` is the average out-degree.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    .. testcode::
       :hide:

       np.random.seed(42)
       gt.seed_rng(42)

    >>> g = gt.random_graph(1000, lambda: (5,5))
    >>> clusts = gt.extended_clustering(g, max_depth=5)
    >>> for i in range(0, 5):
    ...    print(gt.vertex_average(g, clusts[i]))
    ...
    (0.00421, 0.00041685103920811916)
    (0.027226666666666666, 0.0010073522830778825)
    (0.11549166666666667, 0.002081399085417734)
    (0.41280666666666666, 0.0029479221231987185)
    (0.4205716666666666, 0.003156465741141506)

    References
    ----------
    .. [abdo-clustering] A. H. Abdo, A. P. S. de Moura, "Clustering as a
       measure of the local topology of networks", :arxiv:`physics/0605235`
    """

    if g.is_directed() and undirected:
        g = GraphView(g, directed=False)
    if props == None:
        props = []
        for i in range(0, max_depth):
            props.append(g.new_vertex_property("double"))
    _gt.extended_clustering(g._Graph__graph,
                            [_prop("v", g, p) for p in props])
    return props


def motifs(g, k, p=1.0, motif_list=None, return_maps=False):
    r"""
    Count the occurrence of k-size node-induced subgraphs (motifs). A tuple with
    two lists is returned: the list of motifs found, and the list with their
    respective counts.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    k : int
        number of vertices of the motifs
    p : float or float list (optional, default: `1.0`)
        uniform fraction of the motifs to be sampled. If a float list is
        provided, it will be used as the fraction at each depth
        :math:`[1,\dots,k]` in the algorithm. See [wernicke-efficient-2006]_ for
        more details.
    motif_list : list of :class:`~graph_tool.Graph` objects, optional
        If supplied, the algorithms will only search for the motifs in this list
        (or isomorphisms).
    return_maps : bool (optional, default `False`)
        If ``True``, a list will be returned, which provides for each motif graph a
        list of vertex property maps which map the motif to its location in the
        main graph.

    Returns
    -------
    motifs : list of :class:`~graph_tool.Graph` objects
        List of motifs of size k found in the Graph. Graphs are grouped
        according to their isomorphism class, and only one of each class appears
        in this list. The list is sorted according to in-degree sequence,
        out-degree-sequence, and number of edges (in this order).
    counts : list of ints
        The number of times the respective motif in the motifs list was counted
    vertex_maps : list of lists of :class:`~graph_tool.PropertyMap` objects
        List for each motif graph containing the locations in the main
        graph. This is only returned if `return_maps == True`.

    See Also
    --------
    motif_significance: significance profile of motifs
    local_clustering: local clustering coefficient
    global_clustering: global clustering coefficient
    extended_clustering: extended (generalized) clustering coefficient

    Notes
    -----
    This functions implements the ESU and RAND-ESU algorithms described in
    [wernicke-efficient-2006]_.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    .. testcode::
       :hide:

       np.random.seed(42)
       gt.seed_rng(42)

    >>> g = gt.random_graph(1000, lambda: (5,5))
    >>> motifs, counts = gt.motifs(gt.GraphView(g, directed=False), 4)
    >>> print(len(motifs))
    52
    >>> print(counts)
    [29259, 29111, 28420, 28618, 163659, 96696, 96821, 31366, 295, 162, 185, 257, 58, 74, 562, 300, 80, 240, 273, 112, 379, 167, 329, 536, 605, 641, 176, 540, 340, 317, 298, 87, 306, 275, 252, 9, 3, 11, 5, 9, 7, 2, 2, 1, 3, 1, 2, 1, 2, 1, 1, 1]


    References
    ----------
    .. [wernicke-efficient-2006] S. Wernicke, "Efficient detection of network
       motifs", IEEE/ACM Transactions on Computational Biology and
       Bioinformatics (TCBB), Volume 3, Issue 4, Pages 347-359, 2006.
       :doi:`10.1109/TCBB.2006.51`
    .. [induced-subgraph-isomorphism] http://en.wikipedia.org/wiki/Induced_subgraph_isomorphism_problem
    """

    sub_list = []
    directed_motifs = g.is_directed()

    if motif_list is not None:
        directed_motifs = motif_list[0].is_directed()
        for m in motif_list:
            if m.is_directed() != directed_motifs:
                raise ValueError("all motif graphs must be either directed or undirected!")
            if m.num_vertices() != k:
                raise ValueError("all motifs must have the same number of vertices: " + k)
            sub_list.append(m._Graph__graph)

    if directed_motifs != g.is_directed():
        raise ValueError("motifs do not have the same directionality as the graph itself!")

    if type(p) == float:
        pd = [1.0] * (k - 1)
        pd.append(p)
    if type(p) == list:
        pd = [float(x) for x in p]

    hist = []
    vertex_maps = []
    was_directed = g.is_directed()
    _gt.get_motifs(g._Graph__graph, k, sub_list, hist, vertex_maps,
                   return_maps,  pd, True, len(sub_list) == 0,
                   _get_rng())

    # assemble graphs
    temp = []
    for m in sub_list:
        mg = Graph()
        mg._Graph__graph = m
        mg.reindex_edges()
        temp.append(mg)
    sub_list = temp

    if return_maps:
        list_hist = list(zip(sub_list, hist, vertex_maps))
    else:
        list_hist = list(zip(sub_list, hist))
    # sort according to in-degree sequence
    list_hist.sort(key=lambda x: sorted([v.in_degree() for v in x[0].vertices()]))

    # sort according to out-degree sequence
    list_hist.sort(key=lambda x: sorted([v.out_degree() for v in x[0].vertices()]))

    # sort according to ascending number of edges
    list_hist.sort(key=lambda x: x[0].num_edges())

    sub_list = [x[0] for x in list_hist]
    hist = [x[1] for x in list_hist]

    if return_maps:
        vertex_maps = [x[2] for x in list_hist]
        for i, vlist in enumerate(vertex_maps):
            sub = sub_list[i]
            vertex_maps[i] = [PropertyMap(vm, sub, "v") for vm in vlist]
        return sub_list, hist, vertex_maps
    return sub_list, hist


def _graph_sig(g):
    """return the graph signature, i.e., the in and out degree distribution as
    concatenated as a tuple."""
    bins = list(range(0, g.num_vertices() + 1))
    in_dist = vertex_hist(g, "in", bins=bins if g.is_directed() else [0],
                          float_count=False)
    out_dist = vertex_hist(g, "out", bins=bins, float_count=False)
    sig = tuple([(in_dist[1][i], in_dist[0][i]) for \
                 i in range(len(in_dist[0]))] +
                [(out_dist[1][i], out_dist[0][i]) for\
                 i in range(len(out_dist[0]))])
    return sig


def motif_significance(g, k, n_shuffles=100, p=1.0, motif_list=None,
                       threshold=0, self_loops=False, parallel_edges=False,
                       full_output=False, shuffle_model="uncorrelated"):
    r"""
    Obtain the motif significance profile, for subgraphs with k vertices. A
    tuple with two lists is returned: the list of motifs found, and their
    respective z-scores.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    k : int
        Number of vertices of the motifs
    n_shuffles : int (optional, default: 100)
        Number of shuffled networks to consider for the z-score
    p : float or float list (optional, default: 1.0)
        Uniform fraction of the motifs to be sampled. If a float list is
        provided, it will be used as the fraction at each depth
        :math:`[1,\dots,k]` in the algorithm. See [wernicke-efficient-2006]_ for
        more details.
    motif_list : list of :class:`~graph_tool.Graph` objects (optional, default: None)
        If supplied, the algorithms will only search for the motifs in this list
        (isomorphisms)
    threshold : int (optional, default: 0)
        If a given motif count is below this level, it is not considered.
    self_loops : bool (optional, default: False)
        Whether or not the shuffled graphs are allowed to contain self-loops
    parallel_edges : bool (optional, default: False)
        Whether or not the shuffled graphs are allowed to contain parallel
        edges.
    full_output : bool (optional, default: False)
        If set to True, three additional lists are returned: the count
        of each motif, the average count of each motif in the shuffled networks,
        and the standard deviation of the average count of each motif in the
        shuffled networks.
    shuffle_model : string (optional, default: "uncorrelated")
        Shuffle model to use. See :func:`~graph_tool.generation.random_rewire`
        for details.

    Returns
    -------
    motifs : list of :class:`~graph_tool.Graph` objects
        List of motifs of size k found in the Graph. Graphs are grouped
        according to their isomorphism class, and only one of each class appears
        in this list. The list is sorted according to in-degree sequence,
        out-degree-sequence, and number of edges (in this order).
    z-scores : list of floats
        The z-score of the respective motives. See below for the definition of
        the z-score.

    See Also
    --------
    motifs: motif counting or sampling
    local_clustering: local clustering coefficient
    global_clustering: global clustering coefficient
    extended_clustering: extended (generalized) clustering coefficient

    Notes
    -----
    The z-score :math:`z_i` of motif i is defined as

    .. math::
         z_i = \frac{N_i - \left<N^s_i\right>}
         {\sqrt{\left<(N^s_i)^2\right> - \left<N^s_i\right>^2}},

    where :math:`N_i` is the number of times motif i found, and :math:`N^s_i`
    is the count of the same motif but on a shuffled network. It measures how
    many standard deviations is each motif count, in respect to an ensemble of
    randomly shuffled graphs with the same degree sequence.

    The z-scores values are not normalized.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    >>> from numpy import random
    >>> random.seed(10)
    >>> g = gt.random_graph(100, lambda: (3,3))
    >>> motifs, zscores = gt.motif_significance(g, 3)
    >>> print(len(motifs))
    11
    >>> print(zscores)
    [0.04162670293683441, 0.046769781223679897, 0.55888261369689918, 0.82906624302416043, -0.41384527710386287, -0.42070845824900477, -1.3262411347517733, 1.82, -0.13, -0.37, -0.22]
    """

    s_ms, counts = motifs(g, k, p, motif_list)
    if threshold > 0:
        s_ms, counts = list(zip(*[x for x in zip(s_ms, counts) if x[1] > threshold]))
        s_ms = list(s_ms)
        counts = list(counts)
    s_counts = [0] * len(s_ms)
    s_dev = [0] * len(s_ms)

    # group subgraphs by number of edges
    m_e = defaultdict(lambda: [])
    for i in range(len(s_ms)):
        m_e[_graph_sig(s_ms[i])].append(i)

    # get samples
    sg = g.copy()
    for i in range(0, n_shuffles):
        random_rewire(sg, model=shuffle_model, self_loops=self_loops,
                      parallel_edges=parallel_edges)
        m_temp, count_temp = motifs(sg, k, p, motif_list)
        if threshold > 0:
            m_temp, count_temp = list(zip(*[x for x in zip(m_temp, count_temp) \
                                       if x[1] > threshold]))
        for j in range(0, len(m_temp)):
            found = False
            for l in m_e[_graph_sig(m_temp[j])]:
                if isomorphism(s_ms[l], m_temp[j]):
                    found = True
                    s_counts[l] += count_temp[j]
                    s_dev[l] += count_temp[j] ** 2
            if not found:
                s_ms.append(m_temp[j])
                s_counts.append(count_temp[j])
                s_dev.append(count_temp[j] ** 2)
                counts.append(0)
                m_e[_graph_sig(m_temp[j])].append(len(s_ms) - 1)

    s_counts = [x / float(n_shuffles) for x in s_counts]
    s_dev = [max(sqrt(x[0] / float(n_shuffles) - x[1] ** 2), 1) \
              for x in zip(s_dev, s_counts)]

    list_hist = list(zip(s_ms, s_counts, s_dev))
    # sort according to in-degree sequence
    list_hist.sort(key = lambda x: sorted([v.in_degree() for v in x[0].vertices()])),

    # sort according to out-degree sequence
    list_hist.sort(key = lambda x: sorted([v.out_degree() for v in x[0].vertices()]))

    # sort according to ascending number of edges
    list_hist.sort(key = lambda x: x[0].num_edges())

    s_ms, s_counts, s_dev = list(zip(*list_hist))

    zscore = [(x[0] - x[1]) / x[2] for x in zip(counts, s_counts, s_dev)]

    if full_output:
        return s_ms, zscore, counts, s_counts, s_dev
    else:
        return s_ms, zscore
