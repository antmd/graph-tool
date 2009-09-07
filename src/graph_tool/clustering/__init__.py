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
``graph_tool.clustering`` - Clustering coefficients
---------------------------------------------------

Provides algorithms for calculation of clustering coefficients,
aka. transitivity.
"""

from .. dl_import import dl_import
dl_import("import libgraph_tool_clustering as _gt")

from .. core import _degree, _prop, Graph
from .. topology import isomorphism
from .. generation import random_rewire
from itertools import izip
from numpy import *
import sys

__all__ = ["local_clustering", "global_clustering", "extended_clustering",
           "motifs", "motif_significance"]

def local_clustering(g, prop=None, undirected=False):
    r"""
    Return vertex property containing local clustering coefficients for all
    vertices.

    Parameters
    ----------
    g : Graph
        Graph to be used.
    prop : PropertyMap or string, optional
        Vertex property map where results will be stored. If specified, this
        parameter will also be the return value.
    undirected : bool, optional
        Calculate the *undirected* clustering coefficient, if graph is directed
        (this option has no effect if the graph is undirected).

    Returns
    -------
    prop : PropertyMap
        Vertex property containing the clustering coefficients.

    See Also
    --------
    global_clustering: global clustering coefficient
    extended_clustering: extended (generalized) clustering coefficient
    motifs: motif counting

    Notes
    -----
    The local clustering coefficient [1]_ :math:`c_i` is defined as

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
    >>> g = gt.random_graph(1000, lambda: (5,5), seed=42)
    >>> clust = gt.local_clustering(g)
    >>> print gt.vertex_average(g, clust)
    (0.0042016666666666669, 0.00041579094306313759)

    References
    ----------
    .. [1] D. J. Watts and Steven Strogatz, "Collective dynamics of
       'small-world' networks", Nature, vol. 393, pp 440-442, 1998.
       doi:10.1038/30918
    """

    if prop == None:
        prop = g.new_vertex_property("double")
    was_directed = g.is_directed()
    if g.is_directed() and undirected:
        g.set_directed(False)
    try:
       _gt.extended_clustering(g._Graph__graph,
                                [_prop("v", g, prop)])
    finally:
        if was_directed and undirected:
            g.set_directed(True)
    return prop

def global_clustering(g):
    r"""
    Return global clustering coefficients for graphs.

    Parameters
    ----------
    g : Graph
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
    The global clustering coefficient [1]_ :math:`c` is defined as

    .. math::
       c = 3 \times \frac{\text{number of triangles}}
                          {\text{number of connected triples}}

    The implemented algorithm runs in :math:`O(|V|\left< k\right>^3)` time,
    where :math:`\left< k\right>` is the average (total) degree.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    >>> g = gt.random_graph(1000, lambda: (5,5), seed=42)
    >>> print gt.global_clustering(g)
    (0.0083567321834469854, 0.0004417198625120785)

    References
    ----------
    .. [1] M. E. J. Newman, "The structure and function of complex networks",
       SIAM Review, vol. 45, pp. 167-256, 2003
    """

    c =_gt.global_clustering(g._Graph__graph)
    return c

def extended_clustering(g, props=None, max_depth=3, undirected=False):
    r"""
    Return a list of vertex properties containing the extended clustering
    coefficients for all vertices.

    Parameters
    ----------
    g : Graph
        Graph to be used.
    props : list of PropertyMap objects, optional
        list of vertex property maps where results will be stored. If specified,
        this parameter will also be the return value.
    max_depth : int, optional
        Maximum clustering order (default: 3).
    undirected : bool, optional
        Calculate the *undirected* clustering coefficients, if graph is directed
        (this option has no effect if the graph is undirected).

    Returns
    -------
    prop : list of PropertyMap objects
        List of vertex properties containing the clustering coefficients.

    See Also
    --------
    local_clustering: local clustering coefficient
    global_clustering: global clustering coefficient
    motifs: motif counting

    Notes
    -----
    The definition of the extended clustering coefficient :math:`c^d_i` of order
    :math:`d` is defined as

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

    The implemented algorithm runs in :math:`O(|V|\left<
    k\right>^{2+\text{max_depth}})` worst time, where :math:`\left< k\right>` is
    the average out-degree.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    >>> g = gt.random_graph(1000, lambda: (5,5), seed=42)
    >>> clusts = gt.extended_clustering(g, max_depth=5)
    >>> for i in xrange(0, 5):
    ...    print gt.vertex_average(g, clusts[i])
    ...
    (0.0042016666666666669, 0.00041579094306313759)
    (0.02409, 0.00095845344754511225)
    (0.11065333333333334, 0.0019404812805074933)
    (0.40180499999999997, 0.0031866048246265693)
    (0.43999999999999995, 0.0031172994010129269)

    References
    ----------
    .. [1] A. H. Abdo, A. P. S. de Moura, "Clustering as a measure of the local
       topology of networks", arXiv:physics/0605235v4
    """

    was_directed = g.is_directed()
    if g.is_directed() and undirected:
        g.set_directed(False)
    if props == None:
        props = []
        for i in xrange(0, max_depth):
            props.append(g.new_vertex_property("double"))
    try:
       _gt.extended_clustering(g._Graph__graph,
                               [_prop("v", g, p) for p in props])
    finally:
        if was_directed and undirected:
            g.set_directed(True)
    return props


def motifs(g, k, p=1.0, motif_list=None, undirected=None, seed=0):
    r"""
    Count the occurrence of k-size subgraphs (motifs). A tuple with two lists is
    returned: the list of motifs found, and the list with their respective
    counts.

    Parameters
    ----------
    g : Graph
        Graph to be used.
    k : int
        number of vertices of the motifs
    p : float or float list, optional (default: 1.0)
        uniform fraction of the motifs to be sampled. If a float list is
        provided, it will be used as the fraction at each depth
        :math:`[1,\dots,k]` in the algorithm. See [wernicke_efficient_2006]_ for
        more details.
    motif_list : list of Graph objects, optional
        If supplied, the algorithms will only search for the motifs in this list
        (or isomorphisms thereof)
    undirected : bool, optional
        Treat the graph as *undirected*, if graph is directed
        (this option has no effect if the graph is undirected).
    seed : int, optional (default: 0)
        Seed for the random number generator. It the value is 0, a random seed
        is used.

    Returns
    -------
    motifs : list of Graph objects
        List of motifs of size k found in the Graph. Graphs are grouped
        according to their isomorphism class, and only one of each class appears
        in this list. The list is sorted according to in-degree sequence,
        out-degree-sequence, and number of edges (in this order).
    counts : list of ints
        The number of times the respective motif in the motifs list was counted

    See Also
    --------
    motif_significance: significance profile of motifs
    local_clustering: local clustering coefficient
    global_clustering: global clustering coefficient
    extended_clustering: extended (generalized) clustering coefficient

    Notes
    -----
    This functions implements the ESU and RAND-ESU algorithms described in
    [wernicke_efficient_2006]_.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    >>> g = gt.random_graph(1000, lambda: (5,5), seed=42)
    >>> motifs, counts = gt.motifs(g, 4, undirected=True)
    >>> print len(motifs)
    10
    >>> print counts
    [116138, 392344, 389, 445, 2938, 996, 792, 3, 14, 3]

    References
    ----------
    .. [wernicke_efficient_2006] S. Wernicke, "Efficient detection of network
       motifs", IEEE/ACM Transactions on Computational Biology and
       Bioinformatics (TCBB), Volume 3, Issue 4, Pages 347-359, 2006.
    """

    if seed == 0:
        seed = random.randint(0, sys.maxint)

    sub_list = []
    directed_motifs = g.is_directed() if undirected == None else not undirected

    if motif_list != None:
        directed_motifs = motif_list[0].is_directed()
        for m in motif_list:
            if m.is_directed() != directed_motifs:
                raise ValueError("all motif graphs must be either directed or undirected")
            if m.num_vertices() != k:
                raise ValueError("all motifs must have the same number of vertices: " + k)
            sub_list.append(m._Graph__graph)

    if type(p) == float:
        pd = [1.0]*(k-1)
        pd.append(p)
    if type(p) == list:
        pd = [float(x) for x in p]

    hist = []
    was_directed = g.is_directed()
    if g.is_directed() and not directed_motifs:
        g.set_directed(False)
    try:
       _gt.get_motifs(g._Graph__graph, k, sub_list, hist, pd,
                      True, len(sub_list) == 0,
                      seed)
    finally:
        if was_directed and not directed_motifs:
            g.set_directed(True)

    # assemble graphs
    temp = []
    for m in sub_list:
        mg = Graph()
        mg._Graph__graph = m
        temp.append(mg)
    sub_list = temp

    list_hist = zip(sub_list, hist)
    # sort according to in-degree sequence
    list_hist.sort(lambda x,y: cmp(sorted([v.in_degree() for v in x[0].vertices()]),
                                   sorted([v.in_degree() for v in y[0].vertices()])))

    # sort according to out-degree sequence
    list_hist.sort(lambda x,y: cmp(sorted([v.out_degree() for v in x[0].vertices()]),
                                   sorted([v.out_degree() for v in y[0].vertices()])))

    # sort according to ascending number of edges
    list_hist.sort(lambda x,y: cmp(x[0].num_edges(), y[0].num_edges()))

    sub_list = [x[0] for x in list_hist]
    hist = [x[1] for x in list_hist]

    return sub_list, hist

def motif_significance(g, k, n_shuffles=100, p=1.0, motif_list=None,
                       undirected=None, self_loops=False, parallel_edges=False,
                       full_output=False, shuffle_strategy="uncorrelated",
                       seed=0):
    r"""
    Obtain the motif significance profile, for subgraphs with k vertices. A
    tuple with two lists is returned: the list of motifs found, and their
    respective z-scores.

    Parameters
    ----------
    g : Graph
        Graph to be used.
    k : int
        number of vertices of the motifs
    n_shuffles : int, optional (default: 100)
        number of shuffled networks to consider for the z-score
    p : float or float list, optional (default: 1.0)
        uniform fraction of the motifs to be sampled. If a float list is
        provided, it will be used as the fraction at each depth
        :math:`[1,\dots,k]` in the algorithm. See [wernicke_efficient_2006]_ for
        more details.
    motif_list : list of Graph objects, optional
        If supplied, the algorithms will only search for the motifs in this list
        (or isomorphisms thereof)
    undirected : bool, optional
        Treat the graph as *undirected*, if graph is directed
        (this option has no effect if the graph is undirected).
    self_loops : bool, optional (default: False)
        Whether or not the shuffled graphs are allowed to contain self-loops
    parallel_edges : bool, optional (default: False)
        Whether or not the shuffled graphs are allowed to contain parallel
        edges.
    full_output : bool, optional (default: False)
        If set to True, three additional lists are returned: the count
        of each motif, the average count of each motif in the shuffled networks,
        and the standard deviation of the average count of each motif in the
        shuffled networks.
    shuffle_strategy : string, optional (default: "uncorrelated")
        Shuffle strategy to use. Can be either "correlated" or "uncorrelated".
        See random_rewire() for details.
    seed : int, optional (default: 0)
        Seed for the random number generator. It the value is 0, a random seed
        is used.

    Returns
    -------
    motifs : list of Graph objects
        List of motifs of size k found in the Graph. Graphs are grouped
        according to their isomorphism class, and only one of each class appears
        in this list. The list is sorted according to in-degree sequence,
        out-degree-sequence, and number of edges (in this order).
    z-scores : list of floats
        The z-score of the respective motives. See below for the definition of
        a z-score.

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
    many standard deviations is each motif count, in respect to a ensemble of
    randomly shuffled graphs with the same degree sequence.

    The z-scores values are not normalized.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    >>> g = gt.random_graph(100, lambda: (3,3), seed=10)
    >>> motifs, zscores = gt.motif_significance(g, 3)
    >>> print len(motifs)
    12
    >>> print zscores
    [1.6197267471265697, 1.6271265351937325, 0.99350347553112339, -1.4729502635694927, -1.1004922497513825, -1.1181853811005562, 0.92339606894139736, -0.11, -0.10000000000000001, -0.28999999999999998, -0.22, -0.01]
    """

    s_ms, counts = motifs(g, k, p, motif_list, undirected, seed)
    s_counts = [0]*len(s_ms)
    s_dev = [0]*len(s_ms)

    # get samples
    sg = g.copy()
    for i in xrange(0, n_shuffles):
        random_rewire(sg, shuffle_strategy, self_loops=self_loops,
                      parallel_edges=parallel_edges)
        m_temp, count_temp = motifs(sg, k, p, motif_list, undirected, seed)
        for j in xrange(0, len(m_temp)):
            found = False
            for l in xrange(0, len(s_ms)):
                if isomorphism(s_ms[l], m_temp[j]):
                    found = True
                    s_counts[l] += count_temp[j]
                    s_dev[l] += count_temp[j]**2
            if not found:
                s_ms.append(m_temp[j])
                s_counts.append(count_temp[j])
                s_dev.append(count_temp[j]**2)
                counts.append(0)

    s_counts = [ x/float(n_shuffles) for x in s_counts ]
    s_dev = [ max(sqrt(x[0]/float(n_shuffles) - x[1]**2),1) \
              for x in izip(s_dev,s_counts) ]

    list_hist = zip(s_ms, s_counts, s_dev)
    # sort according to in-degree sequence
    list_hist.sort(lambda x,y: cmp(sorted([v.in_degree()\
                                           for v in x[0].vertices()]),
                                   sorted([v.in_degree()\
                                           for v in y[0].vertices()])))

    # sort according to out-degree sequence
    list_hist.sort(lambda x,y: cmp(sorted([v.out_degree()\
                                           for v in x[0].vertices()]),
                                   sorted([v.out_degree()\
                                           for v in y[0].vertices()])))

    # sort according to ascending number of edges
    list_hist.sort(lambda x,y: cmp(x[0].num_edges(), y[0].num_edges()))

    s_ms, s_counts, s_dev = zip(*list_hist)

    zscore = [(x[0] - x[1])/x[2] for x in izip(counts, s_counts, s_dev)]

    if full_output:
        return s_ms, zscore, counts, s_counts, s_dev
    else:
        return s_ms, zscore
