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
from .. misc import isomorphism
from numpy import *
import sys

__all__ = ["local_clustering", "global_clustering", "extended_clustering",
           "motifs"]

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
    (0.0045633333333333333, 0.00041406305209606802)

    References
    ----------
    .. [1] D. J. Watts and Steven Strogatz, "Collective dynamics of
       'small-world' networks", Nature, vol. 393, pp 440-442, 1998.
       doi:10.1038/30918
    """

    if prop == None:
        prop = g.new_vertex_property("double")
    was_directed = g.directed()
    if g.directed() and undirected:
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
    (0.0086380072318200073, 0.00044516087903790925)

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
       c^d_i = \frac{\left|\right\{ \{u,v\}; u,v \in N_i | d_{G(V\diagdown
       \{i\})}(u,v) = d \left\}\right|}{\binom{\left|N_i\right|}{2}},

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
    (0.0045633333333333333, 0.00041406305209606802)
    (0.027705, 0.0010493633929938454)
    (0.11730666666666667, 0.00201118990760307)
    (0.41394666666666663, 0.0030157036105470745)
    (0.41717499999999996, 0.0030272310298907366)

    References
    ----------
    .. [1] A. H. Abdo, A. P. S. de Moura, "Clustering as a measure of the local
       topology of networks", arXiv:physics/0605235v4
    """

    was_directed = g.directed()
    if g.directed() and undirected:
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


def motifs(g, k, p=1.0, motifs=None, undirected=None, seed=0):
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
    motifs : list of Graph objects, optional
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
    >>> motifs, counts = gt.motifs(g, 4)
    >>> print len(motifs)
    19
    >>> print counts
    [1499, 203, 1881, 1022, 8, 7, 5, 18, 6, 10, 4, 14, 10, 3, 9, 6, 15, 5, 27]

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

    if motifs != None:
        directed_motifs = motifs[0].is_directed()
        for m in motifs:
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
