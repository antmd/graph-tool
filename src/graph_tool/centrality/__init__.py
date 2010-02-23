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
``graph_tool.centrality`` - Centrality measures
-----------------------------------------------

This module includes centrality-related algorithms.

Summary
+++++++

.. autosummary::
   :nosignatures:

   pagerank
   betweenness
   central_point_dominance
   eigentrust
   absolute_trust

Contents
++++++++
"""

from .. dl_import import dl_import
dl_import("import libgraph_tool_centrality")

from .. core import _prop
import sys, numpy

__all__ = ["pagerank", "betweenness", "central_point_dominance", "eigentrust",
           "absolute_trust"]

def pagerank(g, damping=0.8, prop=None, epslon=1e-6, max_iter=None,
             ret_iter=False):
    r"""
    Calculate the PageRank of each vertex.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    damping : float, optional (default: 0.8)
        Damping factor.
    prop : :class:`~graph_tool.PropertyMap`, optional (default: None)
        Vertex property map to store the PageRank values.
    epslon : float, optional (default: 1e-6)
        Convergence condition. The iteration will stop if the total delta of all
        vertices are below this value.
    max_iter : int, optional (default: None)
        If supplied, this will limit the total number of iterations.
    ret_iter : bool, optional (default: False)
        If true, the total number of iterations is also returned.

    Returns
    -------
    pagerank : :class:`~graph_tool.PropertyMap`
        A vertex property map containing the PageRank values.

    See Also
    --------
    betweenness: betweenness centrality
    eigentrust: eigentrust centrality
    absolute_trust: absolute trust centrality

    Notes
    -----
    The value of PageRank [pagerank-wikipedia]_ of vertex v :math:`PR(v)` is
    given interactively by the relation:

    .. math::

        PR(v) = \frac{1-d}{N} + d \sum_{w \in \Gamma^{-}(v)}
                \frac{PR (w)}{d^{+}(w)}

    where :math:`\Gamma^{-}(v)` are the in-neighbours of v, :math:`d^{+}(w)` is
    the out-degree of w, and d is a damping factor.

    The implemented algorithm progressively iterates the above condition, until
    it no longer changes, according to the parameter epslon. It has a
    topology-dependent running time.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    >>> from numpy.random import poisson, seed
    >>> seed(42)
    >>> g = gt.random_graph(100, lambda: (poisson(3), poisson(3)))
    >>> pr = gt.pagerank(g)
    >>> print pr.a
    [ 0.63876901  1.13528868  0.31465963  0.55855277  0.2         0.75605741
      0.42628689  0.53066254  0.55004112  0.91717076  0.71164749  0.32015438
      0.67275227  1.08207389  1.14412231  0.9049167   1.32002     1.4692142
      0.76549771  0.71510277  0.23732927  0.40844911  0.2         0.27912876
      0.71309781  0.32015438  1.3376236   0.31352887  0.59346569  0.33381039
      0.67300081  0.73318264  0.65812653  0.73409673  0.93051993  0.83241145
      1.59816568  0.43979363  0.2512247   1.15663357  0.2         0.35977148
      0.72182022  1.01267711  0.76304859  0.49247376  0.49384283  1.8436647
      0.64312224  1.00778243  0.62287633  1.15215387  0.56176895  0.7166227
      0.56506109  0.67104337  0.95570565  0.27996953  0.79975983  0.33631497
      1.09471419  0.33631497  0.2512247   2.09126732  0.68157485  0.2
      0.37140185  0.65619459  1.27370737  0.48383225  1.36125161  0.2
      0.78300573  1.03427279  0.56904755  1.66077917  1.73302035  0.28749261
      0.83143045  1.04969728  0.70090048  0.55991433  0.68440994  0.2
      0.34018009  0.45485484  0.28        1.2015438   2.11850885  1.24990775
      0.59914308  0.59989185  0.73535564  0.78168417  0.55390281  0.38627667
      1.42274704  0.51105348  0.92550979  1.27968065]

    References
    ----------
    .. [pagerank-wikipedia] http://en.wikipedia.org/wiki/Pagerank
    .. [lawrence-pagerank-1998] P. Lawrence, B. Sergey, M. Rajeev, W. Terry,
       "The pagerank citation ranking: Bringing order to the web", Technical
       report, Stanford University, 1998
    """

    if max_iter == None:
        max_iter = 0
    if prop == None:
        prop = g.new_vertex_property("double")
    ic = libgraph_tool_centrality.\
            get_pagerank(g._Graph__graph, _prop("v", g, prop), damping, epslon,
                         max_iter)
    if ret_iter:
        return prop, ic
    else:
        return prop

def betweenness(g, vprop=None, eprop=None, weight=None, norm=True):
    r"""
    Calculate the betweenness centrality for each vertex and edge.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    vprop : :class:`~graph_tool.PropertyMap`, optional (default: None)
        Vertex property map to store the vertex betweenness values.
    eprop : :class:`~graph_tool.PropertyMap`, optional (default: None)
        Edge property map to store the edge betweenness values.
    weight : :class:`~graph_tool.PropertyMap`, optional (default: None)
        Edge property map corresponding to the weight value of each edge.
    norm : bool, optional (default: True)
        Whether or not the betweenness values should be normalized.

    Returns
    -------
    vertex_betweenness : A vertex property map with the vertex betweenness
                         values.
    edge_betweenness : An edge property map with the edge betweenness
                       values.

    See Also
    --------
    central_point_dominance: central point dominance of the graph
    pagerank: PageRank centrality
    eigentrust: eigentrust centrality
    absolute_trust: absolute trust centrality

    Notes
    -----
    Betweenness centrality of a vertex :math:`C_B(v)` is defined as,

    .. math::

        C_B(v)= \sum_{s \neq v \neq t \in V \atop s \neq t}
                \frac{\sigma_{st}(v)}{\sigma_{st}}

    where :math:`\sigma_{st}` is the number of shortest geodesic paths from s to
    t, and :math:`\sigma_{st}(v)` is the number of shortest geodesic paths from
    s to t that pass through a vertex v.  This may be normalised by dividing
    through the number of pairs of vertices not including v, which is
    :math:`(n-1)(n-2)/2`.

    The algorithm used here is defined in [brandes-faster-2001]_, and has a
    complexity of :math:`O(VE)` for unweighted graphs and :math:`O(VE + V(V+E)
    \log V)` for weighted graphs. The space complexity is :math:`O(VE)`.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    >>> from numpy.random import poisson, seed
    >>> seed(42)
    >>> g = gt.random_graph(100, lambda: (poisson(3), poisson(3)))
    >>> vb, eb = gt.betweenness(g)
    >>> print vb.a
    [ 0.03395047  0.07911989  0.00702948  0.02337119  0.          0.02930099
      0.01684377  0.02558675  0.03440095  0.02886187  0.03124262  0.00975953
      0.01307953  0.03938858  0.07266505  0.01313647  0.          0.06450598
      0.0575418   0.00525468  0.00466089  0.01803829  0.          0.00050161
      0.0085034   0.02362432  0.05620574  0.00097157  0.04006816  0.01301474
      0.02154916  0.          0.06009194  0.02780363  0.08963522  0.04049657
      0.06993559  0.02082698  0.00288318  0.03264322  0.          0.03641759
      0.01083859  0.03750864  0.04079359  0.02092599  0.          0.02153655
      0.          0.05674631  0.03861911  0.05473282  0.00904367  0.03249097
      0.00894043  0.0192741   0.03379204  0.02125998  0.0018321   0.0013495
      0.0336502   0.0210088   0.00125318  0.0489189   0.05254974  0.
      0.00432189  0.04866168  0.06444727  0.02508525  0.02533085  0.
      0.05308703  0.02539854  0.02270809  0.044889    0.04766016  0.0086368
      0.01501699  0.          0.03107868  0.0054221   0.          0.
      0.00596081  0.01183977  0.00159761  0.11435876  0.03988501  0.05128991
      0.04558135  0.02303469  0.05092032  0.04700221  0.00927644  0.00841903
      0.          0.03243633  0.04514374  0.05170213]

    References
    ----------
    .. [betweenness-wikipedia] http://en.wikipedia.org/wiki/Centrality#Betweenness_centrality
    .. [brandes-faster-2001] U. Brandes, "A faster algorithm for betweenness
       centrality",  Journal of Mathematical Sociology, 2001
    """
    if vprop == None:
        vprop = g.new_vertex_property("double")
    if eprop == None:
        eprop = g.new_edge_property("double")
    if weight != None and weight.value_type() != eprop.value_type():
        nw = g.new_edge_property(eprop.value_type())
        g.copy_property(weight, nw)
        weight = nw
    libgraph_tool_centrality.\
            get_betweenness(g._Graph__graph, _prop("e", g, weight),
                            _prop("e", g, eprop), _prop("v", g, vprop), norm)
    return vprop, eprop

def central_point_dominance(g, betweenness):
    r"""
    Calculate the central point dominance of the graph, given the betweenness
    centrality of each vertex.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    betweenness : :class:`~graph_tool.PropertyMap`
        Vertex property map with the betweenness centrality values. The values
        must be normalized.

    Returns
    -------
    cp : float
        The central point dominance.

    See Also
    --------
    betweenness: betweenness centrality

    Notes
    -----
    Let :math:`v^*` be the vertex with the largest relative betweenness
    centrality; then, the central point dominance [freeman-set-1977]_ is defined
    as:

    .. math::

        C'_B = \frac{1}{|V|-1} \sum_{v} C_B(v^*) - C_B(v)

    where :math:`C_B(v)` is the normalized betweenness centrality of vertex
    v. The value of :math:`C_B` lies in the range [0,1].

    The algorithm has a complexity of :math:`O(V)`.

    Examples
    --------
    >>> from numpy.random import poisson, seed
    >>> seed(42)
    >>> g = gt.random_graph(100, lambda: (poisson(3), poisson(3)))
    >>> vb, eb = gt.betweenness(g)
    >>> print gt.central_point_dominance(g, vb)
    0.0884414811909

    References
    ----------
    .. [freeman-set-1977] Linton C. Freeman, "A Set of Measures of Centrality
       Based on Betweenness", Sociometry, Vol. 40, No. 1,  pp. 35-41 (1977)
    """

    return libgraph_tool_centrality.\
           get_central_point_dominance(g._Graph__graph,
                                       _prop("v", g, betweenness))


def eigentrust(g, trust_map, vprop=None, norm=False, epslon=1e-6, max_iter=0,
               ret_iter=False):
    r"""
    Calculate the eigentrust centrality of each vertex in the graph.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    trust_map : :class:`~graph_tool.PropertyMap`
        Edge property map with the values of trust associated with each
        edge. The values must lie in the range [0,1].
    vprop : PropertyMap, optional (default: None)
        Vertex property map where the values of eigentrust must be stored.
    norm : bool, optional (default: false)
        Norm eigentrust values so that the total sum equals 1.
    epslon : float, optional (default: 1e-6)
        Convergence condition. The iteration will stop if the total delta of all
        vertices are below this value.
    max_iter : int, optional (default: None)
        If supplied, this will limit the total number of iterations.
    ret_iter : bool, optional (default: False)
        If true, the total number of iterations is also returned.

    Returns
    -------
    eigentrust : A vertex property map containing the eigentrust values.

    See Also
    --------
    betweenness: betweenness centrality
    pagerank: PageRank centrality
    absolute_trust: absolute trust centrality

    Notes
    -----
    The eigentrust [kamvar-eigentrust-2003]_ values :math:`t_i` correspond the
    following limit

    .. math::

        \mathbf{t} = \lim_{n\to\infty} \left(C^T\right)^n \mathbf{c}

    where :math:`c_i = 1/|V|` and the elements of the matrix :math:`C` are the
    normalized trust values:

    .. math::

        c_{ij} = \frac{\max(s_{ij},0)}{\sum_{j} \max(s_{ij}, 0)}

    The algorithm has a topology-dependent complexity.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    >>> from numpy.random import poisson, random, seed
    >>> seed(42)
    >>> g = gt.random_graph(100, lambda: (poisson(3), poisson(3)))
    >>> trust = g.new_edge_property("double")
    >>> trust.get_array()[:] = random(g.num_edges())*42
    >>> t = gt.eigentrust(g, trust, norm=True)
    >>> print t.get_array()
    [  5.51422638e-03   1.12397965e-02   2.34959294e-04   6.32738574e-03
       0.00000000e+00   6.34804836e-03   2.67885424e-03   4.02497751e-03
       1.67943467e-02   6.46196106e-03   1.92402451e-02   9.04032352e-04
       9.70843104e-03   1.40319816e-02   1.04995777e-02   2.86712231e-02
       2.47285894e-02   2.38394469e-02   7.06936059e-03   9.45794717e-03
       2.09970054e-05   1.64768298e-03   0.00000000e+00   1.19346706e-03
       6.88434371e-03   5.36337333e-03   2.08428677e-02   2.85813783e-03
       1.10564670e-02   3.16345060e-04   5.25737238e-03   5.43761445e-03
       7.98048389e-03   7.95939648e-03   2.23891858e-02   5.68630666e-03
       2.09300588e-02   4.28902068e-03   1.70833078e-03   2.37814042e-02
       0.00000000e+00   1.20805010e-03   1.29713483e-02   5.73021992e-03
       8.71093674e-03   7.77661067e-03   8.76489806e-04   2.38519385e-02
       3.53225723e-03   8.46948906e-03   5.09874234e-03   2.44547150e-02
       1.32342629e-02   1.80085559e-03   4.37189381e-03   1.18195253e-02
       1.62748861e-02   1.83200678e-04   1.09745025e-02   1.47544090e-03
       3.34512926e-02   1.58885132e-03   1.13128910e-03   3.04944830e-02
       4.22684975e-03   0.00000000e+00   9.89654274e-04   4.25927156e-03
       2.34516214e-02   4.91370905e-03   2.29366664e-02   0.00000000e+00
       6.83407601e-03   1.60508753e-02   1.62762068e-03   3.94324856e-02
       2.84109571e-02   8.81167727e-04   2.16999908e-02   1.28688125e-02
       1.10825963e-02   2.64915564e-03   2.88711928e-03   0.00000000e+00
       4.24392252e-03   9.38398819e-03   0.00000000e+00   1.74508371e-02
       3.26594153e-02   4.07188867e-02   3.20678152e-03   6.35046287e-03
       8.07061556e-03   5.08505374e-03   3.27300367e-03   3.30989070e-03
       2.30651195e-02   4.20338525e-03   5.04332662e-03   3.58731532e-02]

    References
    ----------
    .. [kamvar-eigentrust-2003] S. D. Kamvar, M. T. Schlosser, H. Garcia-Molina
       "The eigentrust algorithm for reputation management in p2p networks",
       Proceedings of the 12th international conference on World Wide Web,
       Pages: 640 - 651, 2003
    """

    if vprop == None:
        vprop = g.new_vertex_property("double")
    i = libgraph_tool_centrality.\
           get_eigentrust(g._Graph__graph, _prop("e", g, trust_map),
                          _prop("v", g, vprop), epslon, max_iter)
    if norm:
        vprop.get_array()[:] /= sum(vprop.get_array())

    if ret_iter:
        return vprop, i
    else:
        return vprop

def absolute_trust(g, trust_map, source, target = None, vprop=None):
    r"""
    Calculate the absolute trust centrality of each vertex in the graph, from a
    given source.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    trust_map : :class:`~graph_tool.PropertyMap`
        Edge property map with the values of trust associated with each
        edge. The values must lie in the range [0,1].
    source : Vertex
        Source vertex. All trust values are computed relative to this vertex.
    target : Vertex (optional, default: None)
        The only target for which the trust value will be calculated. If left
        unspecified, the trust values for all targets are computed.
    vprop : :class:`~graph_tool.PropertyMap` (optional, default: None)
        A vertex property map where the values of trust for each source
        must be stored.

    Returns
    -------
    absolute_trust : :class:`~graph_tool.PropertyMap` or float
        A vertex property map containing the absolute trust vector from the
        source vertex to the rest of the network. If `target` is specified, the
        result is a single float, with the corresponding trust value for the
        target.

    See Also
    --------
    eigentrust: eigentrust centrality
    betweenness: betweenness centrality
    pagerank: PageRank centrality

    Notes
    -----
    The absolute trust between vertices i and j is defined as

    .. math::

        t_{ij} = \frac{\sum_m A_{m,j} w^2_{G\setminus\{j\}}(i\to m)c_{m,j}}
                 {\sum_m A_{m,j} w_{G\setminus\{j\}}(i\to m)}

    where :math:`A_{ij}` is the adjacency matrix, :math:`c_{ij}` is the direct
    trust from i to j, and :math:`w_G(i\to j)` is the weight of the path with
    maximum weight from i to j, computed as

    .. math::

       w_G(i\to j) = \prod_{e\in i\to j} c_e.

    The algorithm measures the absolute trust by finding the paths with maximum
    weight, using Dijkstra's algorithm, to all in-neighbours of a given
    target. This search needs to be performed repeatedly for every target, since
    it needs to be removed from the graph first. The resulting complexity is
    therefore :math:`O(N^2\log N)` for all targets, and :math:`O(N\log N)` for a
    single target.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    >>> from numpy.random import poisson, random, seed
    >>> seed(42)
    >>> g = gt.random_graph(100, lambda: (poisson(3), poisson(3)))
    >>> trust = g.new_edge_property("double")
    >>> trust.a = random(g.num_edges())
    >>> t = gt.absolute_trust(g, trust, source=g.vertex(0))
    >>> print t.a
    [ 0.16260667  0.04129912  0.13735376  0.19146125  0.          0.09147461
      0.10371912  0.12465511  0.24631221  0.0603916   0.2375385   0.06637879
      0.08897662  0.0800988   0.05250601  0.66759022  0.09368793  0.08275437
      0.13674709  0.15553915  0.01376162  0.417068    0.          0.06096886
      0.08746817  0.39380693  0.09215297  0.09575144  0.15594162  0.04008874
      0.05483972  0.05691086  0.13571077  0.32376012  0.22477937  0.06347962
      0.10445085  0.19447845  0.38007043  0.13810585  0.          0.08451096
      0.06648153  0.18479174  0.13003649  0.14850631  0.00320603  0.1074644
      0.12088162  0.06792678  0.08472666  0.2002143   0.25963204  0.37838425
      0.03089371  0.18389694  0.39420339  0.03348093  0.11483196  0.0656204
      0.14206403  0.07066434  0.25168986  0.07040126  0.04870569  0.
      0.09861349  0.03882069  0.1105267   0.07951823  0.08748441  0.
      0.08393443  0.11121719  0.21903223  0.25529628  0.0414386   0.03695558
      0.17664854  0.05143033  0.11735779  0.06525968  0.19600919  0.          0.1220922
      0.33330041  0.          0.28595961  0.14526678  0.12514885  0.089524
      0.40738962  0.03719195  0.54409979  0.06247424  0.10660136  0.11674385
      0.13218144  0.02214988  0.23215937]
    """

    if vprop == None:
        vprop = g.new_vertex_property("double")

    source = g.vertex_index[source]

    if target == None:
        target = -1
    else:
        target = g.vertex_index[target]

    libgraph_tool_centrality.\
            get_absolute_trust(g._Graph__graph, source, target,
                               _prop("e", g, trust_map), _prop("v", g, vprop))
    if target != -1:
        return vprop.a[target]
    return vprop

