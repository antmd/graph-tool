#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# graph_tool -- a general graph manipulation python module
#
# Copyright (C) 2007-2011 Tiago de Paula Peixoto <tiago@skewed.de>
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
   eigenvector
   eigentrust
   trust_transitivity

Contents
++++++++
"""

from .. dl_import import dl_import
dl_import("import libgraph_tool_centrality")

from .. import _prop, ungroup_vector_property
import sys
import numpy

__all__ = ["pagerank", "betweenness", "central_point_dominance", "eigentrust",
           "eigenvector", "trust_transitivity"]


def pagerank(g, damping=0.85, pers=None, weight=None, prop=None, epsilon=1e-6,
             max_iter=None, ret_iter=False):
    r"""
    Calculate the PageRank of each vertex.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    damping : float, optional (default: 0.85)
        Damping factor.
    pers : :class:`~graph_tool.PropertyMap`, optional (default: None)
        Personalization vector. If omitted, a constant value of :math:`1/N`
        will be used.
    weight : :class:`~graph_tool.PropertyMap`, optional (default: None)
        Edge weights. If omitted, a constant value of 1 will be used.
    prop : :class:`~graph_tool.PropertyMap`, optional (default: None)
        Vertex property map to store the PageRank values.
    epsilon : float, optional (default: 1e-6)
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
    trust_transitivity: pervasive trust transitivity

    Notes
    -----
    The value of PageRank [pagerank-wikipedia]_ of vertex v, :math:`PR(v)`, is
    given iteratively by the relation:

    .. math::

        PR(v) = \frac{1-d}{N} + d \sum_{u \in \Gamma^{-}(v)}
                \frac{PR (u)}{d^{+}(u)}

    where :math:`\Gamma^{-}(v)` are the in-neighbours of v, :math:`d^{+}(w)` is
    the out-degree of w, and d is a damping factor.

    If a personalization property :math:`p(v)` is given, the definition becomes:

    .. math::

        PR(v) = (1-d)p(v) + d \sum_{u \in \Gamma^{-}(v)}
                \frac{PR (u)}{d^{+}(u)}

    If edge weights are also given, the equation is then generalized to:

    .. math::

        PR(v) = (1-d)p(v) + d \sum_{u \in \Gamma^{-}(v)}
                \frac{PR (u) w_{u\to v}}{d^{+}(u)}

    where :math:`d^{+}(u)=\sum_{y}A_{u,y}w_{u\to y}` is redefined to be the sum
    of the weights of the out-going edges from u.

    The implemented algorithm progressively iterates the above equations, until
    it no longer changes, according to the parameter epsilon. It has a
    topology-dependent running time.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    >>> from numpy.random import random, poisson, seed
    >>> seed(42)
    >>> g = gt.random_graph(100, lambda: (poisson(3), poisson(3)))
    >>> pr = gt.pagerank(g)
    >>> print pr.a
    [ 0.00782362  0.01642353  0.00420484  0.0038825   0.0015      0.01145378
      0.00514203  0.00593481  0.00743705  0.00785063  0.00446447  0.00440222
      0.00684158  0.00463226  0.00518308  0.0056288   0.01207045  0.00617264
      0.00958574  0.00817165  0.01041552  0.00508079  0.0015      0.00249411
      0.00842537  0.00293099  0.00873296  0.001755    0.003371    0.00817938
      0.00406813  0.00576584  0.01188752  0.00674565  0.00758134  0.00855306
      0.00975204  0.00823918  0.00209855  0.00753858  0.0015      0.001925
      0.00593262  0.00603431  0.00977679  0.00707922  0.00529399  0.01048882
      0.001755    0.0111949   0.0032813   0.01591077  0.00407595  0.01015827
      0.00383036  0.01024311  0.00714593  0.00379142  0.00955729  0.001925
      0.00737848  0.00352088  0.00654273  0.00676324  0.00353259  0.0015
      0.00809045  0.00864939  0.00626611  0.00632213  0.00939761  0.0015
      0.00584767  0.0077272   0.00688094  0.01010526  0.01071083  0.00550524
      0.0045327   0.00577072  0.00337711  0.00637928  0.01295484  0.0015
      0.00265875  0.003245    0.00203456  0.00969993  0.00908983  0.00759961
      0.00428542  0.00674196  0.0043264   0.01339053  0.00570051  0.00253539
      0.01464169  0.00505055  0.01919599  0.01413612]

    Now with a personalization vector, and edge weights:

    >>> w = g.new_edge_property("double")
    >>> w.a = random(g.num_edges())
    >>> p = g.new_vertex_property("double")
    >>> p.a = random(g.num_vertices())
    >>> p.a /= p.a.sum()
    >>> pr = gt.pagerank(g, pers=p, weight=w)
    >>> print pr.a
    [ 0.01693559  0.01316915  0.00369907  0.00245658  0.00092715  0.01380721
      0.00703909  0.00407121  0.00816254  0.00880131  0.0035886   0.0050914
      0.00815843  0.00624021  0.0069828   0.00647311  0.01260669  0.00884083
      0.01324534  0.01103024  0.01417902  0.00309344  0.00250025  0.00153889
      0.00969556  0.00491575  0.00552323  0.00300698  0.00327355  0.00829017
      0.00274335  0.00440865  0.01436394  0.00671045  0.00788395  0.01092875
      0.0126331   0.00789263  0.00422443  0.00745144  0.00148972  0.00198663
      0.00476339  0.00800871  0.01468149  0.00971962  0.00446663  0.01333257
      0.00085768  0.01044298  0.00286075  0.02119469  0.00406517  0.01317145
      0.00280023  0.0143227   0.00867722  0.00234863  0.01180399  0.00298827
      0.0049022   0.00532752  0.00603759  0.00766617  0.00293739  0.00238803
      0.00863735  0.01110095  0.00660816  0.00170262  0.00884469  0.00300867
      0.00441168  0.00630793  0.00424727  0.00906709  0.0135949   0.00890726
      0.00267835  0.00615783  0.0045653   0.00720592  0.00996495  0.0009367
      0.00233309  0.00265909  0.00211686  0.01277934  0.01284484  0.00625721
      0.00487027  0.00852522  0.00403389  0.01817233  0.00573321  0.0038696
      0.00932334  0.00515806  0.01601592  0.0167547 ]

    References
    ----------
    .. [pagerank-wikipedia] http://en.wikipedia.org/wiki/Pagerank
    .. [lawrence-pagerank-1998] P. Lawrence, B. Sergey, M. Rajeev, W. Terry,
       "The pagerank citation ranking: Bringing order to the web", Technical
       report, Stanford University, 1998
    .. [Langville-survey-2005] A. N. Langville, C. D. Meyer, "A Survey of
       Eigenvector Methods for Web Information Retrieval", SIAM Review, vol. 47,
       no. 1, pp. 135-161, 2005, :DOI:`10.1137/S0036144503424786`
    """

    if max_iter == None:
        max_iter = 0
    if prop == None:
        prop = g.new_vertex_property("double")
    ic = libgraph_tool_centrality.\
            get_pagerank(g._Graph__graph, _prop("v", g, prop),
                         _prop("v", g, pers), _prop("e", g, weight),
                         damping, epsilon, max_iter)
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
    vertex_betweenness : A vertex property map with the vertex betweenness values.
    edge_betweenness : An edge property map with the edge betweenness values.

    See Also
    --------
    central_point_dominance: central point dominance of the graph
    pagerank: PageRank centrality
    eigentrust: eigentrust centrality
    trust_transitivity: pervasive trust transitivity

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
    [  2.65012897e-02   1.04414799e-01   2.73374899e-02   1.52782183e-02
       0.00000000e+00   2.74548352e-02   3.54680121e-02   3.72671558e-02
       2.39732112e-02   2.34942149e-02   2.97950758e-02   4.08351383e-02
       4.31702840e-02   1.90317902e-02   3.66879750e-02   8.65571818e-03
       0.00000000e+00   3.74046494e-02   4.22428130e-02   2.10503176e-02
       1.39558854e-02   8.40349783e-03   0.00000000e+00   4.45784374e-03
       3.38671970e-02   1.72390157e-02   4.82232543e-02   1.03071532e-04
       1.42200266e-02   4.82793598e-02   1.82020235e-02   0.00000000e+00
       7.04969679e-02   2.31267158e-02   6.42817952e-02   3.71139131e-02
       3.81618985e-02   4.06231715e-02   2.16376594e-03   2.44758076e-02
       0.00000000e+00   6.86198722e-03   1.36132952e-02   1.73886977e-02
       2.30213129e-02   4.44999980e-02   0.00000000e+00   1.40589569e-02
       0.00000000e+00   4.74213177e-02   2.65427674e-02   1.05684330e-01
       6.30552365e-03   2.86320444e-02   4.50079022e-03   7.76843152e-02
       2.88642900e-02   3.52207159e-02   2.01852506e-02   9.26784855e-04
       4.35733012e-02   1.84745904e-02   1.35102237e-02   2.69638287e-02
       1.88247064e-02   0.00000000e+00   2.03784688e-02   4.14981678e-02
       1.79538495e-02   1.12983577e-02   3.23765203e-02   0.00000000e+00
       3.99771399e-02   2.85164571e-03   2.18967289e-02   3.96111705e-02
       3.40096863e-02   1.72800650e-02   1.36861815e-02   0.00000000e+00
       1.19328203e-02   1.71726485e-02   0.00000000e+00   0.00000000e+00
       6.33251858e-03   4.64324980e-03   1.33084980e-03   9.89021626e-02
       3.52934995e-02   2.96267777e-02   1.73480268e-02   3.07545000e-02
       2.47891161e-02   3.32486832e-02   7.45403501e-03   1.46792267e-02
       0.00000000e+00   3.35642472e-02   8.78597450e-02   3.94517740e-02]

    References
    ----------
    .. [betweenness-wikipedia] http://en.wikipedia.org/wiki/Centrality#Betweenness_centrality
    .. [brandes-faster-2001] U. Brandes, "A faster algorithm for betweenness
       centrality", Journal of Mathematical Sociology, 2001, :doi:`10.1080/0022250X.2001.9990249`
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
    0.0813233725942

    References
    ----------
    .. [freeman-set-1977] Linton C. Freeman, "A Set of Measures of Centrality
       Based on Betweenness", Sociometry, Vol. 40, No. 1,  pp. 35-41, 1977,
       `http://www.jstor.org/stable/3033543 <http://www.jstor.org/stable/3033543>`_
    """

    return libgraph_tool_centrality.\
           get_central_point_dominance(g._Graph__graph,
                                       _prop("v", g, betweenness))


def eigenvector(g, weight=None, vprop=None, epsilon=1e-6, max_iter=None):
    r"""
    Calculate the eigenvector centrality of each vertex in the graph, as well as
    the largest eigenvalue.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    weights : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Edge property map with the edge weights.
    vprop : :class:`~graph_tool.PropertyMap`, optional (default: ``None``)
        Vertex property map where the values of eigenvector must be stored.
    epsilon : float, optional (default: ``1e-6``)
        Convergence condition. The iteration will stop if the total delta of all
        vertices are below this value.
    max_iter : int, optional (default: ``None``)
        If supplied, this will limit the total number of iterations.

    Returns
    -------
    eigenvalue : float
        The largest eigenvalue of the (weighted) adjacency matrix.
    eigenvector : :class:`~graph_tool.PropertyMap`
        A vertex property map containing the eigenvector values.

    See Also
    --------
    betweenness: betweenness centrality
    pagerank: PageRank centrality
    trust_transitivity: pervasive trust transitivity

    Notes
    -----

    The eigenvector centrality :math:`\mathbf{x}` is the eigenvector of the
    (weighted) adjacency matrix with the largest eigenvalue :math:`\lambda`,
    i.e. it is the solution of

    .. math::

        \mathbf{A}\mathbf{x} = \lambda\mathbf{x},


    where :math:`\mathbf{A}` is the (weighted) adjacency matrix and
    :math:`\lambda` is the largest eigenvalue.

    The algorithm uses the power method which has a topology-dependent complexity of
    :math:`O\left(N\times\frac{-\log\epsilon}{\log|\lambda_1/\lambda_2|}\right)`,
    where :math:`N` is the number of vertices, :math:`\epsilon` is the ``epsilon``
    parameter, and :math:`\lambda_1` and :math:`\lambda_2` are the largest and
    second largest eigenvalues of the (weighted) adjacency matrix, respectively.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    >>> from numpy.random import poisson, random, seed
    >>> seed(42)
    >>> g = gt.random_graph(100, lambda: (poisson(3), poisson(3)))
    >>> w = g.new_edge_property("double")
    >>> w.a = random(g.num_edges()) * 42
    >>> x = gt.eigenvector(g, w)
    >>> print x[0]
    0.0160851991895
    >>> print x[1].a
    [ 0.1376411   0.07207366  0.02727508  0.05805304  0.          0.10690994
      0.04315491  0.01040908  0.02300252  0.08874163  0.04968119  0.06718114
      0.05526028  0.20449371  0.02337425  0.07581173  0.19993899  0.14718912
      0.08464664  0.08474977  0.          0.04843894  0.          0.0089388
      0.16831573  0.00138653  0.11741616  0.          0.13455019  0.03642682
      0.06729803  0.06229526  0.08937098  0.05693976  0.0793375   0.04076743
      0.22176891  0.07717256  0.00518048  0.05722748  0.          0.00055799
      0.04541778  0.06420469  0.06189998  0.08011859  0.05377224  0.29979873
      0.01211309  0.15503588  0.02804072  0.1692873   0.01420732  0.02507
      0.02959899  0.02702304  0.1652933   0.01434992  0.1073001   0.04582697
      0.04618913  0.0220902   0.01421926  0.09891276  0.04522928  0.
      0.00236599  0.07686829  0.03243909  0.00346715  0.1954776   0.
      0.25583217  0.11710921  0.07804282  0.21188464  0.04800656  0.00321866
      0.0552824   0.11204116  0.11420818  0.24071304  0.15451676  0.
      0.00475456  0.10680434  0.17054333  0.18945499  0.15673649  0.03405238
      0.01653319  0.02563015  0.00186129  0.12061027  0.11449362  0.11114196
      0.06779788  0.00595725  0.09127559  0.02380386]

    References
    ----------

    .. [eigenvector-centrality] http://en.wikipedia.org/wiki/Centrality#Eigenvector_centrality
    .. [power-method] http://en.wikipedia.org/wiki/Power_iteration
    .. [langville-survey-2005] A. N. Langville, C. D. Meyer, "A Survey of
       Eigenvector Methods for Web Information Retrieval", SIAM Review, vol. 47,
       no. 1, pp. 135-161, 2005, :DOI:`10.1137/S0036144503424786`


    """

    if vprop == None:
        vprop = g.new_vertex_property("double")
    if max_iter is None:
        max_iter = 0
    ee = libgraph_tool_centrality.\
         get_eigenvector(g._Graph__graph, _prop("e", g, weight),
                         _prop("v", g, vprop), epsilon, max_iter)
    return ee, vprop


def eigentrust(g, trust_map, vprop=None, norm=False, epsilon=1e-6, max_iter=0,
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
    epsilon : float, optional (default: 1e-6)
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
    trust_transitivity: pervasive trust transitivity

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
    [ 0.02100449  0.01735932  0.00227182  0.00342703  0.          0.01739914
      0.00658874  0.00592764  0.00879695  0.01483758  0.00390145  0.00939709
      0.01038803  0.00896039  0.0080222   0.00583084  0.01510505  0.01106463
      0.02048866  0.0179936   0.02196625  0.00604554  0.          0.00038504
      0.01704679  0.00431482  0.00538866  0.          0.00163772  0.02009726
      0.00254747  0.00440903  0.02305541  0.01061566  0.00583414  0.01521545
      0.01894677  0.00941793  0.00259066  0.00454916  0.          0.
      0.00411855  0.01005776  0.029152    0.01500648  0.00797009  0.02057446
      0.          0.02100182  0.00519358  0.02503401  0.00368714  0.02176737
      0.00111934  0.02763714  0.00615445  0.00163793  0.01998869  0.
      0.00831816  0.00692008  0.00439715  0.01287125  0.00534507  0.
      0.00805071  0.02094972  0.00622514  0.00285397  0.01009464  0.
      0.00360911  0.00653993  0.00800227  0.01521205  0.02901848  0.01693622
      0.00323205  0.00748302  0.00443795  0.0076314   0.01147831  0.
      0.00129362  0.00173367  0.00188625  0.02110825  0.01349257  0.00956502
      0.00694694  0.01780551  0.00344632  0.02869166  0.00388418  0.0016279
      0.01691452  0.00783781  0.02795918  0.03327071]

    References
    ----------
    .. [kamvar-eigentrust-2003] S. D. Kamvar, M. T. Schlosser, H. Garcia-Molina
       "The eigentrust algorithm for reputation management in p2p networks",
       Proceedings of the 12th international conference on World Wide Web,
       Pages: 640 - 651, 2003, :doi:`10.1145/775152.775242`
    """

    if vprop == None:
        vprop = g.new_vertex_property("double")
    i = libgraph_tool_centrality.\
           get_eigentrust(g._Graph__graph, _prop("e", g, trust_map),
                          _prop("v", g, vprop), epsilon, max_iter)
    if norm:
        vprop.get_array()[:] /= sum(vprop.get_array())

    if ret_iter:
        return vprop, i
    else:
        return vprop


def trust_transitivity(g, trust_map, source=None, target=None, vprop=None):
    r"""
    Calculate the pervasive trust transitivity between chosen (or all) vertices
    in the graph.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    trust_map : :class:`~graph_tool.PropertyMap`
        Edge property map with the values of trust associated with each
        edge. The values must lie in the range [0,1].
    source : :class:`~graph_tool.Vertex` (optional, default: None)
        Source vertex. All trust values are computed relative to this vertex.
        If left unspecified, the trust values for all sources are computed.
    target : :class:`~graph_tool.Vertex` (optional, default: None)
        The only target for which the trust value will be calculated. If left
        unspecified, the trust values for all targets are computed.
    vprop : :class:`~graph_tool.PropertyMap` (optional, default: None)
        A vertex property map where the values of transitive trust must be
        stored.

    Returns
    -------
    trust_transitivity : :class:`~graph_tool.PropertyMap` or float
        A vertex vector property map containing, for each source vertex, a
        vector with the trust values for the other vertices. If only one of
        `source` or `target` is specified, this will be a single-valued vertex
        property map containing the trust vector from/to the source/target
        vertex to/from the rest of the network. If both `source` and `target`
        are specified, the result is a single float, with the corresponding
        trust value for the target.

    See Also
    --------
    eigentrust: eigentrust centrality
    betweenness: betweenness centrality
    pagerank: PageRank centrality

    Notes
    -----
    The pervasive trust transitivity between vertices i and j is defined as

    .. math::

        t_{ij} = \frac{\sum_m A_{m,j} w^2_{G\setminus\{j\}}(i\to m)c_{m,j}}
                 {\sum_m A_{m,j} w_{G\setminus\{j\}}(i\to m)}

    where :math:`A_{ij}` is the adjacency matrix, :math:`c_{ij}` is the direct
    trust from i to j, and :math:`w_G(i\to j)` is the weight of the path with
    maximum weight from i to j, computed as

    .. math::

       w_G(i\to j) = \prod_{e\in i\to j} c_e.

    The algorithm measures the transitive trust by finding the paths with
    maximum weight, using Dijkstra's algorithm, to all in-neighbours of a given
    target. This search needs to be performed repeatedly for every target, since
    it needs to be removed from the graph first. For each given source, the
    resulting complexity is therefore :math:`O(N^2\log N)` for all targets, and
    :math:`O(N\log N)` for a single target. For a given target, the complexity
    for obtaining the trust from all given sources is :math:`O(kN\log N)`, where
    :math:`k` is the in-degree of the target. Thus, the complexity for obtaining
    the complete trust matrix is :math:`O(EN\log N)`, where :math:`E` is the
    number of edges in the network.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    >>> from numpy.random import poisson, random, seed
    >>> seed(42)
    >>> g = gt.random_graph(100, lambda: (poisson(3), poisson(3)))
    >>> trust = g.new_edge_property("double")
    >>> trust.a = random(g.num_edges())
    >>> t = gt.trust_transitivity(g, trust, source=g.vertex(0))
    >>> print t.a
    [ 1.          0.09649648  0.01375374  0.09864347  0.          0.52668732
      0.02655169  0.05771735  0.25651251  0.13071344  0.1258206   0.13065921
      0.12051013  0.13754053  0.26727787  0.06951245  0.38774441  0.25343023
      0.21297027  0.59232433  0.10843174  0.02810649  0.          0.04000351
      0.13784095  0.06125175  0.04156937  0.          0.05771925  0.04967184
      0.11251086  0.25172931  0.1982562   0.28225643  0.05339001  0.10629504
      0.04440744  0.05815895  0.097983    0.03333347  0.          0.
      0.10845473  0.13751647  0.27567139  0.03946153  0.25063883  0.0755547   0.
      0.25167962  0.33205973  0.08237051  0.12983804  0.02587608  0.09694727
      0.16435599  0.09445501  0.07402817  0.06425702  0.          0.22420236
      0.11284837  0.05567628  0.0561254   0.36563496  0.          0.09358333
      0.06315609  0.3853858   0.01338133  0.08506159  0.          0.23226712
      0.0841518   0.07274848  0.17553984  0.14032908  0.15737553  0.13703351
      0.25035262  0.03570828  0.04341688  0.11955905  0.          0.01757771
      0.04990193  0.10457395  0.41668972  0.04546921  0.04404905  0.24922167
      0.09752267  0.03872946  0.26113888  0.04677363  0.03220735  0.03928181
      0.08696124  0.21697483  0.1388346 ]

    References
    ----------
    .. [richters-trust-2010] Oliver Richters, Tiago P. Peixoto, "Trust
       transitivity in social networks", :arXiv:`1012.1358`, 2010

    """

    if vprop == None:
        vprop = g.new_vertex_property("vector<double>")

    if target == None:
        target = -1
    else:
        target = g.vertex_index[target]

    if source == None:
        source = -1
    else:
        source = g.vertex_index[source]

    libgraph_tool_centrality.\
            get_trust_transitivity(g._Graph__graph, source, target,
                                   _prop("e", g, trust_map),
                                   _prop("v", g, vprop))
    if target != -1 or source != -1:
        vprop = ungroup_vector_property(vprop, [0])[0]
    if target != -1 and source != -1:
        return vprop.a[target]
    return vprop
