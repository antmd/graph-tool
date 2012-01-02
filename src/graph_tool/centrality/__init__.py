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
    [ 0.00865316  0.0054067   0.00406312  0.00426668  0.0015      0.00991696
      0.00550065  0.00936397  0.00347917  0.00731864  0.00689843  0.00286274
      0.00508731  0.01020047  0.00562247  0.00584915  0.02457086  0.00438568
      0.0057385   0.00621745  0.001755    0.0045073   0.0015      0.00225167
      0.00698342  0.00206302  0.01094466  0.001925    0.00710093  0.00519877
      0.00460646  0.00994648  0.01005248  0.00904629  0.00676221  0.00789208
      0.00933103  0.00301154  0.00264951  0.00842812  0.0015      0.00191034
      0.00594069  0.00884372  0.00453417  0.00388987  0.00317433  0.0086067
      0.00385394  0.00672702  0.00258411  0.01468262  0.00454     0.00381159
      0.00402607  0.00451133  0.00480966  0.00811557  0.00571949  0.00317433
      0.00856838  0.00280517  0.00280563  0.00906324  0.00614421  0.0015
      0.00292034  0.00479769  0.00552694  0.00604799  0.0115922   0.0015
      0.00676183  0.00695336  0.01023352  0.01737541  0.00451443  0.00197688
      0.00553866  0.00486233  0.0078653   0.00867599  0.01248092  0.0015
      0.00399605  0.00399605  0.00881571  0.00638008  0.01056944  0.00353724
      0.00249869  0.00684919  0.00241374  0.01061397  0.00673569  0.00590937
      0.01004638  0.00331612  0.00926359  0.00460809]

    Now with a personalization vector, and edge weights:

    >>> w = g.new_edge_property("double")
    >>> w.a = random(g.num_edges())
    >>> p = g.new_vertex_property("double")
    >>> p.a = random(g.num_vertices())
    >>> p.a /= p.a.sum()
    >>> pr = gt.pagerank(g, pers=p, weight=w)
    >>> print pr.a
    [ 0.00712999  0.00663336  0.00685722  0.00402663  0.00092715  0.01021926
      0.00269502  0.0073301   0.00449892  0.00582793  0.00580542  0.00275149
      0.00676363  0.01157972  0.00486918  0.00616345  0.02506695  0.00607967
      0.00553375  0.00359075  0.00293808  0.00362247  0.00250025  0.00186946
      0.00895516  0.00318147  0.01489786  0.00312436  0.0074751   0.0040342
      0.006254    0.00687051  0.0098073   0.01076278  0.00887077  0.00806759
      0.00969532  0.00252648  0.00278688  0.00972144  0.00148972  0.00215428
      0.00713602  0.00559849  0.00495517  0.00457118  0.00323767  0.01257406
      0.00120179  0.00514838  0.00130655  0.01724465  0.00343819  0.00420962
      0.00297617  0.00588287  0.00657206  0.00775082  0.00758217  0.00433776
      0.00576829  0.00464595  0.00307274  0.00585795  0.00745881  0.00238803
      0.00230431  0.00437046  0.00492464  0.00275414  0.01524646  0.00300867
      0.00816665  0.00548853  0.00874738  0.01871498  0.00216776  0.00245196
      0.00308878  0.00646323  0.01287978  0.00911384  0.01628604  0.0009367
      0.00222119  0.00864202  0.01199119  0.01126539  0.01086846  0.00309224
      0.0020319   0.00659422  0.00226965  0.0134399   0.01094141  0.00732916
      0.00489314  0.0030402   0.00783914  0.00278588]

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
    [ 0.04889806  0.07181892  0.0256799   0.02885791  0.          0.05060927
      0.04490836  0.03763462  0.02033383  0.03163202  0.02641248  0.03171598
      0.03771112  0.02194663  0.0374907   0.01072567  0.          0.03079281
      0.05409258  0.00163434  0.00051978  0.01045902  0.          0.00796784
      0.0494527   0.00647576  0.03708252  0.00304503  0.0663657   0.03903257
      0.03305169  0.          0.07787098  0.03938866  0.08577116  0.020183
      0.06024004  0.01004935  0.0443127   0.06397736  0.          0.00363548
      0.01742486  0.03216543  0.01918144  0.02059159  0.          0.01476213
      0.          0.0466751   0.01072612  0.10288046  0.00563973  0.03850413
      0.00629595  0.01292137  0.0537963   0.04454985  0.01227018  0.00729488
      0.02092959  0.02308238  0.00712703  0.02193975  0.03823342  0.
      0.00995364  0.04023839  0.0312708   0.0111312   0.00228516  0.
      0.09659583  0.01327402  0.05792071  0.08606828  0.0143541   0.00221604
      0.02144698  0.          0.04023879  0.00715758  0.          0.
      0.02348452  0.00760922  0.01486521  0.08132792  0.0382674   0.03078318
      0.00430209  0.01772787  0.02280666  0.0373011   0.03077511  0.02871265
      0.          0.01044655  0.04415432  0.04447525]

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
    0.0766473408634

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
    weight : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
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
    vprop : :class:`~graph_tool.PropertyMap`, optional (default: ``None``)
        Vertex property map where the values of eigentrust must be stored.
    norm : bool, optional (default:  ``False``)
        Norm eigentrust values so that the total sum equals 1.
    epsilon : float, optional (default: ``1e-6``)
        Convergence condition. The iteration will stop if the total delta of all
        vertices are below this value.
    max_iter : int, optional (default: ``None``)
        If supplied, this will limit the total number of iterations.
    ret_iter : bool, optional (default: ``False``)
        If true, the total number of iterations is also returned.

    Returns
    -------
    eigentrust : :class:`~graph_tool.PropertyMap`
        A vertex property map containing the eigentrust values.

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
    >>> trust.a = random(g.num_edges())*42
    >>> t = gt.eigentrust(g, trust, norm=True)
    >>> print t.a
    [  1.12095562e-02   3.97280231e-03   1.31675503e-02   9.61282478e-03
       0.00000000e+00   1.73295741e-02   3.53395497e-03   1.06203582e-02
       1.36906165e-03   8.64587777e-03   1.12049516e-02   3.18891993e-03
       9.28265221e-03   2.25294315e-02   3.24795656e-03   9.16555333e-03
       5.68412465e-02   6.79686311e-03   6.37474649e-03   6.04696712e-03
       0.00000000e+00   8.51131034e-03   0.00000000e+00   1.09336777e-03
       1.49885187e-02   1.09327367e-04   3.73928902e-02   0.00000000e+00
       1.74638522e-02   8.21101864e-03   5.79876899e-03   1.34905262e-02
       1.71525132e-02   2.25425503e-02   1.04184903e-02   1.05537922e-02
       1.34096247e-02   2.82760533e-03   4.31713918e-04   7.39114668e-03
       0.00000000e+00   2.21328121e-05   8.79050007e-03   7.08148889e-03
       5.88651144e-03   7.45401425e-03   5.66098580e-03   2.80738199e-02
       2.41472197e-03   1.00673881e-02   2.29910658e-03   3.23790630e-02
       3.02136064e-03   2.25030440e-03   3.53325357e-03   6.90672383e-03
       1.01692058e-02   1.03783022e-02   1.22476413e-02   4.82453065e-03
       1.15878890e-02   3.41943633e-03   1.57958469e-03   6.56648121e-03
       1.28152141e-02   0.00000000e+00   1.29192164e-03   9.35867476e-03
       3.89329603e-03   1.78002682e-03   2.81987911e-02   0.00000000e+00
       1.74943514e-02   6.24079508e-03   1.57572103e-02   3.77119257e-02
       4.78552984e-03   3.30463136e-04   5.60118687e-03   5.75656186e-03
       2.65412905e-02   1.59663210e-02   2.88844192e-02   0.00000000e+00
       7.87754853e-04   1.76957899e-02   3.19907905e-02   1.94650690e-02
       1.32052233e-02   3.57577093e-03   7.09968545e-04   8.70787481e-03
       1.24901391e-04   2.61215462e-02   2.25923034e-02   1.10928239e-02
       9.39210737e-03   5.61073138e-04   1.59987179e-02   3.02799309e-03]

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
    [  1.00000000e+00   9.59916062e-02   4.27717883e-02   7.70755875e-02
       0.00000000e+00   2.04476926e-01   5.55315822e-02   2.82854665e-02
       5.08479257e-02   1.68128402e-01   3.28567434e-02   7.39525583e-02
       1.34463196e-01   8.83740756e-02   1.79990535e-01   7.08809615e-02
       6.37757645e-02   7.24187957e-02   4.83082241e-02   9.90676983e-02
       0.00000000e+00   6.50497060e-02   0.00000000e+00   1.77344948e-02
       1.08677897e-01   1.00958718e-03   4.49524961e-02   0.00000000e+00
       1.64902280e-01   4.31492976e-02   2.19446085e-01   3.00890381e-02
       6.86750847e-02   2.72460575e-02   3.57314594e-02   4.87776483e-02
       4.11748930e-01   7.91396467e-02   2.54835127e-03   3.01711432e-01
       0.00000000e+00   4.14406224e-04   4.24794624e-02   9.14096554e-02
       4.17528677e-01   3.79112573e-02   1.16489950e-01   5.18112902e-02
       8.49111259e-03   5.26399996e-02   2.45690139e-02   7.51435125e-02
       5.62381854e-02   2.90115777e-02   2.72543383e-02   1.46877163e-01
       7.81446822e-02   1.24417763e-02   1.01337976e-01   9.92776442e-02
       3.14622176e-02   1.20097319e-01   3.30335980e-02   4.61757040e-02
       1.01085599e-01   0.00000000e+00   4.44660446e-03   6.31066845e-02
       1.94702084e-02   8.45343379e-04   4.82190327e-02   0.00000000e+00
       6.60346087e-02   7.44581695e-02   6.19535229e-02   1.82072422e-01
       1.45366611e-02   2.59020075e-02   2.52208295e-02   6.80519730e-02
       6.74671969e-02   1.14198914e-01   5.12493343e-02   0.00000000e+00
       6.33427008e-03   1.42290348e-01   6.90459437e-02   1.00565411e-01
       5.88966867e-02   3.28157280e-02   2.80046903e-02   2.41520032e-01
       8.45879329e-04   6.76633672e-02   6.05080467e-02   9.12575826e-02
       1.97789973e-02   6.40885493e-02   4.80934526e-02   1.28787181e-02]

    References
    ----------
    .. [richters-trust-2010] Oliver Richters and Tiago P. Peixoto, "Trust
       Transitivity in Social Networks," PLoS ONE 6, no. 4:
       e1838 (2011), :doi:`10.1371/journal.pone.0018384`

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
