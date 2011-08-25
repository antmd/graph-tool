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
    [ 0.00867754  0.00729246  0.00363279  0.00668265  0.0015      0.00859964
      0.00449637  0.00961946  0.01295288  0.00882362  0.00719256  0.00280697
      0.00518114  0.01047904  0.00569656  0.00519058  0.00759745  0.00700835
      0.00870244  0.00522561  0.00233159  0.00236035  0.0015      0.00255374
      0.00872139  0.00227483  0.00686341  0.001755    0.00488567  0.01045994
      0.00393206  0.00988283  0.01376133  0.00721883  0.01429166  0.00752748
      0.01846797  0.00674401  0.00412138  0.00842639  0.0015      0.00233159
      0.00306271  0.01902149  0.0099247   0.00428981  0.00215072  0.01123842
      0.00236035  0.00768803  0.00463719  0.01130437  0.00392423  0.00491263
      0.00899519  0.00680983  0.0091988   0.00459334  0.00809094  0.00881614
      0.01381946  0.00489171  0.00425249  0.01957383  0.00708763  0.0015
      0.00418613  0.00607306  0.01287535  0.00639268  0.01578391  0.0015
      0.01541987  0.00860721  0.01378758  0.0173314   0.00775072  0.00247939
      0.00524088  0.00686587  0.00436895  0.00755964  0.00708251  0.0015
      0.00226438  0.00184085  0.00555171  0.01159494  0.01297596  0.00460887
      0.00406717  0.00578091  0.00548516  0.01197071  0.00674202  0.01011666
      0.01072786  0.00646937  0.01430012  0.01483996]

    Now with a personalization vector, and edge weights:

    >>> w = g.new_edge_property("double")
    >>> w.a = random(g.num_edges())
    >>> p = g.new_vertex_property("double")
    >>> p.a = random(g.num_vertices())
    >>> p.a /= p.a.sum()
    >>> pr = gt.pagerank(g, pers=p, weight=w)
    >>> print pr.a
    [ 0.00761942  0.00761689  0.00418938  0.00947758  0.00092715  0.00349991
      0.00811226  0.00448968  0.01209889  0.01384828  0.00600036  0.00221745
      0.00432908  0.01036427  0.00536132  0.00692364  0.00575216  0.00750936
      0.00924536  0.00461255  0.00422277  0.00055639  0.00250025  0.00289125
      0.00925617  0.003356    0.00642017  0.00298276  0.00571097  0.01115541
      0.00452616  0.01670105  0.01788592  0.00580217  0.01350007  0.00837655
      0.01535733  0.00497981  0.00436008  0.01324374  0.00148972  0.00287379
      0.00408663  0.02785282  0.00790422  0.00491795  0.00070143  0.00789247
      0.00033551  0.00777089  0.00278393  0.00801468  0.00452296  0.00378295
      0.00642244  0.00698618  0.01069855  0.0019177   0.00742151  0.00872767
      0.01868187  0.00442359  0.00593616  0.01517386  0.00712472  0.00238803
      0.00468324  0.0024983   0.011788    0.00577489  0.01242015  0.00300867
      0.01390361  0.00796192  0.00822753  0.01062897  0.00815637  0.00332914
      0.00911336  0.00915715  0.00945334  0.00880299  0.00758402  0.0009367
      0.00378413  0.00174124  0.00283594  0.00929262  0.01090867  0.00460206
      0.00341061  0.00699703  0.00232131  0.01244958  0.00731098  0.01288061
      0.00820259  0.00430521  0.01633379  0.0119308 ]

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
    [ 0.02412512  0.08233748  0.01789612  0.03997773  0.          0.03476439
      0.03490215  0.02812729  0.05385124  0.01614861  0.01894254  0.03552189
      0.01648793  0.02743878  0.02243743  0.0052126   0.          0.02648145
      0.05045875  0.01670867  0.00027069  0.00235053  0.          0.00424986
      0.02153982  0.01635659  0.03177692  0.00152088  0.02099129  0.05311383
      0.00802715  0.          0.06129706  0.03148129  0.10265395  0.02762436
      0.05323276  0.028209    0.01641328  0.03547918  0.          0.00998142
      0.01043084  0.06810444  0.01435047  0.02884138  0.          0.02336079
      0.          0.09098673  0.02911358  0.05909676  0.01314448  0.0304931
      0.01315283  0.03795536  0.02756845  0.020655    0.00268182  0.0151298
      0.05597887  0.04095171  0.00853146  0.06176456  0.03165429  0.          0.0114488
      0.05435076  0.04046437  0.0122299   0.03417369  0.          0.05979601
      0.01502672  0.05908044  0.09567821  0.01206665  0.01259951  0.02381255
      0.          0.01840359  0.00708731  0.          0.          0.01274922
      0.00385357  0.          0.11470908  0.04903743  0.03336991  0.01270614
      0.00448173  0.02163365  0.06685948  0.01032924  0.02400925  0.          0.0300603
      0.03220004  0.03932736]

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
    0.0888588565184

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
    0.0144966901652
    >>> print x[1].a
    [ 0.06930562  0.15259089  0.02186899  0.16605389  0.          0.03329417
      0.15389839  0.05956241  0.02364241  0.14555449  0.12223866  0.04284837
      0.02930074  0.12362235  0.08914056  0.06613553  0.04127041  0.06894616
      0.11768146  0.10576438  0.03184658  0.01097733  0.          0.00244834
      0.12373674  0.04521477  0.06729261  0.          0.0119752   0.09613809
      0.05371652  0.0760725   0.26342747  0.06901843  0.12261332  0.07428021
      0.16074381  0.04316101  0.01382325  0.07976272  0.          0.02358784
      0.04087066  0.22404933  0.05037626  0.07311321  0.00118498  0.10162138
      0.00695398  0.15475802  0.03601454  0.12771193  0.00204537  0.00850904
      0.06338004  0.11785637  0.171928    0.00577511  0.08163299  0.02358024
      0.25112631  0.05618482  0.00868661  0.15998601  0.11513557  0.
      0.07307991  0.06984483  0.1043859   0.15786829  0.12612902  0.
      0.11247494  0.06210605  0.07547538  0.09774456  0.01932     0.01132741
      0.07459586  0.14437367  0.06477981  0.13211211  0.04886609  0.
      0.09043424  0.00644994  0.0112379   0.10472704  0.2077468   0.04146317
      0.04745921  0.07442587  0.00757032  0.18575827  0.09975438  0.11682436
      0.22233397  0.08072367  0.13527766  0.16050936]

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
    [ 0.01047704  0.00916079  0.00329728  0.02183509  0.          0.00570515
      0.01795154  0.00745325  0.01556352  0.02103068  0.00766733  0.0017898
      0.00477604  0.02075469  0.00596795  0.00777273  0.00486571  0.00922539
      0.00947171  0.00808096  0.00159634  0.00039577  0.          0.00170102
      0.01221233  0.00105577  0.00499953  0.          0.00104692  0.02124238
      0.00460677  0.02377163  0.03682221  0.00788995  0.01588396  0.01090949
      0.01738658  0.00629269  0.00252249  0.01494837  0.          0.00118236
      0.00086828  0.04746846  0.01541295  0.00777287  0.00014366  0.01079173
      0.00025071  0.0110701   0.0041122   0.00931269  0.0003396   0.00085525
      0.00715489  0.01178172  0.01968206  0.00086696  0.01155932  0.01291423
      0.03650986  0.00409654  0.00635913  0.02706258  0.01329375  0.          0.004781
      0.00393831  0.0159743   0.01175119  0.01924371  0.          0.01624875
      0.01213196  0.01327894  0.02122698  0.01120801  0.00261826  0.01197119
      0.01071532  0.01065604  0.00891656  0.00394593  0.          0.00519856
      0.00013063  0.00206336  0.01471883  0.02029324  0.00793793  0.00560927
      0.00944371  0.00205155  0.01610597  0.009506    0.0211548   0.01688137
      0.00731456  0.0175582   0.02643281]

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
    [ 1.          0.05067451  0.0289009   0.08405839  0.          0.03545782
      0.07538031  0.22379499  0.05854478  0.07125532  0.05587198  0.02518448
      0.02466963  0.15362921  0.06018343  0.05567524  0.27358358  0.05275904
      0.03729045  0.09492006  0.05238415  0.00737096  0.          0.00510552
      0.08817629  0.03785795  0.02446122  0.          0.00738919  0.05119031
      0.34527065  0.03114924  0.07197484  0.47425602  0.03472084  0.04662672
      0.05560849  0.02364994  0.01215894  0.31788601  0.          0.03879942
      0.01857147  0.06889803  0.04377159  0.08730759  0.00246369  0.39098931
      0.0046694   0.04502054  0.02463968  0.02645761  0.00707666  0.00472783
      0.02984028  0.1148617   0.0572029   0.00675625  0.02800698  0.12556046
      0.08124247  0.05626964  0.06027068  0.04269226  0.07518783  0.
      0.02654427  0.013547    0.35126413  0.05338214  0.02683736  0.
      0.03967174  0.15764245  0.03434911  0.06712879  0.01081633  0.01438191
      0.22453002  0.13426137  0.11622003  0.18696878  0.02270863  0.
      0.09096154  0.00622639  0.01807092  0.06891202  0.05429632  0.03123112
      0.04239384  0.06347529  0.00980117  0.08103904  0.06022526  0.09774797
      0.0442361   0.03511221  0.09291923  0.10395899]

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
