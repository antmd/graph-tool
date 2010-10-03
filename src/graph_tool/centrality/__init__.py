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
   trust_transitivity

Contents
++++++++
"""

from .. dl_import import dl_import
dl_import("import libgraph_tool_centrality")

from .. core import _prop, ungroup_vector_property
import sys
import numpy

__all__ = ["pagerank", "betweenness", "central_point_dominance", "eigentrust",
           "trust_transitivity"]


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
    trust_transitivity: pervasive trust transitivity

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
    [ 0.87011681  1.73449398  0.47587866  0.4534494   0.2         1.26596887
      0.60964865  0.68064477  0.8137542   0.86269096  0.51833002  0.49194604
      0.74875795  0.52831993  0.601438    0.63921165  1.32489495  0.68360746
      1.02608206  0.90903761  1.1026286   0.56290713  0.2         0.30840086
      0.90726785  0.35583967  0.95582862  0.232       0.41090313  0.88734742
      0.47424296  0.66138242  1.26313184  0.7459428   0.84110051  0.9497316
      1.0589998   0.94412292  0.26433617  0.86197354  0.2         0.25333333
      0.65974242  0.69889305  1.02798531  0.77618244  0.57905885  1.12828577
      0.232       1.18366748  0.38929224  1.72424164  0.47966878  1.0931673
      0.45937603  1.09479766  0.80274459  0.44782081  1.04618114  0.25333333
      0.82295953  0.40210109  0.72779393  0.75075946  0.41742276  0.2
      0.8984279   0.92941713  0.69682427  0.69340983  1.02679348  0.2
      0.67750539  0.85622403  0.77232588  1.09093307  1.14410169  0.59413937
      0.54456339  0.64371752  0.40275133  0.72976606  1.40446885  0.2
      0.31831299  0.3734494   0.2562224   1.05807688  1.02419007  0.82747632
      0.49646186  0.72960178  0.48621114  1.42147072  0.65622314  0.31664379
      1.55387576  0.58439879  2.03922765  1.47802266]

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
    0.0813233725942

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
    [ 0.01610395  0.03518828  0.00387335  0.00506519  0.          0.02120586
      0.00328345  0.00514034  0.00361398  0.01331587  0.00626757  0.00788882
      0.01599836  0.00607798  0.00879484  0.01028104  0.01742029  0.00522399
      0.0206618   0.0098984   0.00918508  0.01344131  0.          0.00047679
      0.01760032  0.00078869  0.01045936  0.          0.00387405  0.01761267
      0.00730843  0.00514523  0.01708638  0.0084908   0.01237811  0.01401104
      0.0209564   0.0132232   0.00031255  0.01400855  0.          0.          0.0077233
      0.00479587  0.01646928  0.01499744  0.01901516  0.00843277  0.
      0.01764526  0.00243523  0.01726375  0.01272935  0.0163525   0.00382533
      0.02037745  0.00758792  0.00350063  0.01303079  0.          0.02086308
      0.00062028  0.00841231  0.00983605  0.00327547  0.          0.01016667
      0.0170241   0.00782474  0.00516862  0.02394048  0.          0.00747778
      0.00792131  0.01495136  0.01513948  0.02287957  0.00788276  0.0053207
      0.00145811  0.00183203  0.0033493   0.01627589  0.          0.00476343
      0.00937439  0.00200381  0.01400712  0.02135004  0.00549685  0.00230923
      0.01426992  0.01083921  0.03439618  0.00514281  0.00114438  0.02259093
      0.00672266  0.02753108  0.01859351]

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
    source : Vertex (optional, default: None)
        Source vertex. All trust values are computed relative to this vertex.
        If left unspecified, the trust values for all sources are computed.
    target : Vertex (optional, default: None)
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
    The absolute trust between vertices i and j is defined as

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
    [ 1.          0.15271582  0.07130332  0.10597708  0.          0.58940763
      0.04233924  0.03619048  0.04137002  0.05926363  0.06584407  0.06315985
      0.22301815  0.02671845  0.10566551  0.08018763  0.57668762  0.08440303
      0.17612948  0.37579015  0.0415804   0.19919108  0.          0.0141547
      0.14901031  0.00910391  0.02680543  0.          0.0887711   0.0296914
      0.09800672  0.06421615  0.16420105  0.10226839  0.08667606  0.07944174
      0.17174637  0.10932321  0.0137295   0.09342906  0.          0.
      0.11065065  0.03725047  0.23554212  0.10971862  0.54564134  0.0462946   0.
      0.24820041  0.15281463  0.09449931  0.22419781  0.03108608  0.10964166
      0.08642532  0.03495468  0.05656444  0.04045297  0.          0.13789871
      0.0197414   0.05512572  0.08297112  0.21448002  0.          0.08649514
      0.0718887   0.16546776  0.04108292  0.11710843  0.          0.12518596
      0.04797708  0.02275816  0.10413969  0.1294644   0.08656727  0.28371423
      0.1036658   0.01575087  0.02023104  0.067158    0.          0.03241519
      0.19613692  0.05684533  0.29652909  0.03038526  0.02423028  0.01695595
      0.0759531   0.17360708  0.51113999  0.03714076  0.03167552  0.04359062
      0.0267188   0.47605313  0.06471942]
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
