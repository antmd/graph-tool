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
           "trust_transitivity"]


def pagerank(g, damping=0.8, prop=None, epsilon=1e-6, max_iter=None,
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

        PR(v) = \frac{1-d}{N} + d \sum_{w \in \Gamma^{-}(v)}
                \frac{PR (w)}{d^{+}(w)}

    where :math:`\Gamma^{-}(v)` are the in-neighbours of v, :math:`d^{+}(w)` is
    the out-degree of w, and d is a damping factor.

    The implemented algorithm progressively iterates the above condition, until
    it no longer changes, according to the parameter epsilon. It has a
    topology-dependent running time.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    >>> from numpy.random import poisson, seed
    >>> seed(42)
    >>> g = gt.random_graph(100, lambda: (poisson(3), poisson(3)))
    >>> pr = gt.pagerank(g)
    >>> print pr.a
    [ 0.0087012   0.01734503  0.0047588   0.00453451  0.002       0.01265973
      0.0060965   0.00680647  0.00813758  0.00862694  0.00518331  0.00491948
      0.00748761  0.00528322  0.00601439  0.00639214  0.013249    0.0068361
      0.01026087  0.00909041  0.01102634  0.0056291   0.002       0.00308401
      0.00907272  0.0035584   0.00955833  0.00232     0.00410904  0.00887352
      0.00474244  0.00661384  0.01263138  0.00745946  0.00841104  0.00949735
      0.01059004  0.00944125  0.00264336  0.00861976  0.002       0.00253333
      0.00659745  0.00698895  0.01027991  0.00776186  0.00579061  0.01128291
      0.00232     0.01183673  0.00389293  0.01724249  0.0047967   0.01093172
      0.00459377  0.01094803  0.00802747  0.00447822  0.01046185  0.00253333
      0.00822962  0.00402102  0.00727797  0.00750763  0.00417424  0.002
      0.00898431  0.00929422  0.00696827  0.00693413  0.01026798  0.002
      0.00677507  0.00856227  0.00772329  0.01090938  0.01144107  0.00594142
      0.00544564  0.0064372   0.00402752  0.00729768  0.01404475  0.002
      0.00318314  0.00373451  0.00256223  0.01058081  0.01024193  0.0082748
      0.00496463  0.00729605  0.00486213  0.01421478  0.00656225  0.00316644
      0.01553884  0.005844    0.02039237  0.01478031]

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
            get_pagerank(g._Graph__graph, _prop("v", g, prop), damping, epsilon,
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
