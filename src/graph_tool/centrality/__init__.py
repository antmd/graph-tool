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
    >>> print pr.get_array()
    [ 0.99988081  0.39997616  0.80428057  0.43237369  0.2         0.75830329
      0.41447482  1.56621542  0.30841665  0.86432715  0.79374139  0.54573086
      0.89372179  0.93590145  0.25159724  1.12033843  0.2         0.98486039
      0.28191404  0.88133806  0.31166878  1.73878838  0.6903469   0.94100349
      0.25159724  0.32248278  1.03788472  0.58022932  0.35009064  0.94542317
      0.85751934  0.69608227  1.11373543  1.13477707  0.2         0.71559888
      0.30461189  0.2         1.02871995  1.14657561  0.2         0.25031945
      0.51841423  0.44709022  0.75239816  0.76551737  0.25638281  1.51657252
      0.30841665  0.59707408  0.34179258  1.0590272   2.16427996  0.51196274
      1.2264604   1.71578696  0.85838961  0.41931136  0.96797602  0.61882367
      1.07826603  0.2984934   1.1305187   0.75006564  0.48066231  1.61759314
      0.73870051  1.08374044  0.38258693  0.98112013  0.2         0.25590818
      1.17500568  1.2288973   0.29613246  1.45937444  0.39997616  1.18311783
      0.67063807  0.39229458  0.72314004  0.88473325  0.32859279  0.40656244
      0.51754349  0.5315028   0.55196274  0.2335463   1.56357203  0.91464458
      0.46999727  1.06779933  0.4852867   0.48933035  0.58997931  0.52883683
      0.79385874  0.59244805  0.99896399  1.0470592 ]

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
    >>> print vb.get_array()
    [ 0.06129648  0.02004734  0.04305659  0.01071136  0.          0.0252028
      0.00679622  0.06981881  0.00541371  0.02462107  0.05328111  0.0107051
      0.05981227  0.          0.01315561  0.00131498  0.          0.01883264
      0.01663386  0.03195175  0.01942617  0.13693745  0.01378875  0.00962001
      0.01325009  0.04685362  0.03839758  0.03395201  0.02160984  0.01727593
      0.0478231   0.          0.03826993  0.05124999  0.          0.
      0.00705917  0.          0.02190356  0.04505211  0.          0.00676419
      0.00110802  0.00169839  0.08733666  0.10546473  0.          0.12058932
      0.          0.00907921  0.02182859  0.08865455  0.          0.0418017
      0.03500162  0.07492683  0.03856307  0.04300598  0.02173347  0.00488363
      0.03739852  0.01113193  0.04386369  0.02994719  0.03383728  0.
      0.09230395  0.05449223  0.02507715  0.04944675  0.          0.00215935
      0.04371057  0.01749238  0.00104315  0.04688928  0.00444627  0.0178016
      0.01358585  0.02193068  0.03184527  0.05640358  0.00214389  0.03922583
      0.02195544  0.02613584  0.02246488  0.00066481  0.0755375   0.03142692
      0.04533332  0.03188087  0.04227853  0.03926328  0.00810412  0.02888085
      0.0455241   0.01373183  0.07029039  0.04382892]

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
    0.108411171667

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
        edge. The values must not lie in the range [0,1].
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
    [  1.78295032e-02   1.10159977e-03   8.27504534e-03   3.34579667e-03
       0.00000000e+00   9.28795883e-03   7.56225537e-03   2.03772288e-02
       6.87447577e-04   8.87085111e-03   2.84707349e-03   2.55095571e-03
       7.65302351e-03   5.06044724e-02   3.98617107e-04   1.02897822e-02
       0.00000000e+00   6.76980749e-03   6.91342330e-04   1.13998018e-02
       1.91846222e-03   3.74940757e-02   8.65907932e-03   5.76596060e-03
       1.11786939e-05   8.20855949e-04   9.45056085e-03   1.76099276e-02
       2.67746802e-03   1.03182164e-02   1.80748361e-02   8.49781556e-03
       7.89442825e-03   1.11838761e-02   0.00000000e+00   4.37095317e-03
       2.50451228e-05   0.00000000e+00   6.04054677e-03   1.51361293e-02
       0.00000000e+00   1.62557422e-04   1.02859153e-03   3.38079641e-03
       3.06115271e-03   2.96226918e-03   7.40021010e-05   1.64096932e-02
       1.12026631e-03   3.33521569e-03   1.77214999e-03   6.62472745e-03
       3.17014482e-02   1.93793538e-03   5.24056364e-02   4.04200200e-02
       2.96053927e-02   2.06294202e-03   2.93045979e-02   1.87688605e-03
       1.13962350e-02   6.94033709e-03   1.57347756e-02   3.97987237e-03
       1.15994824e-03   1.81252731e-02   2.06848985e-02   3.73314296e-03
       1.27163202e-03   1.08081901e-02   0.00000000e+00   2.25590063e-04
       8.55970439e-03   4.15387826e-02   8.61792076e-05   6.48435253e-02
       5.61799591e-03   4.69096686e-02   4.24627753e-03   9.16721227e-04
       4.86522362e-03   4.42735866e-03   5.50595265e-04   3.12087221e-03
       8.75442087e-03   4.25588041e-03   2.91851609e-03   1.80331544e-06
       2.89281502e-02   1.75099401e-02   1.14704807e-02   3.30940821e-02
       2.84005465e-03   4.92435108e-03   4.34713976e-03   2.72336599e-03
       9.37679329e-03   8.64912360e-03   3.96113432e-03   1.07637051e-02]

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
        A vertex which is used the as the source for gathering trust values.
    target : Vertex (optional, default: None)
        A vertex which is used the as the only target for which the trust value
        will be calculated. If left unspecified, the trust values for all
        targets are computed.
    vprop : :class:`~graph_tool.PropertyMap`, optional (default: None)
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
    [ 0.05927703  0.06133836  0.          0.05630559  0.          0.03317174
      0.03488483  0.15920558  0.16940159  0.09716039  0.1485169   0.0120287
      0.03787312  0.37284274  0.00646336  0.0084941   0.0379645   0.07997339
      0.10733769  0.10053845  0.00283938  0.05224064  0.          0.16523684
      0.0393326   0.25853808  0.14682555  0.03254906  0.12124144  0.0118341
      0.18110839  0.18513216  0.05031324  0.04484457  0.17197674  0.08569659
      0.17523371  0.22435776  0.33916191  0.07980329  0.          0.
      0.09750183  0.09811054  0.14574289  0.0085499   0.34593499  0.03151408
      0.083739    0.05409947  0.09161205  0.19921201  0.10647812  0.21597253
      0.06266044  0.8738786   0.11239455  0.09493216  0.19073287  0.11968616
      0.13409125  0.00626821  0.05857625  0.05917779  0.05673643  0.
      0.02682173  0.00355514  0.17475858  0.15113517  0.13247358  0.
      0.04003866  0.00997401  0.11126411  0.07400706  0.11247583  0.10125886
      0.16028191  0.04300862  0.03259707  0.0225482   0.05538721  0.
      0.06715919  0.0701153   0.02999368  0.04675702  0.06310919  0.01722603
      0.18455906  0.08034113  0.00376382  0.10041304  0.3437539   0.10530238
      0.11654855  0.09495419  0.05317485  0.10727767]
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

