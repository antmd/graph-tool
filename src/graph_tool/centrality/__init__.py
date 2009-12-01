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
    >>> print pr.a
    [ 0.89482844  1.37847566  0.24        1.30716676  0.2         0.70397009
      0.40205781  0.74783725  1.37167015  0.66836587  0.5868133   0.47968714
      1.52225854  1.07388611  0.76316432  0.39214247  0.9302883   0.86455762
      0.77546264  1.87740317  0.25482139  0.29902553  0.2         0.24756383
      0.97205301  0.29727392  1.34742309  0.30905457  0.55032542  0.56654712
      0.40895463  0.77928729  0.73227413  0.59911926  1.39946277  0.72793699
      2.27008393  0.88929335  0.48636962  0.73070609  0.2         0.232
      0.96857512  2.97683022  0.58581032  0.80217847  0.37896569  0.93866821
      0.27337672  0.98201842  0.48551839  1.22651796  0.73263045  0.43013228
      1.00971133  0.72075953  0.66715456  0.58705749  0.74286661  0.37785867
      1.8475279   0.26432925  0.33994628  0.97319326  0.78104447  0.2
      0.33333761  0.51756267  0.47811583  0.85905246  1.46428623  0.2
      1.70687671  1.0107342   0.94504737  1.29858046  2.19707395  0.55931282
      0.85129509  1.09493368  1.22168331  0.64108136  0.70690188  0.2
      0.31736266  0.42372513  0.79429328  1.44749664  1.20741669  0.65763236
      0.40895463  0.62628812  0.32671006  0.85626447  0.59925496  0.3399879
      0.81215046  0.71506902  2.25678844  1.04882679]

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
    [ 0.03047981  0.07396685  0.00270882  0.044637    0.          0.03259048
      0.0243547   0.04265909  0.06274696  0.01778475  0.03502657  0.02692273
      0.05170277  0.05522454  0.02303023  0.0038858   0.          0.04852871
      0.02398655  0.00232365  0.          0.01064643  0.          0.01105872
      0.03564021  0.0222059   0.05170383  0.00140447  0.03935299  0.02644813
      0.01831885  0.          0.0453981   0.04552396  0.1242787   0.04983878
      0.07248363  0.04676976  0.03481327  0.04473583  0.          0.0027417
      0.01061048  0.0470108   0.01059109  0.05290495  0.          0.02541583
      0.          0.04012033  0.02616307  0.09056515  0.01640322  0.01599007
      0.02784563  0.05008998  0.03788222  0.03028745  0.01097982  0.00178571
      0.05804645  0.01015181  0.0061582   0.0255485   0.05504439  0.
      0.00179516  0.03367643  0.00304982  0.02333254  0.00843039  0.
      0.05947385  0.01936996  0.0521946   0.04928937  0.03955121  0.01360865
      0.02942447  0.          0.05149102  0.01054765  0.          0.
      0.00537915  0.01251828  0.01097982  0.06667564  0.04090169  0.02161779
      0.02941671  0.01793679  0.02360528  0.02638257  0.0062989   0.00946123
      0.          0.02255701  0.05081734  0.04846652]

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
    0.0980212339559

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
    [  1.04935746e-02   2.82745068e-02   0.00000000e+00   1.81121002e-02
       0.00000000e+00   3.70898521e-03   1.00108703e-03   1.29620638e-02
       1.71874047e-02   7.07523828e-03   8.29873222e-03   1.79259666e-03
       4.08925756e-02   1.55855653e-02   2.92256968e-03   1.71520782e-03
       5.04335865e-03   1.25678184e-02   1.92903241e-02   2.46642649e-02
       1.76431290e-04   1.85066489e-04   0.00000000e+00   4.52686439e-04
       7.13943855e-03   2.36002975e-03   1.44366165e-02   4.39632543e-04
       7.50316671e-03   8.13521884e-03   3.98083843e-03   1.04883920e-02
       7.42099689e-03   2.46651355e-03   2.08148781e-02   8.02104873e-03
       2.59366573e-02   2.11125347e-02   7.45781416e-03   6.62338254e-03
       0.00000000e+00   0.00000000e+00   1.72521147e-02   4.74346499e-02
       8.10593668e-03   2.27229702e-02   2.21525586e-03   6.24223052e-03
       2.59753300e-03   9.15181124e-03   3.67310718e-03   1.18998211e-02
       1.66177496e-02   6.44748287e-03   8.01978992e-03   1.48621102e-02
       6.65606246e-03   3.39887550e-03   1.20188240e-02   3.51012614e-03
       2.79661104e-02   7.90103914e-05   1.18015521e-03   8.17179744e-03
       1.05694658e-02   0.00000000e+00   4.49123443e-04   9.80728243e-04
       2.70933271e-03   1.61865322e-02   2.13504124e-02   0.00000000e+00
       1.17773123e-02   4.63490203e-03   1.79331966e-02   1.46366115e-02
       3.26856602e-02   4.31126006e-03   1.68787878e-02   2.02752156e-02
       1.48203062e-02   1.17346898e-03   7.87933309e-03   0.00000000e+00
       1.13274458e-03   2.25418313e-03   1.27966643e-02   2.46154526e-02
       7.15248968e-03   8.35660945e-03   3.88259360e-03   5.95428313e-03
       1.16751480e-04   5.78637193e-03   6.50575506e-03   1.47111816e-03
       1.22855215e-02   1.34294277e-02   4.03141738e-02   2.77313687e-02]

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

