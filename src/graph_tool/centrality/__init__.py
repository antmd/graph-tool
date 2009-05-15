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
    g : Graph
        Graph to be used.
    damping : float, optional (default: 0.8)
        Damping factor.
    prop : ProperyMap, optional (default: None)
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
    pagerank : A vertex property map containing the PageRank values.

    See Also
    --------
    betweenness: betweenness centrality
    eigentrust: eigentrust centrality
    absolute_trust: absolute trust centrality

    Notes
    -----
    The value of PageRank [pagerank_wikipedia]_ of vertex v :math:`PR(v)` is
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
    >>> g = gt.random_graph(100, lambda: (poisson(3), poisson(3)), seed=42)
    >>> pr = gt.pagerank(g)
    >>> print pr.get_array()
    [ 1.01514315  0.60117439  0.32514372  0.28        0.2         1.54971179
      0.28        1.0236911   0.33123536  0.4778296   0.62078363  1.25377064
      0.49213262  1.70011842  0.30671734  0.56424761  0.86810689  1.68765055
      0.49551575  0.72837655  0.39240949  1.43802363  0.51563806  0.41983927
      0.37857787  0.45875573  0.97033399  0.38531927  0.54001665  0.89328562
      0.52122532  0.94064256  1.39911631  0.64663655  1.23521006  0.71722741
      0.59460778  0.2         0.63239854  1.86292923  0.2         0.31277737
      0.74650027  0.32415672  0.47975325  1.11611173  0.53433883  0.63352435
      0.23822967  0.93151021  0.5440643   0.69188579  0.97489471  0.51216733
      1.31721331  1.32808547  0.39894203  0.50384137  0.75225633  0.28220146
      1.10818407  0.58685184  1.26437262  0.67929902  0.69678112  1.34428502
      0.61651094  0.43008378  0.7905129   1.35318411  0.2         0.2
      1.6584374   0.98009079  0.27200222  0.3413639   0.23822967  0.27963213
      1.22498499  0.34097559  0.50749002  1.21145838  0.50430676  0.50218939
      0.74232491  0.5335867   0.27254191  0.36031317  0.65344358  0.96712961
      0.53252883  0.86479464  0.59958851  0.82703737  0.68722079  0.52036384
      0.65299724  0.42291513  0.81729152  1.7586996 ]

    References
    ----------
    .. [pagerank_wikipedia] http://en.wikipedia.org/wiki/Pagerank
    .. [lawrence_pagerank_1998] P. Lawrence, B. Sergey, M. Rajeev, W. Terry,
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
    g : Graph
        Graph to be used.
    vprop : ProperyMap, optional (default: None)
        Vertex property map to store the vertex betweenness values.
    eprop : ProperyMap, optional (default: None)
        Edge property map to store the edge betweenness values.
    weight : ProperyMap, optional (default: None)
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

    The algorithm used here is defined in [brandes_faster_2001]_, and has a
    complexity of :math:`O(VE)` for unweighted graphs and :math:`O(VE + V(V+E)
    \log V)` for weighted graphs. The space complexity is :math:`O(VE)`.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    >>> from numpy.random import poisson, seed
    >>> seed(42)
    >>> g = gt.random_graph(100, lambda: (poisson(3), poisson(3)), seed=42)
    >>> vb, eb = gt.betweenness(g)
    >>> print vb.get_array()
    [ 0.03536033  0.03251351  0.00813873  0.00496977  0.          0.08339989
      0.00948258  0.05751528  0.00236377  0.00868464  0.04443961  0.04691023
      0.01768388  0.          0.01130552  0.01277964  0.04223144  0.05040177
      0.01202611  0.0012722   0.00828095  0.11598601  0.01864867  0.01412404
      0.03343004  0.01772387  0.04780278  0.01351748  0.03616999  0.09074218
      0.          0.          0.03901368  0.02526396  0.07471888  0.00219886
      0.          0.          0.01062083  0.07973799  0.          0.01410051
      0.02025676  0.          0.00988767  0.07519014  0.          0.06380861
      0.          0.01954769  0.04576145  0.04151243  0.          0.04198926
      0.0462918   0.07353227  0.00606605  0.02597097  0.02566416  0.00196642
      0.06240786  0.02996611  0.03252566  0.01451141  0.05149852  0.
      0.03582571  0.04600123  0.03776439  0.03326425  0.          0.
      0.11568858  0.01361223  0.00515358  0.007151    0.00241302  0.00271168
      0.01780978  0.01867583  0.02020758  0.01254292  0.00054971  0.00698211
      0.02359226  0.0385241   0.00157871  0.00576513  0.04173662  0.03233332
      0.0208791   0.02286212  0.04366053  0.03701801  0.02142117  0.03099565
      0.02555676  0.03365458  0.03542124  0.06174975]

    References
    ----------
    .. [betweenness_wikipedia] http://en.wikipedia.org/wiki/Centrality#Betweenness_centrality
    .. [brandes_faster_2001] U. Brandes, "A faster algorithm for betweenness
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
    g : Graph
        Graph to be used.
    betweenness : ProperyMap
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
    centrality; then, the central point dominance [freeman_set_1977]_ is defined
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
    >>> g = gt.random_graph(100, lambda: (poisson(3), poisson(3)), seed=42)
    >>> vb, eb = gt.betweenness(g)
    >>> print gt.central_point_dominance(g, vb)
    0.0902382147799

    References
    ----------
    .. [freeman_set_1977] Linton C. Freeman, "A Set of Measures of Centrality
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
    g : Graphs
        Graph to be used.
    trust_map : ProperyMap
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
    The eigentrust _[kamvar_eigentrust_2003] values :math:`t_i` correspond the
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
    >>> g = gt.random_graph(100, lambda: (poisson(3), poisson(3)), seed=42)
    >>> trust = g.new_edge_property("double")
    >>> trust.get_array()[:] = random(g.num_edges())*42
    >>> t = gt.eigentrust(g, trust, norm=True)
    >>> print t.get_array()
    [  1.86170852e-02   3.54163528e-03   6.09712602e-04   0.00000000e+00
       0.00000000e+00   3.49545179e-02   0.00000000e+00   2.59814288e-02
       8.41396546e-04   4.78599541e-03   1.01228999e-02   1.43178181e-02
       2.24766294e-03   1.80046830e-02   3.55639433e-03   4.24974765e-03
       1.11631004e-02   3.12332820e-02   6.70174456e-03   1.09689726e-02
       5.42202976e-03   2.51547994e-02   6.87197775e-03   3.90316493e-03
       2.81858126e-03   6.26514036e-03   1.12322993e-02   4.35905749e-03
       1.86938930e-02   1.93055029e-02   3.25522183e-03   9.48081499e-03
       1.84882500e-02   8.17367673e-03   4.02113149e-02   1.07092572e-02
       1.02184616e-02   0.00000000e+00   4.21126174e-03   3.97005433e-02
       0.00000000e+00   6.23025347e-04   1.92797472e-02   5.22705075e-04
       4.07751175e-03   2.11704089e-02   7.49484415e-03   8.10935540e-03
       9.47352873e-05   1.74518912e-02   1.18865927e-02   8.49808309e-03
       8.07449129e-03   6.04464513e-03   1.31497182e-02   1.61277706e-02
       3.45965628e-03   9.28003800e-03   5.81189874e-03   2.67273946e-03
       1.33359267e-02   3.99664807e-03   1.45641237e-02   2.06551771e-03
       1.89334085e-02   2.44376969e-02   7.44521415e-03   6.35266998e-03
       9.90439343e-03   2.61315207e-02   0.00000000e+00   0.00000000e+00
       4.08351424e-02   1.21805039e-02   3.45041723e-03   1.84601840e-03
       1.09623699e-03   2.37115682e-03   1.70221593e-02   4.57709422e-03
       4.21193747e-03   2.26493986e-02   3.92636239e-03   2.42441556e-03
       7.41276227e-03   7.01899189e-03   3.30982461e-03   4.18470116e-04
       8.46801514e-03   9.05050341e-03   5.09784610e-03   3.20304076e-02
       6.71276214e-03   5.26109355e-03   5.29170118e-03   3.46248974e-03
       1.10436337e-02   2.20158077e-03   1.26859707e-02   2.25728004e-02]

    References
    ----------
    .. [kamvar_eigentrust_2003] S. D. Kamvar, M. T. Schlosser, H. Garcia-Molina
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

def absolute_trust(g, trust_map, source=None, vprop=None, epslon=0.001,
                   max_iter=None, reversed=False, seed=None, ret_iter=False):
    r"""
    Samples the absolute trust centrality of each vertex in the graph, or only
    for a given source, if one is provided.

    Parameters
    ----------
    g : Graphs
        Graph to be used.
    trust_map : ProperyMap
        Edge property map with the values of trust associated with each
        edge. The values must lie in the range [0,1].
    source : Vertex, optional (default: None)
        A vertex which is used the as the sole source for gathering trust
        values, instead of all the vertices in the graph.
    vprop : PropertyMap, optional (default: None)
        Vertex property map where the values of eigentrust must be stored.
    epslon : float, optional (default: 0.001)
        Convergence condition. The iteration will stop if the total delta of all
        vertices are below this value.
    max_iter : int, optional (default: None)
        If supplied, this will limit the total number of iterations.
    reversed : bool, optional (default: False)
        Calculates the "reversed" trust instead: The direction of the edges are
        inverted, but the path weighting is preserved in the original direction
        (see Notes below).
    seed : int, optional (default: None)
         The initializing seed for the random number generator. If not supplied
         a different random value will be chosen each time.
    ret_iter : bool, optional (default: False)
        If true, the total number of iterations is also returned.

    Returns
    -------
    absolute_trust : PropertyMap
        A vertex property map containing the absolute trust vector from the
        corresponding vertex to the rest of the network. Each e lement i of the
        vector is the trust value of the vertex with index i, from the given
        vertex.

        If the parameter "source" is specified, the values of the
        property map are scalars, instead of vectors.

    See Also
    --------
    eigentrust: eigentrust centrality
    betweenness: betweenness centrality
    pagerank: PageRank centrality

    Notes
    -----
    The absolute trust between vertices i and j is defined as

    .. math::

        t_{ij} = \frac{1}{\sum_{\{i\to j\}}w_{\{i\to j\}}}\sum_{\{i\to j\}}
                 w_{\{i\to j\}} \prod_{e\in \{i\to j\}}c_e

    where the sum is taken over all paths from i to j (without loops),
    :math:`c_e` is the direct trust value associated with edge e, and
    :math:`w_{\{i\to j\}}` is the weight of a given path, which is defined as

    .. math::

       w_{\{i\to j\}} = \prod_{e\in \{i\to j\}}\{\delta_{t(e),j}(1-c_e) + c_e\},

    such that the direct trust of the last edge on the path is not considered.

    The algorithm progressively samples all possible paths, until the trust
    values converge, and has a topology-dependent complexity.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    >>> from numpy.random import poisson, random, seed
    >>> seed(42)
    >>> g = gt.random_graph(100, lambda: (poisson(3), poisson(3)), seed=42)
    >>> trust = g.new_edge_property("double")
    >>> trust.get_array()[:] = random(g.num_edges())
    >>> t = gt.absolute_trust(g, trust)
    >>> print array(t[g.vertex(10)])
    [ 0.          0.          0.          0.          0.          0.          0.
      0.          0.          0.          0.          0.          0.          0.
      0.          0.          0.          0.          0.          0.          0.
      0.          0.          0.          0.          0.          0.          0.
      0.          0.          0.          0.          0.          0.
      0.83889809  0.          0.          0.          0.          0.26116067
      0.          0.          0.          0.          0.          0.          0.
      0.          0.          0.          0.54758746  0.          0.          0.
      0.          0.          0.          0.          0.          0.          0.
      0.          0.          0.          0.          0.          0.          0.
      0.          0.04949872  0.          0.          0.          0.          0.
      0.          0.          0.          0.05011676  0.          0.          0.
      0.          0.          0.          0.43799071  0.          0.          0.
      0.          0.          0.04054452  0.          0.          0.          0.
      0.          0.          0.          0.17447756]
    """

    if seed != 0:
        seed = numpy.random.randint(0, sys.maxint)
    if vprop == None:
        if source == None:
            vprop = g.new_vertex_property("vector<double>")
        else:
            vprop = g.new_vertex_property("double")

    if source != None:
        vprop_temp = vprop
        vprop = g.new_vertex_property("vector<double>")
        source = g.vertex_index[source]
    else:
        source = -1

    if max_iter == None:
        max_iter = 0

    if reversed:
        g.stash_filter(reversed=True)
        g.set_reversed(True)

    ic = libgraph_tool_centrality.\
            get_absolute_trust(g._Graph__graph, source,
                               _prop("e", g, trust_map), _prop("v", g, vprop),
                               epslon, max_iter, reversed, seed)
    if reversed:
        g.pop_filter(reversed=True)

    if source != -1:
        vprop_temp.get_array()[:] = numpy.array(vprop[g.vertex(source)])
        vprop = vprop_temp

    if ret_iter:
        return vprop, ic
    else:
        return vprop

