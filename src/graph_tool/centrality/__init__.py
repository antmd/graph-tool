#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# graph_tool -- a general graph manipulation python module
#
# Copyright (C) 2006-2015 Tiago de Paula Peixoto <tiago@skewed.de>
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
   closeness
   eigenvector
   katz
   hits
   eigentrust
   trust_transitivity

Contents
++++++++
"""

from __future__ import division, absolute_import, print_function

from .. dl_import import dl_import
dl_import("from . import libgraph_tool_centrality")

from .. import _prop, ungroup_vector_property
from .. topology import shortest_distance
import sys
import numpy
import numpy.linalg

__all__ = ["pagerank", "betweenness", "central_point_dominance", "closeness",
           "eigentrust", "eigenvector", "katz", "hits", "trust_transitivity"]


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
        Vertex property map to store the PageRank values. If supplied, it will
        be used uninitialized.
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
    eigenvector: eigenvector centrality
    hits: hubs and authority centralities
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

    .. testsetup:: pagerank

       import matplotlib

    .. doctest:: pagerank

       >>> g = gt.collection.data["polblogs"]
       >>> g = gt.GraphView(g, vfilt=gt.label_largest_component(g))
       >>> pr = gt.pagerank(g)
       >>> gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=pr,
       ...               vertex_size=gt.prop_to_size(pr, mi=5, ma=15),
       ...               vorder=pr, vcmap=matplotlib.cm.gist_heat,
       ...               output="polblogs_pr.pdf")
       <...>

    .. testcode:: pagerank
       :hide:

       gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=pr,
                     vertex_size=gt.prop_to_size(pr, mi=5, ma=15),
                     vorder=pr, vcmap=matplotlib.cm.gist_heat,
                     output="polblogs_pr.png")


    .. figure:: polblogs_pr.*
       :align: center

       PageRank values of the a political blogs network of [adamic-polblogs]_.

    Now with a personalization vector, and edge weights:

    .. doctest:: pagerank

       >>> d = g.degree_property_map("total")
       >>> periphery = d.a <= 2
       >>> p = g.new_vertex_property("double")
       >>> p.a[periphery] = 100
       >>> pr = gt.pagerank(g, pers=p)
       >>> gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=pr,
       ...               vertex_size=gt.prop_to_size(pr, mi=5, ma=15),
       ...               vorder=pr, vcmap=matplotlib.cm.gist_heat,
       ...               output="polblogs_pr_pers.pdf")
       <...>

    .. testcode:: pagerank
       :hide:

       gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=pr,
                     vertex_size=gt.prop_to_size(pr, mi=5, ma=15),
                     vcmap=matplotlib.cm.gist_heat,
                     vorder=pr, output="polblogs_pr_pers.png")


    .. figure:: polblogs_pr_pers.*
       :align: center

       Personalized PageRank values of the a political blogs network of
       [adamic-polblogs]_, where vertices with very low degree are given
       artificially high scores.

    References
    ----------
    .. [pagerank-wikipedia] http://en.wikipedia.org/wiki/Pagerank
    .. [lawrence-pagerank-1998] P. Lawrence, B. Sergey, M. Rajeev, W. Terry,
       "The pagerank citation ranking: Bringing order to the web", Technical
       report, Stanford University, 1998
    .. [Langville-survey-2005] A. N. Langville, C. D. Meyer, "A Survey of
       Eigenvector Methods for Web Information Retrieval", SIAM Review, vol. 47,
       no. 1, pp. 135-161, 2005, :DOI:`10.1137/S0036144503424786`
    .. [adamic-polblogs] L. A. Adamic and N. Glance, "The political blogosphere
       and the 2004 US Election", in Proceedings of the WWW-2005 Workshop on the
       Weblogging Ecosystem (2005). :DOI:`10.1145/1134271.1134277`
    """

    if max_iter == None:
        max_iter = 0
    if prop == None:
        prop = g.new_vertex_property("double")
        N = len(prop.fa)
        prop.fa = pers.fa[:N] if pers is not None else 1. / g.num_vertices()
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
    eigenvector: eigenvector centrality
    hits: hubs and authority centralities
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

    .. testsetup:: betweenness

       import matplotlib

    .. doctest:: betweenness

       >>> g = gt.collection.data["polblogs"]
       >>> g = gt.GraphView(g, vfilt=gt.label_largest_component(g))
       >>> vp, ep = gt.betweenness(g)
       >>> gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=vp,
       ...               vertex_size=gt.prop_to_size(vp, mi=5, ma=15),
       ...               edge_pen_width=gt.prop_to_size(ep, mi=0.5, ma=5),
       ...               vcmap=matplotlib.cm.gist_heat,
       ...               vorder=vp, output="polblogs_betweenness.pdf")
       <...>

    .. testcode:: betweenness
       :hide:

       gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=vp,
                     vertex_size=gt.prop_to_size(vp, mi=5, ma=15),
                     edge_pen_width=gt.prop_to_size(ep, mi=0.5, ma=5),
                     vcmap=matplotlib.cm.gist_heat,
                     vorder=vp, output="polblogs_betweenness.png")


    .. figure:: polblogs_betweenness.*
       :align: center

       Betweenness values of the a political blogs network of [adamic-polblogs]_.

    References
    ----------
    .. [betweenness-wikipedia] http://en.wikipedia.org/wiki/Centrality#Betweenness_centrality
    .. [brandes-faster-2001] U. Brandes, "A faster algorithm for betweenness
       centrality", Journal of Mathematical Sociology, 2001, :doi:`10.1080/0022250X.2001.9990249`
    .. [adamic-polblogs] L. A. Adamic and N. Glance, "The political blogosphere
       and the 2004 US Election", in Proceedings of the WWW-2005 Workshop on the
       Weblogging Ecosystem (2005). :DOI:`10.1145/1134271.1134277`
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

def closeness(g, weight=None, source=None, vprop=None, norm=True, harmonic=False):
    r"""
    Calculate the closeness centrality for each vertex.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    weight : :class:`~graph_tool.PropertyMap`, optional (default: None)
        Edge property map corresponding to the weight value of each edge.
    source : :class:`~graph_tool.Vertex`, optional (default: ``None``)
        If specified, the centrality is computed for this vertex alone.
    vprop : :class:`~graph_tool.PropertyMap`, optional (default: ``None``)
        Vertex property map to store the vertex centrality values.
    norm : bool, optional (default: ``True``)
        Whether or not the centrality values should be normalized.
    harmonic : bool, optional (default: ``False``)
        If true, the sum of the inverse of the distances will be computed,
        instead of the inverse of the sum.

    Returns
    -------
    vertex_closeness : :class:`~graph_tool.PropertyMap`
        A vertex property map with the vertex closeness values.

    See Also
    --------
    central_point_dominance: central point dominance of the graph
    pagerank: PageRank centrality
    eigentrust: eigentrust centrality
    eigenvector: eigenvector centrality
    hits: hubs and authority centralities
    trust_transitivity: pervasive trust transitivity

    Notes
    -----
    The closeness centrality of a vertex :math:`i` is defined as,

    .. math::

        c_i = \frac{1}{\sum_j d_{ij}}

    where :math:`d_{ij}` is the (possibly directed and/or weighted) distance
    from :math:`i` to :math:`j`. In case there is no path between the two
    vertices, here the distance is taken to be zero.

    If ``harmonic == True``, the definition becomes

    .. math::

        c_i = \sum_j\frac{1}{d_{ij}},

    but now, in case there is no path between the two vertices, we take
    :math:`d_{ij} \to\infty` such that :math:`1/d_{ij}=0`.

    If ``norm == True``, the values of :math:`c_i` are normalized by
    :math:`n_i-1` where :math:`n_i` is the size of the (out-) component of
    :math:`i`. If ``harmonic == True``, they are instead simply normalized by
    :math:`N-1`.

    The algorithm complexity of :math:`O(N(N + E))` for unweighted graphs and
    :math:`O(N(N+E) \log N)` for weighted graphs. If the option ``source`` is
    specified, this drops to :math:`O(N + E)` and :math:`O((N+E)\log N)`
    respectively.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------

    .. testsetup:: closeness

       import matplotlib

    .. doctest:: closeness

       >>> g = gt.collection.data["polblogs"]
       >>> g = gt.GraphView(g, vfilt=gt.label_largest_component(g))
       >>> c = gt.closeness(g)
       >>> gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=c,
       ...               vertex_size=gt.prop_to_size(c, mi=5, ma=15),
       ...               vcmap=matplotlib.cm.gist_heat,
       ...               vorder=c, output="polblogs_closeness.pdf")
       <...>

    .. testcode:: closeness
       :hide:

       gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=c,
                     vertex_size=gt.prop_to_size(c, mi=5, ma=15),
                     vcmap=matplotlib.cm.gist_heat,
                     vorder=c, output="polblogs_closeness.png")


    .. figure:: polblogs_closeness.*
       :align: center

       Closeness values of the a political blogs network of [adamic-polblogs]_.

    References
    ----------
    .. [closeness-wikipedia] https://en.wikipedia.org/wiki/Closeness_centrality
    .. [opsahl-node-2010] Opsahl, T., Agneessens, F., Skvoretz, J., "Node
       centrality in weighted networks: Generalizing degree and shortest
       paths". Social Networks 32, 245-251, 2010 :DOI:`10.1016/j.socnet.2010.03.006`
    .. [adamic-polblogs] L. A. Adamic and N. Glance, "The political blogosphere
       and the 2004 US Election", in Proceedings of the WWW-2005 Workshop on the
       Weblogging Ecosystem (2005). :DOI:`10.1145/1134271.1134277`

    """
    if source is None:
        if vprop == None:
            vprop = g.new_vertex_property("double")
        libgraph_tool_centrality.\
            closeness(g._Graph__graph, _prop("e", g, weight),
                      _prop("v", g, vprop), harmonic, norm)
        return vprop
    else:
        max_dist = g.num_vertices() + 1
        dist = shortest_distance(g, source=source, weights=weight,
                                 max_dist=max_dist)
        dists = dist.fa[(dist.fa < max_dist) * (dist.fa > 0)]
        if harmonic:
            c = (1. / dists).sum()
            if norm:
                c /= g.num_vertices() - 1
        else:
            c = 1. / dists.sum()
            if norm:
                c *= len(dists)
        return c


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

    >>> g = gt.collection.data["polblogs"]
    >>> g = gt.GraphView(g, vfilt=gt.label_largest_component(g))
    >>> vp, ep = gt.betweenness(g)
    >>> print(gt.central_point_dominance(g, vp))
    0.11610685614353008

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
        Vertex property map where the values of eigenvector must be stored. If
        provided, it will be used uninitialized.
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
    hits: hubs and authority centralities
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

    .. testsetup:: eigenvector

       np.random.seed(42)
       import matplotlib

    .. doctest:: eigenvector

       >>> g = gt.collection.data["polblogs"]
       >>> g = gt.GraphView(g, vfilt=gt.label_largest_component(g))
       >>> w = g.new_edge_property("double")
       >>> w.a = np.random.random(len(w.a)) * 42
       >>> ee, x = gt.eigenvector(g, w)
       >>> print(ee)
       724.3027459221508
       >>> gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=x,
       ...               vertex_size=gt.prop_to_size(x, mi=5, ma=15),
       ...               vcmap=matplotlib.cm.gist_heat,
       ...               vorder=x, output="polblogs_eigenvector.pdf")
       <...>

    .. testcode:: eigenvector
       :hide:

       gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=x,
                     vertex_size=gt.prop_to_size(x, mi=5, ma=15),
                     vcmap=matplotlib.cm.gist_heat,
                     vorder=x, output="polblogs_eigenvector.png")


    .. figure:: polblogs_eigenvector.*
       :align: center

       Eigenvector values of the a political blogs network of
       [adamic-polblogs]_, with random weights attributed to the edges.

    References
    ----------

    .. [eigenvector-centrality] http://en.wikipedia.org/wiki/Centrality#Eigenvector_centrality
    .. [power-method] http://en.wikipedia.org/wiki/Power_iteration
    .. [langville-survey-2005] A. N. Langville, C. D. Meyer, "A Survey of
       Eigenvector Methods for Web Information Retrieval", SIAM Review, vol. 47,
       no. 1, pp. 135-161, 2005, :DOI:`10.1137/S0036144503424786`
    .. [adamic-polblogs] L. A. Adamic and N. Glance, "The political blogosphere
       and the 2004 US Election", in Proceedings of the WWW-2005 Workshop on the
       Weblogging Ecosystem (2005). :DOI:`10.1145/1134271.1134277`

    """

    if vprop is None:
        vprop = g.new_vertex_property("double")
        vprop.fa = 1. / g.num_vertices()
    if max_iter is None:
        max_iter = 0
    ee = libgraph_tool_centrality.\
         get_eigenvector(g._Graph__graph, _prop("e", g, weight),
                         _prop("v", g, vprop), epsilon, max_iter)
    return ee, vprop


def katz(g, alpha=0.01, beta=None, weight=None, vprop=None, epsilon=1e-6,
         max_iter=None, norm=True):
    r"""
    Calculate the Katz centrality of each vertex in the graph.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    weight : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Edge property map with the edge weights.
    alpha : float, optional (default: ``0.01``)
        Free parameter :math:`\alpha`. This must be smaller than the inverse of
        the largest eigenvalue of the adjacency matrix.
    beta : :class:`~graph_tool.PropertyMap`, optional (default: ``None``)
        Vertex property map where the local personalization values. If not
        provided, the global value of 1 will be used.
    vprop : :class:`~graph_tool.PropertyMap`, optional (default: ``None``)
        Vertex property map where the values of eigenvector must be stored. If
        provided, it will be used uninitialized.
    epsilon : float, optional (default: ``1e-6``)
        Convergence condition. The iteration will stop if the total delta of all
        vertices are below this value.
    max_iter : int, optional (default: ``None``)
        If supplied, this will limit the total number of iterations.
    norm : bool, optional (default: ``True``)
        Whether or not the centrality values should be normalized.

    Returns
    -------
    centrality : :class:`~graph_tool.PropertyMap`
        A vertex property map containing the Katz centrality values.

    See Also
    --------
    betweenness: betweenness centrality
    pagerank: PageRank centrality
    eigenvector: eigenvector centrality
    hits: hubs and authority centralities
    trust_transitivity: pervasive trust transitivity

    Notes
    -----

    The Katz centrality :math:`\mathbf{x}` is the solution of the nonhomogeneous
    linear system

    .. math::

        \mathbf{x} = \alpha\mathbf{A}\mathbf{x} + \mathbf{\beta},


    where :math:`\mathbf{A}` is the (weighted) adjacency matrix and
    :math:`\mathbf{\beta}` is the personalization vector (if not supplied,
    :math:`\mathbf{\beta} = \mathbf{1}` is assumed).

    The algorithm uses successive iterations of the equation above, which has a
    topology-dependent convergence complexity.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    .. testsetup:: katz

       np.random.seed(42)
       import matplotlib

    .. doctest:: katz

       >>> g = gt.collection.data["polblogs"]
       >>> g = gt.GraphView(g, vfilt=gt.label_largest_component(g))
       >>> w = g.new_edge_property("double")
       >>> w.a = np.random.random(len(w.a))
       >>> x = gt.katz(g, weight=w)
       >>> gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=x,
       ...               vertex_size=gt.prop_to_size(x, mi=5, ma=15),
       ...               vcmap=matplotlib.cm.gist_heat,
       ...               vorder=x, output="polblogs_katz.pdf")
       <...>

    .. testcode:: katz
       :hide:

       gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=x,
                     vertex_size=gt.prop_to_size(x, mi=5, ma=15),
                     vcmap=matplotlib.cm.gist_heat,
                     vorder=x, output="polblogs_katz.png")


    .. figure:: polblogs_katz.*
       :align: center

       Katz centrality values of the a political blogs network of
       [adamic-polblogs]_, with random weights attributed to the edges.

    References
    ----------

    .. [katz-centrality] http://en.wikipedia.org/wiki/Katz_centrality
    .. [katz-new] L. Katz, "A new status index derived from sociometric analysis",
       Psychometrika 18, Number 1, 39-43, 1953, :DOI:`10.1007/BF02289026`
    .. [adamic-polblogs] L. A. Adamic and N. Glance, "The political blogosphere
       and the 2004 US Election", in Proceedings of the WWW-2005 Workshop on the
       Weblogging Ecosystem (2005). :DOI:`10.1145/1134271.1134277`
    """

    if vprop is None:
        vprop = g.new_vertex_property("double")
    if max_iter is None:
        max_iter = 0
    libgraph_tool_centrality.\
         get_katz(g._Graph__graph, _prop("e", g, weight), _prop("v", g, vprop),
                  _prop("v", g, beta), float(alpha), epsilon, max_iter)
    if norm:
        vprop.fa = vprop.fa / numpy.linalg.norm(vprop.fa)
    return vprop


def hits(g, weight=None, xprop=None, yprop=None, epsilon=1e-6, max_iter=None):
    r"""
    Calculate the authority and hub centralities of each vertex in the graph.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    weight : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Edge property map with the edge weights.
    xprop : :class:`~graph_tool.PropertyMap`, optional (default: ``None``)
        Vertex property map where the authority centrality must be stored.
    yprop : :class:`~graph_tool.PropertyMap`, optional (default: ``None``)
        Vertex property map where the hub centrality must be stored.
    epsilon : float, optional (default: ``1e-6``)
        Convergence condition. The iteration will stop if the total delta of all
        vertices are below this value.
    max_iter : int, optional (default: ``None``)
        If supplied, this will limit the total number of iterations.

    Returns
    -------
    eig : `float`
        The largest eigenvalue of the cocitation matrix.
    x : :class:`~graph_tool.PropertyMap`
        A vertex property map containing the authority centrality values.
    y : :class:`~graph_tool.PropertyMap`
        A vertex property map containing the hub centrality values.

    See Also
    --------
    betweenness: betweenness centrality
    eigenvector: eigenvector centrality
    pagerank: PageRank centrality
    trust_transitivity: pervasive trust transitivity

    Notes
    -----

    The Hyperlink-Induced Topic Search (HITS) centrality assigns hub
    (:math:`\mathbf{y}`) and authority (:math:`\mathbf{x}`) centralities to the
    vertices, following:

    .. math::

        \begin{align}
            \mathbf{x} &= \alpha\mathbf{A}\mathbf{y} \\
            \mathbf{y} &= \beta\mathbf{A}^T\mathbf{x}
        \end{align}


    where :math:`\mathbf{A}` is the (weighted) adjacency matrix and
    :math:`\lambda = 1/(\alpha\beta)` is the largest eigenvalue of the
    cocitation matrix, :math:`\mathbf{A}\mathbf{A}^T`. (Without loss of
    generality, we set :math:`\beta=1` in the algorithm.)

    The algorithm uses the power method which has a topology-dependent complexity of
    :math:`O\left(N\times\frac{-\log\epsilon}{\log|\lambda_1/\lambda_2|}\right)`,
    where :math:`N` is the number of vertices, :math:`\epsilon` is the ``epsilon``
    parameter, and :math:`\lambda_1` and :math:`\lambda_2` are the largest and
    second largest eigenvalues of the (weighted) cocitation matrix, respectively.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------

    .. testsetup:: hits

       import matplotlib

    .. doctest:: hits

       >>> g = gt.collection.data["polblogs"]
       >>> g = gt.GraphView(g, vfilt=gt.label_largest_component(g))
       >>> ee, x, y = gt.hits(g)
       >>> gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=x,
       ...               vertex_size=gt.prop_to_size(x, mi=5, ma=15),
       ...               vcmap=matplotlib.cm.gist_heat,
       ...               vorder=x, output="polblogs_hits_auths.pdf")
       <...>
       >>> gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=y,
       ...               vertex_size=gt.prop_to_size(y, mi=5, ma=15),
       ...               vcmap=matplotlib.cm.gist_heat,
       ...               vorder=y, output="polblogs_hits_hubs.pdf")
       <...>

    .. testcode:: hits
       :hide:

       gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=x,
                     vertex_size=gt.prop_to_size(x, mi=5, ma=15),
                     vcmap=matplotlib.cm.gist_heat,
                     vorder=x, output="polblogs_hits_auths.png")
       gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=y,
                     vertex_size=gt.prop_to_size(y, mi=5, ma=15),
                     vcmap=matplotlib.cm.gist_heat,
                     vorder=y, output="polblogs_hits_hubs.png")


    .. figure:: polblogs_hits_auths.*
       :align: center

       HITS authority values of the a political blogs network of
       [adamic-polblogs]_.

    .. figure:: polblogs_hits_hubs.*
       :align: center

       HITS hub values of the a political blogs network of [adamic-polblogs]_.

    References
    ----------

    .. [hits-algorithm] http://en.wikipedia.org/wiki/HITS_algorithm
    .. [kleinberg-authoritative] J. Kleinberg, "Authoritative sources in a
       hyperlinked environment", Journal of the ACM 46 (5): 604-632, 1999,
       :DOI:`10.1145/324133.324140`.
    .. [power-method] http://en.wikipedia.org/wiki/Power_iteration
    .. [adamic-polblogs] L. A. Adamic and N. Glance, "The political blogosphere
       and the 2004 US Election", in Proceedings of the WWW-2005 Workshop on the
       Weblogging Ecosystem (2005). :DOI:`10.1145/1134271.1134277`
    """

    if xprop is None:
        xprop = g.new_vertex_property("double")
    if yprop is None:
        yprop = g.new_vertex_property("double")
    if max_iter is None:
        max_iter = 0
    l = libgraph_tool_centrality.\
         get_hits(g._Graph__graph, _prop("e", g, weight), _prop("v", g, xprop),
                  _prop("v", g, yprop), epsilon, max_iter)
    return 1. / l, xprop, yprop


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

    .. testsetup:: eigentrust

       np.random.seed(42)
       import matplotlib

    .. doctest:: eigentrust

       >>> g = gt.collection.data["polblogs"]
       >>> g = gt.GraphView(g, vfilt=gt.label_largest_component(g))
       >>> w = g.new_edge_property("double")
       >>> w.a = np.random.random(len(w.a)) * 42
       >>> t = gt.eigentrust(g, w)
       >>> gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=t,
       ...               vertex_size=gt.prop_to_size(t, mi=5, ma=15),
       ...               vcmap=matplotlib.cm.gist_heat,
       ...               vorder=t, output="polblogs_eigentrust.pdf")
       <...>

    .. testcode:: eigentrust
       :hide:

       gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=t,
                     vertex_size=gt.prop_to_size(t, mi=5, ma=15),
                     vcmap=matplotlib.cm.gist_heat,
                     vorder=t, output="polblogs_eigentrust.png")


    .. figure:: polblogs_eigentrust.*
       :align: center

       Eigentrust values of the a political blogs network of
       [adamic-polblogs]_, with random weights attributed to the edges.


    References
    ----------
    .. [kamvar-eigentrust-2003] S. D. Kamvar, M. T. Schlosser, H. Garcia-Molina
       "The eigentrust algorithm for reputation management in p2p networks",
       Proceedings of the 12th international conference on World Wide Web,
       Pages: 640 - 651, 2003, :doi:`10.1145/775152.775242`
    .. [adamic-polblogs] L. A. Adamic and N. Glance, "The political blogosphere
       and the 2004 US Election", in Proceedings of the WWW-2005 Workshop on the
       Weblogging Ecosystem (2005). :DOI:`10.1145/1134271.1134277`
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
    .. testsetup:: trust_transitivity

       np.random.seed(42)
       import matplotlib

    .. doctest:: trust_transitivity

       >>> g = gt.collection.data["polblogs"]
       >>> g = gt.GraphView(g, vfilt=gt.label_largest_component(g))
       >>> g = gt.Graph(g, prune=True)
       >>> w = g.new_edge_property("double")
       >>> w.a = np.random.random(len(w.a))
       >>> g.vp["label"][g.vertex(42)]
       'blogforamerica.com'
       >>> t = gt.trust_transitivity(g, w, source=g.vertex(42))
       >>> gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=t,
       ...               vertex_size=gt.prop_to_size(t, mi=5, ma=15),
       ...               vcmap=matplotlib.cm.gist_heat,
       ...               vorder=t, output="polblogs_trust_transitivity.pdf")
       <...>

    .. testcode:: trust_transitivity
       :hide:

       gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=t,
                     vertex_size=gt.prop_to_size(t, mi=5, ma=15),
                     vcmap=matplotlib.cm.gist_heat,
                     vorder=t, output="polblogs_trust_transitivity.png")


    .. figure:: polblogs_trust_transitivity.*
       :align: center

       Trust transitivity values from source vertex 42 of the a political blogs
       network of [adamic-polblogs]_, with random weights attributed to the
       edges.

    References
    ----------
    .. [richters-trust-2010] Oliver Richters and Tiago P. Peixoto, "Trust
       Transitivity in Social Networks," PLoS ONE 6, no. 4:
       e1838 (2011), :doi:`10.1371/journal.pone.0018384`
    .. [adamic-polblogs] L. A. Adamic and N. Glance, "The political blogosphere
       and the 2004 US Election", in Proceedings of the WWW-2005 Workshop on the
       Weblogging Ecosystem (2005). :DOI:`10.1145/1134271.1134277`

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
