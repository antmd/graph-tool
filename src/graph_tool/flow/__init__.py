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
``graph_tool.flow`` - Maximum flow algorithms
---------------------------------------------

Summary
+++++++

.. autosummary::
   :nosignatures:

   edmonds_karp_max_flow
   push_relabel_max_flow
   kolmogorov_max_flow
   max_cardinality_matching

Contents
++++++++
"""

from .. dl_import import dl_import
dl_import("import libgraph_tool_flow")

from .. core import _prop, _check_prop_scalar, _check_prop_writable
__all__ = ["edmonds_karp_max_flow", "push_relabel_max_flow",
           "kolmogorov_max_flow", "max_cardinality_matching"]


def edmonds_karp_max_flow(g, source, target, capacity, residual=None):
    r"""
    Calculate maximum flow on the graph with Edmonds-Karp algorithm.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    source : Vertex
        The source vertex.
    target : Vertex
        The target (or "sink") vertex.
    capacity : :class:`~graph_tool.PropertyMap`
        Edge property map with the edge capacities.
    residual : :class:`~graph_tool.PropertyMap` (optional, default: none)
        Edge property map where the residuals should be stored.

    Returns
    -------
    residual : :class:`~graph_tool.PropertyMap`
        Edge property map with the residual capacities (capacity - flow).

    Notes
    -----
    The algorithm is due to [edmonds-theoretical-1972]_, though we are using the
    variation called the "labeling algorithm" described in
    [ravindra-network-1993]_.

    This algorithm provides a very simple and easy to implement solution to the
    maximum flow problem. However, there are several reasons why this algorithm
    is not as good as the push_relabel_max_flow() or the kolmogorov_max_flow()
    algorithm.

    - In the non-integer capacity case, the time complexity is :math:`O(V E^2)`
      which is worse than the time complexity of the push-relabel algorithm
      :math:`O(V^2E^{1/2})` for all but the sparsest of graphs.

    - In the integer capacity case, if the capacity bound U is very large then
      the algorithm will take a long time.

    Examples
    --------
    >>> from numpy.random import seed, random
    >>> seed(43)
    >>> g = gt.random_graph(100, lambda: (2,2))
    >>> c = g.new_edge_property("double")
    >>> c.a = random(len(c.a))
    >>> res = gt.edmonds_karp_max_flow(g, g.vertex(0), g.vertex(1), c)
    >>> res.a = c.a - res.a  # the actual flow
    >>> print res.a[0:g.num_edges()]
    [ 0.13339096  0.24058962  0.          0.          0.          0.          0.
      0.          0.31026953  0.          0.          0.          0.          0.
      0.10719866  0.          0.          0.          0.          0.          0.
      0.          0.13339096  0.          0.          0.          0.10719866
      0.          0.          0.          0.          0.10719866  0.          0.
      0.          0.          0.          0.          0.          0.          0.
      0.          0.          0.24058962  0.          0.          0.          0.
      0.          0.          0.          0.          0.          0.          0.
      0.          0.          0.          0.10719866  0.          0.          0.
      0.          0.          0.          0.          0.31026953  0.          0.
      0.          0.          0.44366049  0.          0.          0.          0.
      0.          0.          0.          0.          0.          0.
      0.10719866  0.          0.          0.          0.          0.          0.
      0.          0.          0.          0.          0.          0.          0.
      0.          0.          0.          0.          0.          0.          0.
      0.          0.10719866  0.          0.          0.          0.          0.
      0.          0.          0.          0.          0.          0.          0.
      0.          0.          0.          0.          0.          0.          0.
      0.          0.          0.          0.          0.          0.          0.
      0.          0.          0.          0.          0.          0.          0.
      0.          0.          0.          0.10719866  0.          0.          0.
      0.          0.10719866  0.          0.44366049  0.          0.          0.
      0.          0.          0.          0.          0.          0.          0.
      0.          0.          0.          0.          0.          0.          0.
      0.          0.          0.          0.          0.          0.10719866
      0.          0.          0.31026953  0.          0.31026953  0.          0.
      0.          0.          0.          0.          0.          0.          0.
      0.          0.          0.          0.          0.          0.          0.
      0.          0.          0.          0.          0.          0.          0.        ]

    References
    ----------
    .. [boost-edmonds-karp] http://www.boost.org/libs/graph/doc/edmonds_karp_max_flow.html
    .. [edmonds-theoretical-1972] Jack Edmonds and Richard M. Karp, "Theoretical
       improvements in the algorithmic efficiency for network flow problems.
       Journal of the ACM", 1972 19:248-264
    .. [ravindra-network-1993] Ravindra K. Ahuja and Thomas L. Magnanti and
       James B. Orlin,"Network Flows: Theory, Algorithms, and Applications".
       Prentice Hall, 1993.
    """

    _check_prop_scalar(capacity, "capacity")
    if residual == None:
        residual = g.new_edge_property(capacity.value_type())
    _check_prop_scalar(residual, "residual")
    _check_prop_writable(residual, "residual")

    libgraph_tool_flow.\
           edmonds_karp_max_flow(g._Graph__graph, int(source), int(target),
                                 _prop("e", g, capacity),
                                 _prop("e", g, residual))
    return residual


def push_relabel_max_flow(g, source, target, capacity, residual=None):
    r"""
    Calculate maximum flow on the graph with push-relabel algorithm.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    source : Vertex
        The source vertex.
    target : Vertex
        The target (or "sink") vertex.
    capacity : :class:`~graph_tool.PropertyMap`
        Edge property map with the edge capacities.
    residual : :class:`~graph_tool.PropertyMap` (optional, default: none)
        Edge property map where the residuals should be stored.

    Returns
    -------
    residual : :class:`~graph_tool.PropertyMap`
        Edge property map with the residual capacities (capacity - flow).

    Notes
    -----
    The algorithm is defined in [goldberg-new-1985]_. The complexity is
    :math:`O(V^3)`.

    Examples
    --------
    >>> from numpy.random import seed, random
    >>> seed(43)
    >>> g = gt.random_graph(100, lambda: (2,2))
    >>> c = g.new_edge_property("double")
    >>> c.a = random(len(c.a))
    >>> res = gt.push_relabel_max_flow(g, g.vertex(0), g.vertex(1), c)
    >>> res.a = c.a - res.a  # the actual flow
    >>> print res.a[0:g.num_edges()]
    [ 0.13339096  0.24058962  0.          0.01254176  0.          0.          0.
      0.          0.31026953  0.          0.00745848  0.00335702  0.03601961
      0.01254176  0.0794977   0.          0.00745848  0.          0.06036354
      0.01254176  0.00335702  0.          0.14084944  0.          0.03601961
      0.00745848  0.08285472  0.          0.          0.06372057  0.
      0.10384163  0.          0.00335702  0.          0.          0.          0.
      0.          0.02434393  0.          0.          0.02434393  0.2480481   0.
      0.          0.          0.          0.0368857   0.          0.06372057
      0.          0.03601961  0.          0.00335702  0.01254176  0.02434393
      0.          0.10719866  0.04347809  0.01254176  0.          0.01254176
      0.          0.          0.00335702  0.31026953  0.03601961  0.          0.
      0.          0.45111897  0.00335702  0.          0.          0.
      0.00335702  0.          0.00335702  0.06036354  0.00745848  0.01589879
      0.0794977   0.          0.06372057  0.          0.00508328  0.          0.
      0.06036354  0.06036354  0.03601961  0.0368857   0.          0.
      0.00335702  0.          0.          0.          0.03601961  0.01589879
      0.01254176  0.          0.01254176  0.03601961  0.          0.0368857
      0.03601961  0.01589879  0.03601961  0.          0.03601961  0.          0.
      0.          0.00335702  0.01589879  0.0368857   0.03601961  0.03601961
      0.01254176  0.01254176  0.          0.00335702  0.02434393  0.00745848
      0.          0.03937663  0.06036354  0.          0.          0.          0.
      0.          0.          0.          0.02434393  0.          0.03601961
      0.          0.          0.11974042  0.          0.          0.03601961
      0.01254176  0.11465713  0.09465689  0.46800442  0.03601961  0.
      0.00335702  0.          0.          0.03601961  0.          0.
      0.00745848  0.          0.          0.03601961  0.          0.
      0.09974018  0.00335702  0.01254176  0.00335702  0.01254176  0.
      0.01254176  0.10719866  0.          0.          0.          0.31026953
      0.          0.31026953  0.          0.01254176  0.          0.          0.
      0.02000024  0.0108155   0.03601961  0.          0.03601961  0.03937663
      0.00335702  0.01254176  0.          0.          0.          0.01688546
      0.          0.          0.06036354  0.03601961  0.          0.        ]

    References
    ----------
    .. [boost-push-relabel] http://www.boost.org/libs/graph/doc/push_relabel_max_flow.html
    .. [goldberg-new-1985] A. V. Goldberg, "A New Max-Flow Algorithm",  MIT
       Tehnical report MIT/LCS/TM-291, 1985.
    """

    _check_prop_scalar(capacity, "capacity")
    if residual == None:
        residual = g.new_edge_property(capacity.value_type())
    _check_prop_scalar(residual, "residual")
    _check_prop_writable(residual, "residual")

    libgraph_tool_flow.\
           push_relabel_max_flow(g._Graph__graph, int(source), int(target),
                                 _prop("e", g, capacity),
                                 _prop("e", g, residual))
    return residual


def kolmogorov_max_flow(g, source, target, capacity, residual=None):
    r"""Calculate maximum flow on the graph with Kolmogorov algorithm.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    source : Vertex
        The source vertex.
    target : Vertex
        The target (or "sink") vertex.
    capacity : :class:`~graph_tool.PropertyMap`
        Edge property map with the edge capacities.
    residual : :class:`~graph_tool.PropertyMap` (optional, default: none)
        Edge property map where the residuals should be stored.

    Returns
    -------
    residual : :class:`~graph_tool.PropertyMap`
        Edge property map with the residual capacities (capacity - flow).

    Notes
    -----

    The algorithm is defined in [kolmogorov-graph-2003]_ and
    [boykov-experimental-2004]_. The worst case complexity is
    :math:`O(EV^2|C|)`, where :math:`|C|` is the minimum cut (but typically
    performs much better).

    For a more detailed description, see [boost-kolmogorov]_.

    Examples
    --------
    >>> from numpy.random import seed, random
    >>> seed(43)
    >>> g = gt.random_graph(100, lambda: (2,2))
    >>> c = g.new_edge_property("double")
    >>> c.a = random(len(c.a))
    >>> res = gt.push_relabel_max_flow(g, g.vertex(0), g.vertex(1), c)
    >>> res.a = c.a - res.a  # the actual flow
    >>> print res.a[0:g.num_edges()]
    [ 0.13339096  0.24058962  0.          0.01254176  0.          0.          0.
      0.          0.31026953  0.          0.00745848  0.00335702  0.03601961
      0.01254176  0.0794977   0.          0.00745848  0.          0.06036354
      0.01254176  0.00335702  0.          0.14084944  0.          0.03601961
      0.00745848  0.08285472  0.          0.          0.06372057  0.
      0.10384163  0.          0.00335702  0.          0.          0.          0.
      0.          0.02434393  0.          0.          0.02434393  0.2480481   0.
      0.          0.          0.          0.0368857   0.          0.06372057
      0.          0.03601961  0.          0.00335702  0.01254176  0.02434393
      0.          0.10719866  0.04347809  0.01254176  0.          0.01254176
      0.          0.          0.00335702  0.31026953  0.03601961  0.          0.
      0.          0.45111897  0.00335702  0.          0.          0.
      0.00335702  0.          0.00335702  0.06036354  0.00745848  0.01589879
      0.0794977   0.          0.06372057  0.          0.00508328  0.          0.
      0.06036354  0.06036354  0.03601961  0.0368857   0.          0.
      0.00335702  0.          0.          0.          0.03601961  0.01589879
      0.01254176  0.          0.01254176  0.03601961  0.          0.0368857
      0.03601961  0.01589879  0.03601961  0.          0.03601961  0.          0.
      0.          0.00335702  0.01589879  0.0368857   0.03601961  0.03601961
      0.01254176  0.01254176  0.          0.00335702  0.02434393  0.00745848
      0.          0.03937663  0.06036354  0.          0.          0.          0.
      0.          0.          0.          0.02434393  0.          0.03601961
      0.          0.          0.11974042  0.          0.          0.03601961
      0.01254176  0.11465713  0.09465689  0.46800442  0.03601961  0.
      0.00335702  0.          0.          0.03601961  0.          0.
      0.00745848  0.          0.          0.03601961  0.          0.
      0.09974018  0.00335702  0.01254176  0.00335702  0.01254176  0.
      0.01254176  0.10719866  0.          0.          0.          0.31026953
      0.          0.31026953  0.          0.01254176  0.          0.          0.
      0.02000024  0.0108155   0.03601961  0.          0.03601961  0.03937663
      0.00335702  0.01254176  0.          0.          0.          0.01688546
      0.          0.          0.06036354  0.03601961  0.          0.        ]

    References
    ----------
    .. [boost-kolmogorov] http://www.boost.org/libs/graph/doc/kolmogorov_max_flow.html
    .. [kolmogorov-graph-2003] Vladimir Kolmogorov, "Graph Based Algorithms for
       Scene Reconstruction from Two or More Views", PhD thesis, Cornell
        University, September 2003.
    .. [boykov-experimental-2004] Yuri Boykov and Vladimir Kolmogorov, "An
       Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy
       Minimization", Vision In IEEE Transactions on Pattern Analysis and
       Machine Intelligence, vol. 26, no. 9, pp. 1124-1137, Sept. 2004.
    """
    _check_prop_scalar(capacity, "capacity")
    if residual == None:
        residual = g.new_edge_property(capacity.value_type())
    _check_prop_scalar(residual, "residual")
    _check_prop_writable(residual, "residual")

    libgraph_tool_flow.\
           kolmogorov_max_flow(g._Graph__graph, int(source), int(target),
                                 _prop("e", g, capacity),
                                 _prop("e", g, residual))
    return residual


def max_cardinality_matching(g, match=None):
    r"""Find the maximum cardinality matching in the graph.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    match : :class:`~graph_tool.PropertyMap` (optional, default: none)
        Edge property map where the matching will be specified.

    Returns
    -------
    match : :class:`~graph_tool.PropertyMap`
        Boolean edge property map where the matching is specified.
    is_maximal : bool
        True if the matching is indeed maximal, or False otherwise.

    Notes
    -----
    A *matching* is a subset of the edges of a graph such that no two edges
    share a common vertex. A *maximum cardinality matching* has maximum size
    over all matchings in the graph.

    For a more detailed description, see [boost-max-matching]_.

    Examples
    --------
    >>> from numpy.random import seed, random
    >>> seed(43)
    >>> g = gt.random_graph(100, lambda: (2,2))
    >>> res = gt.max_cardinality_matching(g)
    >>> print res[1]
    True
    >>> gt.graph_draw(g, ecolor=res[0], output="max_card_match.png")
    <...>

    .. figure:: max_card_match.png
        :align: center

        Edges belonging to the matching are in red.

    References
    ----------
    .. [boost-max-matching] http://www.boost.org/libs/graph/doc/maximum_matching.html
    """
    if match == None:
        match = g.new_edge_property("bool")
    _check_prop_scalar(match, "match")
    _check_prop_writable(match, "match")

    g.stash_filter(directed=True)
    g.set_directed(False)
    check = libgraph_tool_flow.\
            max_cardinality_matching(g._Graph__graph, _prop("e", g, match))
    g.pop_filter(directed=True)
    return match, check
