#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# graph_tool -- a general graph manipulation python module
#
# Copyright (C) 2007-2012 Tiago de Paula Peixoto <tiago@skewed.de>
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
   boykov_kolmogorov_max_flow

Contents
++++++++

The following network will be used as an example throughout the documentation.

.. testcode::

    from numpy.random import seed, random
    from scipy.linalg import norm
    seed(42)
    points = random((400, 2)) * 10
    points[0] = [0, 0]
    points[1] = [10, 10]
    g, pos = gt.triangulation(points, type="delaunay")
    g.set_directed(True)
    edges = list(g.edges())
    # reciprocate edges
    for e in edges:
        g.add_edge(e.target(), e.source())
    # The capacity will be defined as the inverse euclidian distance
    cap = g.new_edge_property("double")
    for e in g.edges():
        cap[e] = min(1.0 / norm(pos[e.target()].a - pos[e.source()].a), 10)
    g.edge_properties["cap"] = cap
    g.vertex_properties["pos"] = pos
    g.save("flow-example.xml.gz")
    cap.a /= cap.a.max() / 10
    gt.graph_draw(g, pos=pos, edge_pen_width=cap, output="flow-example.pdf")


.. figure:: flow-example.*
    :align: center

    Example network used in the documentation below. The edge capacities are
    represented by the edge width.

"""

from __future__ import division, absolute_import, print_function

from .. dl_import import dl_import
dl_import("from . import libgraph_tool_flow")

from .. import _prop, _check_prop_scalar, _check_prop_writable
__all__ = ["edmonds_karp_max_flow", "push_relabel_max_flow",
           "boykov_kolmogorov_max_flow"]


def edmonds_karp_max_flow(g, source, target, capacity, residual=None):
    r"""
    Calculate maximum flow on the graph with the Edmonds-Karp algorithm.

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
    is not as good as the push_relabel_max_flow() or the
    boykov_kolmogorov_max_flow() algorithm.

    - In the non-integer capacity case, the time complexity is :math:`O(V E^2)`
      which is worse than the time complexity of the push-relabel algorithm
      :math:`O(V^2E^{1/2})` for all but the sparsest of graphs.

    - In the integer capacity case, if the capacity bound U is very large then
      the algorithm will take a long time.

    Examples
    --------
    >>> g = gt.load_graph("flow-example.xml.gz")
    >>> cap = g.edge_properties["cap"]
    >>> src, tgt = g.vertex(0), g.vertex(1)
    >>> res = gt.edmonds_karp_max_flow(g, src, tgt, cap)
    >>> res.a = cap.a - res.a  # the actual flow
    >>> max_flow = sum(res[e] for e in tgt.in_edges())
    >>> print(max_flow)
    6.92759897841
    >>> pos = g.vertex_properties["pos"]
    >>> res.a /= res.a.max() / 10
    >>> gt.graph_draw(g, pos=pos, edge_pen_width=res, output="example-edmonds-karp.pdf")
    <...>

    .. figure:: example-edmonds-karp.*
        :align: center

        Edge flows obtained with the Edmonds-Karp algorithm. The source and
        target are on the lower left and upper right corners, respectively. The
        edge flows are represented by the edge width.


    References
    ----------
    .. [boost-edmonds-karp] http://www.boost.org/libs/graph/doc/edmonds_karp_max_flow.html
    .. [edmonds-theoretical-1972] Jack Edmonds and Richard M. Karp, "Theoretical
       improvements in the algorithmic efficiency for network flow problems.
       Journal of the ACM", 19:248-264, 1972 :doi:`10.1145/321694.321699`
    .. [ravindra-network-1993] Ravindra K. Ahuja and Thomas L. Magnanti and
       James B. Orlin,"Network Flows: Theory, Algorithms, and Applications".
       Prentice Hall, 1993.
    """

    _check_prop_scalar(capacity, "capacity")
    if residual == None:
        residual = g.new_edge_property(capacity.value_type())
    _check_prop_scalar(residual, "residual")
    _check_prop_writable(residual, "residual")

    if not g.is_directed():
        raise ValueError("The graph provided must be directed!")

    libgraph_tool_flow.\
               edmonds_karp_max_flow(g._Graph__graph, int(source), int(target),
                                     _prop("e", g, capacity),
                                     _prop("e", g, residual))
    return residual


def push_relabel_max_flow(g, source, target, capacity, residual=None):
    r"""
    Calculate maximum flow on the graph with the push-relabel algorithm.

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
    >>> g = gt.load_graph("flow-example.xml.gz")
    >>> cap = g.edge_properties["cap"]
    >>> src, tgt = g.vertex(0), g.vertex(1)
    >>> res = gt.push_relabel_max_flow(g, src, tgt, cap)
    >>> res.a = cap.a - res.a  # the actual flow
    >>> max_flow = sum(res[e] for e in tgt.in_edges())
    >>> print(max_flow)
    6.92759897841
    >>> pos = g.vertex_properties["pos"]
    >>> res.a /= res.a.max() / 10
    >>> gt.graph_draw(g, pos=pos, edge_pen_width=res, output="example-push-relabel.pdf")
    <...>

    .. figure:: example-push-relabel.*
        :align: center


        Edge flows obtained with the push-relabel algorithm. The source and
        target are on the lower left and upper right corners, respectively. The
        edge flows are represented by the edge width.

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

    if not g.is_directed():
        raise ValueError("The graph provided must be directed!")

    libgraph_tool_flow.\
               push_relabel_max_flow(g._Graph__graph, int(source), int(target),
                                     _prop("e", g, capacity),
                                     _prop("e", g, residual))

    return residual


def boykov_kolmogorov_max_flow(g, source, target, capacity, residual=None):
    r"""Calculate maximum flow on the graph with the Boykov-Kolmogorov algorithm.

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
    >>> g = gt.load_graph("flow-example.xml.gz")
    >>> cap = g.edge_properties["cap"]
    >>> src, tgt = g.vertex(0), g.vertex(1)
    >>> res = gt.boykov_kolmogorov_max_flow(g, src, tgt, cap)
    >>> res.a = cap.a - res.a  # the actual flow
    >>> max_flow = sum(res[e] for e in tgt.in_edges())
    >>> print(max_flow)
    6.92759897841
    >>> pos = g.vertex_properties["pos"]
    >>> res.a /= res.a.max() / 10
    >>> gt.graph_draw(g, pos=pos, edge_pen_width=res, output="example-kolmogorov.pdf")
    <...>

    .. figure:: example-kolmogorov.*
        :align: center

        Edge flows obtained with the Boykov-Kolmogorov algorithm. The source and
        target are on the lower left and upper right corners, respectively. The
        edge flows are represented by the edge width.

    References
    ----------
    .. [boost-kolmogorov] http://www.boost.org/libs/graph/doc/boykov_kolmogorov_max_flow.html
    .. [kolmogorov-graph-2003] Vladimir Kolmogorov, "Graph Based Algorithms for
       Scene Reconstruction from Two or More Views", PhD thesis, Cornell
       University, September 2003.
    .. [boykov-experimental-2004] Yuri Boykov and Vladimir Kolmogorov, "An
       Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy
       Minimization in Vision", IEEE Transactions on Pattern Analysis and
       Machine Intelligence, vol. 26, no. 9, pp. 1124-1137, Sept. 2004.
       :doi:`10.1109/TPAMI.2004.60`
    """
    _check_prop_scalar(capacity, "capacity")
    if residual == None:
        residual = g.new_edge_property(capacity.value_type())
    _check_prop_scalar(residual, "residual")
    _check_prop_writable(residual, "residual")

    if not g.is_directed():
        raise ValueError("The graph provided must be directed!")

    libgraph_tool_flow.\
               kolmogorov_max_flow(g._Graph__graph, int(source), int(target),
                                   _prop("e", g, capacity),
                                   _prop("e", g, residual))
    return residual

