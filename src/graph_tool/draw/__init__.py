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
``graph_tool.draw`` - Graph drawing and layout
----------------------------------------------

Summary
+++++++

Layout algorithms
=================

.. autosummary::
   :nosignatures:

   sfdp_layout
   fruchterman_reingold_layout
   arf_layout
   radial_tree_layout
   random_layout
   get_hierarchy_control_points

Graph drawing
=============

.. autosummary::
   :nosignatures:

   graph_draw
   graphviz_draw
   prop_to_size


Low-level graph drawing
^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :nosignatures:

   cairo_draw
   interactive_window
   GraphWidget
   GraphWindow

Contents
++++++++
"""

from __future__ import division, absolute_import, print_function

from .. import GraphView, _check_prop_vector, group_vector_property, \
     ungroup_vector_property, infect_vertex_property, _prop, _get_rng
from .. topology import max_cardinality_matching, max_independent_vertex_set, \
    label_components, pseudo_diameter, shortest_distance
from .. community import condensation_graph
from .. stats import label_parallel_edges
from .. generation import predecessor_tree
import numpy.random
from numpy import sqrt
import sys

from .. dl_import import dl_import
dl_import("from . import libgraph_tool_layout")


__all__ = ["graph_draw", "graphviz_draw", "fruchterman_reingold_layout",
           "arf_layout", "sfdp_layout", "random_layout", "radial_tree_layout",
           "cairo_draw", "prop_to_size", "get_hierarchy_control_points",
           "default_cm"]


def random_layout(g, shape=None, pos=None, dim=2):
    r"""Performs a random layout of the graph.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    shape : tuple or list (optional, default: ``None``)
        Rectangular shape of the bounding area. The size of this parameter must
        match `dim`, and each element can be either a pair specifying a range,
        or a single value specifying a range starting from zero. If None is
        passed, a square of linear size :math:`\sqrt{N}` is used.
    pos : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Vector vertex property maps where the coordinates should be stored.
    dim : int (optional, default: ``2``)
        Number of coordinates per vertex.

    Returns
    -------
    pos : :class:`~graph_tool.PropertyMap`
        A vector-valued vertex property map with the coordinates of the
        vertices.

    Notes
    -----
    This algorithm has complexity :math:`O(V)`.

    Examples
    --------
    .. testcode::
       :hide:

       np.random.seed(42)
       gt.seed_rng(42)

    >>> g = gt.random_graph(100, lambda: (3, 3))
    >>> shape = [[50, 100], [1, 2], 4]
    >>> pos = gt.random_layout(g, shape=shape, dim=3)
    >>> pos[g.vertex(0)].a
    array([ 68.72700594,   1.03142919,   2.56812658])

    """

    if pos == None:
        pos = g.new_vertex_property("vector<double>")
    _check_prop_vector(pos, name="pos")

    pos = ungroup_vector_property(pos, list(range(0, dim)))

    if shape == None:
        shape = [sqrt(g.num_vertices())] * dim

    for i in range(dim):
        if hasattr(shape[i], "__len__"):
            if len(shape[i]) != 2:
                raise ValueError("The elements of 'shape' must have size 2.")
            r = [min(shape[i]), max(shape[i])]
        else:
            r = [min(shape[i], 0), max(shape[i], 0)]
        d = r[1] - r[0]

        # deal with filtering
        p = pos[i].ma
        p[:] = numpy.random.random(len(p)) * d + r[0]

    pos = group_vector_property(pos)
    return pos


def fruchterman_reingold_layout(g, weight=None, a=None, r=1., scale=None,
                                circular=False, grid=True, t_range=None,
                                n_iter=100, pos=None):
    r"""Calculate the Fruchterman-Reingold spring-block layout of the graph.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    weight : :class:`PropertyMap` (optional, default: ``None``)
        An edge property map with the respective weights.
    a : float (optional, default: :math:`V`)
        Attracting force between adjacent vertices.
    r : float (optional, default: 1.0)
        Repulsive force between vertices.
    scale : float (optional, default: :math:`\sqrt{V}`)
        Total scale of the layout (either square side or radius).
    circular : bool (optional, default: ``False``)
        If ``True``, the layout will have a circular shape. Otherwise the shape
        will be a square.
    grid : bool (optional, default: ``True``)
        If ``True``, the repulsive forces will only act on vertices which are on
        the same site on a grid. Otherwise they will act on all vertex pairs.
    t_range : tuple of floats (optional, default: ``(scale / 10, scale / 1000)``)
        Temperature range used in annealing. The temperature limits the
        displacement at each iteration.
    n_iter : int (optional, default: ``100``)
        Total number of iterations.
    pos : :class:`PropertyMap` (optional, default: ``None``)
        Vector vertex property maps where the coordinates should be stored. If
        provided, this will also be used as the initial position of the
        vertices.

    Returns
    -------
    pos : :class:`~graph_tool.PropertyMap`
        A vector-valued vertex property map with the coordinates of the
        vertices.

    Notes
    -----
    This algorithm is defined in [fruchterman-reingold]_, and has
    complexity :math:`O(\text{n-iter}\times V^2)` if `grid=False` or
    :math:`O(\text{n-iter}\times (V + E))` otherwise.

    Examples
    --------
    .. testcode::
       :hide:

       np.random.seed(42)
       gt.seed_rng(42)

    >>> g = gt.price_network(300)
    >>> pos = gt.fruchterman_reingold_layout(g, n_iter=1000)
    >>> gt.graph_draw(g, pos=pos, output="graph-draw-fr.pdf")
    <...>

    .. testcode::
       :hide:

       gt.graph_draw(g, pos=pos, output="graph-draw-fr.png")

    .. figure:: graph-draw-fr.*
        :align: center

        Fruchterman-Reingold layout of a Price network.

    References
    ----------
    .. [fruchterman-reingold] Fruchterman, Thomas M. J.; Reingold, Edward M.
       "Graph Drawing by Force-Directed Placement". Software - Practice & Experience
       (Wiley) 21 (11): 1129-1164. (1991) :doi:`10.1002/spe.4380211102`
    """

    if pos == None:
        pos = random_layout(g, dim=2)
    _check_prop_vector(pos, name="pos", floating=True)

    if a is None:
        a = float(g.num_vertices())

    if scale is None:
        scale = sqrt(g.num_vertices())

    if t_range is None:
        t_range = (scale / 10, scale / 1000)

    ug = GraphView(g, directed=False)
    libgraph_tool_layout.fruchterman_reingold_layout(ug._Graph__graph,
                                                     _prop("v", g, pos),
                                                     _prop("e", g, weight),
                                                     a, r, not circular, scale,
                                                     grid, t_range[0],
                                                     t_range[1], n_iter)
    return pos


def arf_layout(g, weight=None, d=0.5, a=10, dt=0.001, epsilon=1e-6,
               max_iter=1000, pos=None, dim=2):
    r"""Calculate the ARF spring-block layout of the graph.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    weight : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        An edge property map with the respective weights.
    d : float (optional, default: ``0.5``)
        Opposing force between vertices.
    a : float (optional, default: ``10``)
        Attracting force between adjacent vertices.
    dt : float (optional, default: ``0.001``)
        Iteration step size.
    epsilon : float (optional, default: ``1e-6``)
        Convergence criterion.
    max_iter : int (optional, default: ``1000``)
        Maximum number of iterations. If this value is ``0``, it runs until
        convergence.
    pos : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Vector vertex property maps where the coordinates should be stored.
    dim : int (optional, default: ``2``)
        Number of coordinates per vertex.

    Returns
    -------
    pos : :class:`~graph_tool.PropertyMap`
        A vector-valued vertex property map with the coordinates of the
        vertices.

    Notes
    -----
    This algorithm is defined in [geipel-self-organization-2007]_, and has
    complexity :math:`O(V^2)`.

    Examples
    --------
    .. testcode::
       :hide:

       np.random.seed(42)
       gt.seed_rng(42)

    >>> g = gt.price_network(300)
    >>> pos = gt.arf_layout(g, max_iter=0)
    >>> gt.graph_draw(g, pos=pos, output="graph-draw-arf.pdf")
    <...>

    .. testcode::
       :hide:

       gt.graph_draw(g, pos=pos, output="graph-draw-arf.png")

    .. figure:: graph-draw-arf.*
        :align: center

        ARF layout of a Price network.

    References
    ----------
    .. [geipel-self-organization-2007] Markus M. Geipel, "Self-Organization
       applied to Dynamic Network Layout", International Journal of Modern
       Physics C vol. 18, no. 10 (2007), pp. 1537-1549,
       :doi:`10.1142/S0129183107011558`, :arxiv:`0704.1748v5`
    .. _arf: http://www.sg.ethz.ch/research/graphlayout
    """

    if pos is None:
        pos = random_layout(g, dim=dim)
    _check_prop_vector(pos, name="pos", floating=True)

    ug = GraphView(g, directed=False)
    libgraph_tool_layout.arf_layout(ug._Graph__graph, _prop("v", g, pos),
                                    _prop("e", g, weight), d, a, dt, max_iter,
                                    epsilon, dim)
    return pos


def _coarse_graph(g, vweight, eweight, mivs=False, groups=None):
    if groups is None:
        if mivs:
            mivs = max_independent_vertex_set(g, high_deg=True)
            u = GraphView(g, vfilt=mivs, directed=False)
            c = label_components(u)[0]
            c.fa += 1
            u = GraphView(g, directed=False)
            infect_vertex_property(u, c,
                                   list(range(1, c.fa.max() + 1)))
            c = g.own_property(c)
        else:
            mivs = None
            m = max_cardinality_matching(GraphView(g, directed=False),
                                         heuristic=True, weight=eweight,
                                         minimize=False)
            u = GraphView(g, efilt=m, directed=False)
            c = label_components(u)[0]
            c = g.own_property(c)
            u = GraphView(g, directed=False)
    else:
        mivs = None
        c = groups
    cg, cc, vcount, ecount = condensation_graph(g, c, vweight, eweight)[:4]
    return cg, cc, vcount, ecount, c, mivs


def _propagate_pos(g, cg, c, cc, cpos, delta, mivs):
    pos = g.new_vertex_property(cpos.value_type())

    if mivs is not None:
        g = GraphView(g, vfilt=mivs)
    libgraph_tool_layout.propagate_pos(g._Graph__graph,
                                       cg._Graph__graph,
                                       _prop("v", g, c),
                                       _prop("v", cg, cc),
                                       _prop("v", g, pos),
                                       _prop("v", cg, cpos),
                                       delta if mivs is None else 0,
                                       _get_rng())

    if mivs is not None:
        g = g.base
        u = GraphView(g, directed=False)
        try:
            libgraph_tool_layout.propagate_pos_mivs(u._Graph__graph,
                                                    _prop("v", u, mivs),
                                                    _prop("v", u, pos),
                                                    delta, _get_rng())
        except ValueError:
            graph_draw(u, mivs, vertex_fillcolor=mivs)
    return pos


def _avg_edge_distance(g, pos):
    libgraph_tool_layout.sanitize_pos(g._Graph__graph, _prop("v", g, pos))
    ad = libgraph_tool_layout.avg_dist(g._Graph__graph, _prop("v", g, pos))
    if numpy.isnan(ad) or ad == 0:
        ad = 1.
    return ad


def coarse_graphs(g, method="hybrid", mivs_thres=0.9, ec_thres=0.75,
                  weighted_coarse=False, eweight=None, vweight=None,
                  groups=None, verbose=False):
    cg = [[g, None, None, None, None, None]]
    if weighted_coarse:
        cg[-1][2], cg[-1][3] = vweight, eweight
    mivs = not (method in ["hybrid", "ec"])
    while True:
        u = _coarse_graph(cg[-1][0], cg[-1][2], cg[-1][3], mivs, groups)
        groups = None
        thres = mivs_thres if mivs else ec_thres
        if u[0].num_vertices() >= thres * cg[-1][0].num_vertices():
            if method == "hybrid" and not mivs:
                mivs = True
            else:
                break
        if u[0].num_vertices() <= 2:
            break
        cg.append(u)
        if verbose:
            print("Coarse level (%s):" % ("MIVS" if mivs else "EC"), end=' ')
            print(len(cg), " num vertices:", end=' ')
            print(u[0].num_vertices())
    cg.reverse()
    Ks = []
    pos = random_layout(cg[0][0], dim=2)
    for i in range(len(cg)):
        if i == 0:
            u = cg[i][0]
            K = _avg_edge_distance(u, pos)
            if K == 0:
                K = 1.
            Ks.append(K)
            continue
        if weighted_coarse:
            gamma = 1.
        else:
            #u = cg[i - 1][0]
            #w = cg[i][0]
            #du = pseudo_diameter(u)[0]
            #dw = pseudo_diameter(w)[0]
            #gamma = du / float(max(dw, du))
            gamma = 0.75
        Ks.append(Ks[-1] * gamma)

    for i in range(len(cg)):
        u, cc, vcount, ecount, c, mivs = cg[i]
        yield u, pos, Ks[i], vcount, ecount

        if verbose:
            print("avg edge distance:", _avg_edge_distance(u, pos))

        if i < len(cg) - 1:
            if verbose:
                print("propagating...", end=' ')
                print(mivs.a.sum() if mivs is not None else "")
            pos = _propagate_pos(cg[i + 1][0], u, c, cc, pos,
                                 Ks[i] / 1000., mivs)

def coarse_graph_stack(g, c, coarse_stack, eweight=None, vweight=None,
                       weighted_coarse=True, verbose=False):
    cg = [[g, c, None, None]]
    if weighted_coarse:
        cg[-1][2], cg[-1][3] = vweight, eweight
    for u in coarse_stack:
        c = u.vp["b"]
        vcount = u.vp["count"]
        ecount = u.ep["count"]
        cg.append((u, c, vcount, ecount))
        if verbose:
            print("Coarse level:", end=' ')
            print(len(cg), " num vertices:", end=' ')
            print(u.num_vertices())
    cg.reverse()
    Ks = []
    pos = random_layout(cg[0][0], dim=2)
    for i in range(len(cg)):
        if i == 0:
            u = cg[i][0]
            K = _avg_edge_distance(u, pos)
            if K == 0:
                K = 1.
            Ks.append(K)
            continue
        if weighted_coarse:
            gamma = 1.
        else:
            #u = cg[i - 1][0]
            #w = cg[i][0]
            #du = pseudo_diameter(u)[0]
            #dw = pseudo_diameter(w)[0]
            #gamma = du / float(max(dw, du))
            gamma = 0.75
        Ks.append(Ks[-1] * gamma)

    for i in range(len(cg)):
        u, c, vcount, ecount = cg[i]
        yield u, pos, Ks[i], vcount, ecount

        if verbose:
            print("avg edge distance:", _avg_edge_distance(u, pos))

        if i < len(cg) - 1:
            if verbose:
                print("propagating...")
            pos = _propagate_pos(cg[i + 1][0], u, c, u.vertex_index.copy("int"),
                                 pos, Ks[i] / 1000., None)


def sfdp_layout(g, vweight=None, eweight=None, pin=None, groups=None, C=0.2,
                K=None, p=2., theta=0.6, max_level=11, gamma=1., mu=0., mu_p=1.,
                init_step=None, cooling_step=0.95, adaptive_cooling=True,
                epsilon=1e-2, max_iter=0, pos=None, multilevel=None,
                coarse_method="hybrid", mivs_thres=0.9, ec_thres=0.75,
                coarse_stack=None, weighted_coarse=False, verbose=False):
    r"""Obtain the SFDP spring-block layout of the graph.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    vweight : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        A vertex property map with the respective weights.
    eweight : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        An edge property map with the respective weights.
    pin : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        A vertex property map with boolean values, which, if given,
        specify the vertices which will not have their positions modified.
    groups : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        A vertex property map with group assignments. Vertices belonging to the
        same group will be put close together.
    C : float (optional, default: ``0.2``)
        Relative strength of repulsive forces.
    K : float (optional, default: ``None``)
        Optimal edge length. If not provided, it will be taken to be the average
        edge distance in the initial layout.
    p : float (optional, default: ``2``)
        Repulsive force exponent.
    theta : float (optional, default: ``0.6``)
        Quadtree opening parameter, a.k.a. Barnes-Hut opening criterion.
    max_level : int (optional, default: ``11``)
        Maximum quadtree level.
    gamma : float (optional, default: ``1.0``)
        Strength of the attractive force between connected components, or group
        assignments.
    mu : float (optional, default: ``0.0``)
        Strength of the attractive force between vertices of the same connected
        component, or group assignment.
    mu_p : float (optional, default: ``1.0``)
        Scaling exponent of the attractive force between vertices of the same
        connected component, or group assignment.
    init_step : float (optional, default: ``None``)
        Initial update step. If not provided, it will be chosen automatically.
    cooling_step : float (optional, default: ``0.95``)
        Cooling update step.
    adaptive_cooling : bool (optional, default: ``True``)
        Use an adaptive cooling scheme.
    epsilon : float (optional, default: ``0.01``)
        Relative convergence criterion.
    max_iter : int (optional, default: ``0``)
        Maximum number of iterations. If this value is ``0``, it runs until
        convergence.
    pos : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Initial vertex layout. If not provided, it will be randomly chosen.
    multilevel : bool (optional, default: ``None``)
        Use a multilevel layout algorithm. If ``None`` is given, it will be
        activated based on the size of the graph.
    coarse_method : str (optional, default: ``"hybrid"``)
        Coarsening method used if ``multilevel == True``. Allowed methods are
        ``"hybrid"``, ``"mivs"`` and ``"ec"``.
    mivs_thres : float (optional, default: ``0.9``)
        If the relative size of the MIVS coarse graph is above this value, the
        coarsening stops.
    ec_thres : float (optional, default: ``0.75``)
        If the relative size of the EC coarse graph is above this value, the
        coarsening stops.
    weighted_coarse : bool (optional, default: ``False``)
        Use weighted coarse graphs.
    verbose : bool (optional, default: ``False``)
        Provide verbose information.

    Returns
    -------
    pos : :class:`~graph_tool.PropertyMap`
        A vector-valued vertex property map with the coordinates of the
        vertices.

    Notes
    -----
    This algorithm is defined in [hu-multilevel-2005]_, and has
    complexity :math:`O(V\log V)`.

    Examples
    --------
    .. testcode::
       :hide:

       np.random.seed(42)
       gt.seed_rng(42)

    >>> g = gt.price_network(3000)
    >>> pos = gt.sfdp_layout(g)
    >>> gt.graph_draw(g, pos=pos, output="graph-draw-sfdp.pdf")
    <...>

    .. testcode::
       :hide:

       gt.graph_draw(g, pos=pos, output="graph-draw-sfdp.png")

    .. figure:: graph-draw-sfdp.*
        :align: center

        SFDP layout of a Price network.

    References
    ----------
    .. [hu-multilevel-2005] Yifan Hu, "Efficient and High Quality Force-Directed
       Graph", Mathematica Journal, vol. 10, Issue 1, pp. 37-71, (2005)
       http://www.mathematica-journal.com/issue/v10i1/graph_draw.html
    """

    if pos is None:
        pos = random_layout(g, dim=2)
    _check_prop_vector(pos, name="pos", floating=True)

    g = GraphView(g, directed=False)

    if pin is not None:
        if pin.value_type() != "bool":
            raise ValueError("'pin' property must be of type 'bool'.")
    else:
        pin = g.new_vertex_property("bool")

    if K is None:
        K = _avg_edge_distance(g, pos)

    if init_step is None:
        init_step = 2 * max(_avg_edge_distance(g, pos), K)

    if multilevel is None:
        multilevel = g.num_vertices() > 1000

    if multilevel:
        if eweight is not None or vweight is not None:
            weighted_coarse = True
        if coarse_stack is None:
            cgs = coarse_graphs(g, method=coarse_method,
                                mivs_thres=mivs_thres,
                                ec_thres=ec_thres,
                                weighted_coarse=weighted_coarse,
                                eweight=eweight,
                                vweight=vweight,
                                groups=groups,
                                verbose=verbose)
        else:
            cgs = coarse_graph_stack(g, coarse_stack[0], coarse_stack[1],
                                     eweight=eweight, vweight=vweight,
                                     verbose=verbose)
        for count, (u, pos, K, vcount, ecount) in enumerate(cgs):
            if verbose:
                print("Positioning level:", count, u.num_vertices(), end=' ')
                print("with K =", K, "...")
                count += 1
            #graph_draw(u, pos)
            pos = sfdp_layout(u, pos=pos,
                              vweight=vcount if weighted_coarse else None,
                              eweight=ecount if weighted_coarse else None,
                              groups=None if u.num_vertices() < g.num_vertices() else groups,
                              C=C, K=K, p=p,
                              theta=theta, gamma=gamma, mu=mu, mu_p=mu_p,
                              epsilon=epsilon,
                              max_iter=max_iter,
                              cooling_step=cooling_step,
                              adaptive_cooling=False,
                              # init_step=max(2 * K,
                              #               _avg_edge_distance(u, pos)),
                              multilevel=False,
                              verbose=False)
            #graph_draw(u, pos)
        return pos

    if g.num_vertices() <= 1:
        return pos
    if g.num_vertices() == 2:
        vs = [g.vertex(0, False), g.vertex(1, False)]
        pos[vs[0]] = [0, 0]
        pos[vs[1]] = [1, 1]
        return pos
    if g.num_vertices() <= 50:
        max_level = 0
    if groups is None:
        groups = label_components(g)[0]
    elif groups.value_type() != "int32_t":
        raise ValueError("'groups' property must be of type 'int32_t'.")
    libgraph_tool_layout.sanitize_pos(g._Graph__graph, _prop("v", g, pos))
    libgraph_tool_layout.sfdp_layout(g._Graph__graph, _prop("v", g, pos),
                                     _prop("v", g, vweight),
                                     _prop("e", g, eweight),
                                     _prop("v", g, pin),
                                     (C, K, p, gamma, mu, mu_p, _prop("v", g, groups)),
                                     theta, init_step, cooling_step, max_level,
                                     epsilon, max_iter, not adaptive_cooling,
                                     verbose, _get_rng())
    return pos

def radial_tree_layout(g, root, rel_order=None, weighted=False, r=1.):
    r"""Computes a radial layout of the graph according to the minimum spanning
    tree centered at the ``root`` vertex.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    root : :class:`~graph_tool.Vertex` or ``int``
        The root of the radial tree.
    rel_order : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Relative order of the nodes at the lowest branch.
    weighted : ``bool`` (optional, default: ``False``)
        If true, the angle between the child branches will be computed according
        to weight of the entire sub-branches.
    r : ``float`` (optional, default: ``1.``)
        Layer spacing.

    Returns
    -------
    pos : :class:`~graph_tool.PropertyMap`
        A vector-valued vertex property map with the coordinates of the
        vertices.

    Notes
    -----
    This algorithm has complexity :math:`O(V + E)`.

    Examples
    --------
    .. testcode::
       :hide:

       np.random.seed(42)
       gt.seed_rng(42)

    >>> g = gt.price_network(1000)
    >>> pos = gt.radial_tree_layout(g, g.vertex(0))
    >>> gt.graph_draw(g, pos=pos, output="graph-draw-radial.pdf")
    <...>

    .. testcode::
       :hide:

       gt.graph_draw(g, pos=pos, output="graph-draw-radial.png")

    .. figure:: graph-draw-radial.*
        :align: center

        Radial tree layout of a Price network.

    """

    levels, pred_map = shortest_distance(GraphView(g, directed=False), root,
                                         pred_map=True)
    t = predecessor_tree(g, pred_map)
    pos = t.new_vertex_property("vector<double>")
    levels = t.own_property(levels)
    if rel_order is None:
        rel_order = g.vertex_index.copy("int")

    libgraph_tool_layout.get_radial(t._Graph__graph,
                                    _prop("v", g, pos),
                                    _prop("v", g, levels),
                                    _prop("v", g, rel_order),
                                    int(root), weighted, r)
    return g.own_property(pos)

try:
    from .cairo_draw import graph_draw, cairo_draw, get_hierarchy_control_points, default_cm
except ImportError:
    pass

try:
    from .cairo_draw import GraphWidget, GraphWindow, \
        interactive_window
    __all__ += ["interactive_window", "GraphWidget", "GraphWindow"]
except ImportError:
    pass

try:
   from .graphviz_draw import graphviz_draw
except ImportError:
   pass

def prop_to_size(prop, mi=0, ma=5, log=False, power=0.5):
    r"""Convert property map values to be more useful as a vertex size, or edge
    width. The new values are taken to be

    .. math::

        y = mi + (ma - mi) \left(\frac{x_i - min(x)} {max(x) - min(x)}\right)^\text{power}

    If `log=True`, the natural logarithm of the property values is used instead.

    """
    prop = prop.copy(value_type="double")
    if log:
        vals = numpy.log(prop.fa)
    else:
        vals = prop.fa

    delta = vals.max() - vals.min()
    if delta == 0:
        delta = 1
    prop.fa = mi + (ma - mi) * ((vals - vals.min()) / delta) ** power
    return prop
