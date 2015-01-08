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
``graph_tool.stats`` - Miscellaneous statistics
-----------------------------------------------

Summary
+++++++

.. autosummary::
   :nosignatures:

   vertex_hist
   edge_hist
   vertex_average
   edge_average
   label_parallel_edges
   remove_parallel_edges
   label_self_loops
   remove_self_loops
   remove_labeled_edges
   distance_histogram

Contents
++++++++

"""

from __future__ import division, absolute_import, print_function

from .. dl_import import dl_import
dl_import("from . import libgraph_tool_stats")

from .. import _degree, _prop, _get_rng, GraphView
from numpy import *
import numpy
import sys

__all__ = ["vertex_hist", "edge_hist", "vertex_average", "edge_average",
           "label_parallel_edges", "remove_parallel_edges",
           "label_self_loops", "remove_self_loops", "remove_labeled_edges",
           "distance_histogram"]


def vertex_hist(g, deg, bins=[0, 1], float_count=True):
    """
    Return the vertex histogram of the given degree type or property.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    deg : string or :class:`~graph_tool.PropertyMap`
        Degree or property to be used for the histogram. It can be either "in",
        "out" or "total", for in-, out-, or total degree of the vertices. It can
        also be a vertex property map.
    bins : list of bins (optional, default: [0, 1])
        List of bins to be used for the histogram. The values given represent
        the edges of the bins (i.e. lower and upper bounds). If the list
        contains two values, this will be used to automatically create an
        appropriate bin range, with a constant width given by the second value,
        and starting from the first value.
    float_count : bool (optional, default: True)
        If True, the counts in each histogram bin will be returned as floats. If
        False, they will be returned as integers.

    Returns
    -------
    counts : :class:`~numpy.ndarray`
        The bin counts.
    bins : :class:`~numpy.ndarray`
        The bin edges.

    See Also
    --------
    edge_hist: Edge histograms.
    vertex_average: Average of vertex properties, degrees.
    edge_average: Average of edge properties.
    distance_histogram : Shortest-distance histogram.

    Notes
    -----
    The algorithm runs in :math:`O(|V|)` time.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    .. testsetup::

       import numpy.random
       numpy.random.seed(42)
       gt.seed_rng(42)

    >>> from numpy.random import poisson
    >>> g = gt.random_graph(1000, lambda: (poisson(5), poisson(5)))
    >>> print(gt.vertex_hist(g, "out"))
    [array([   7.,   33.,   91.,  145.,  165.,  164.,  152.,  115.,   62.,
             29.,   28.,    6.,    1.,    1.,    0.,    1.]), array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16], dtype=uint64)]
    """

    ret = libgraph_tool_stats.\
          get_vertex_histogram(g._Graph__graph, _degree(g, deg),
                               [float(x) for x in bins])
    return [array(ret[0], dtype="float64") if float_count else ret[0], ret[1]]


def edge_hist(g, eprop, bins=[0, 1], float_count=True):
    """
    Return the edge histogram of the given property.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    eprop : :class:`~graph_tool.PropertyMap`
        Edge property to be used for the histogram.
    bins : list of bins (optional, default: [0, 1])
        List of bins to be used for the histogram. The values given represent
        the edges of the bins (i.e. lower and upper bounds). If the list
        contains two values, this will be used to automatically create an
        appropriate bin range, with a constant width given by the second value,
        and starting from the first value.
    float_count : bool (optional, default: True)
        If True, the counts in each histogram bin will be returned as floats. If
        False, they will be returned as integers.

    Returns
    -------
    counts : :class:`~numpy.ndarray`
        The bin counts.
    bins : :class:`~numpy.ndarray`
        The bin edges.

    See Also
    --------
    vertex_hist : Vertex histograms.
    vertex_average : Average of vertex properties, degrees.
    edge_average : Average of edge properties.
    distance_histogram : Shortest-distance histogram.

    Notes
    -----
    The algorithm runs in :math:`O(|E|)` time.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    .. testsetup::

       import numpy.random
       numpy.random.seed(42)
       gt.seed_rng(42)

    >>> from numpy import arange
    >>> from numpy.random import random
    >>> g = gt.random_graph(1000, lambda: (5, 5))
    >>> eprop = g.new_edge_property("double")
    >>> eprop.get_array()[:] = random(g.num_edges())
    >>> print(gt.edge_hist(g, eprop, linspace(0, 1, 11)))
    [array([ 501.,  441.,  478.,  480.,  506.,  494.,  507.,  535.,  499.,  559.]), array([ 0. ,  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1. ])]

    """

    ret = libgraph_tool_stats.\
          get_edge_histogram(g._Graph__graph, _prop("e", g, eprop),
                             [float(x) for x in bins])
    return [array(ret[0], dtype="float64") if float_count else ret[0], ret[1]]


def vertex_average(g, deg):
    """
    Return the average of the given degree or vertex property.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    deg : string or :class:`~graph_tool.PropertyMap`
        Degree or property to be used for the histogram. It can be either "in",
        "out" or "total", for in-, out-, or total degree of the vertices. It can
        also be a vertex property map.

    Returns
    -------
    average : float
        The average of the given degree or property.
    std : float
        The standard deviation of the average.

    See Also
    --------
    vertex_hist : Vertex histograms.
    edge_hist : Edge histograms.
    edge_average : Average of edge properties.
    distance_histogram : Shortest-distance histogram.

    Notes
    -----
    The algorithm runs in :math:`O(|V|)` time.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    .. testsetup::

       import numpy.random
       numpy.random.seed(42)
       gt.seed_rng(42)

    >>> from numpy.random import poisson
    >>> g = gt.random_graph(1000, lambda: (poisson(5), poisson(5)))
    >>> print(gt.vertex_average(g, "in"))
    (4.975, 0.0686758691244603)
    """

    ret = libgraph_tool_stats.\
          get_vertex_average(g._Graph__graph, _degree(g, deg))
    return ret


def edge_average(g, eprop):
    """
    Return the average of the given degree or vertex property.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    eprop : :class:`~graph_tool.PropertyMap`
        Edge property to be used for the histogram.

    Returns
    -------
    average : float
        The average of the given property.
    std : float
        The standard deviation of the average.

    See Also
    --------
    vertex_hist : Vertex histograms.
    edge_hist : Edge histograms.
    vertex_average : Average of vertex degree, properties.
    distance_histogram : Shortest-distance histogram.

    Notes
    -----
    The algorithm runs in :math:`O(|E|)` time.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    .. testsetup::

       import numpy.random
       numpy.random.seed(42)
       gt.seed_rng(42)

    >>> from numpy import arange
    >>> from numpy.random import random
    >>> g = gt.random_graph(1000, lambda: (5, 5))
    >>> eprop = g.new_edge_property("double")
    >>> eprop.get_array()[:] = random(g.num_edges())
    >>> print(gt.edge_average(g, eprop))
    (0.4989741369720412, 0.004101065927783255)
    """

    ret = libgraph_tool_stats.\
          get_edge_average(g._Graph__graph, _prop("e", g, eprop))
    return ret


def remove_labeled_edges(g, label):
    """Remove every edge `e` such that `label[e] != 0`."""
    u = GraphView(g, directed=True, reversed=g.is_reversed(),
                  skip_properties=True)
    libgraph_tool_stats.\
          remove_labeled_edges(u._Graph__graph, _prop("e", g, label))


def label_parallel_edges(g, mark_only=False, count_all=False, eprop=None):
    r"""Label edges which are parallel, i.e, have the same source and target
    vertices. For each parallel edge set :math:`PE`, the labelling starts from 0
    to :math:`|PE|-1`. (If `count_all==True`, the range is 0 to :math:`|PE|`
    instead). If `mark_only==True`, all parallel edges are simply marked with
    the value 1. If the `eprop` parameter is given
    (a :class:`~graph_tool.PropertyMap`), the labelling is stored there."""
    if eprop is None:
        if mark_only:
            eprop = g.new_edge_property("bool")
        else:
            eprop = g.new_edge_property("int32_t")
    libgraph_tool_stats.\
          label_parallel_edges(g._Graph__graph, _prop("e", g, eprop),
                               mark_only, count_all)
    return eprop


def remove_parallel_edges(g):
    """Remove all parallel edges from the graph. Only one edge from each
    parallel edge set is left."""
    eprop = label_parallel_edges(g)
    remove_labeled_edges(g, eprop)


def label_self_loops(g, mark_only=False, eprop=None):
    """Label edges which are self-loops, i.e, the source and target vertices are
    the same. For each self-loop edge set :math:`SL`, the labelling starts from 0
    to :math:`|SL|-1`. If `mark_only == True`, self-loops are labeled with 1
    and others with 0. If the `eprop` parameter is given
    (a :class:`~graph_tool.PropertyMap`), the labelling is stored there."""

    if eprop is None:
        if mark_only:
            eprop = g.new_edge_property("bool")
        else:
            eprop = g.new_edge_property("int32_t")
    libgraph_tool_stats.\
          label_self_loops(g._Graph__graph, _prop("e", g, eprop),
                           mark_only)
    return eprop


def remove_self_loops(g):
    """Remove all self-loops edges from the graph."""
    eprop = label_self_loops(g)
    remove_labeled_edges(g, eprop)


def distance_histogram(g, weight=None, bins=[0, 1], samples=None,
                       float_count=True):
    r"""
    Return the shortest-distance histogram for each vertex pair in the graph.

    Parameters
    ----------
    g : :class:`Graph`
        Graph to be used.
    weight : :class:`~graph_tool.PropertyMap` (optional, default: None)
        Edge weights.
    bins : list of bins (optional, default: [0, 1])
        List of bins to be used for the histogram. The values given represent
        the edges of the bins (i.e. lower and upper bounds). If the list
        contains two values, this will be used to automatically create an
        appropriate bin range, with a constant width given by the second value,
        and starting from the first value.
    samples : int (optional, default: None)
        If supplied, the distances will be randomly sampled from a number of
        source vertices given by this parameter. It `samples == None` (default),
        all pairs are used.
    float_count : bool (optional, default: True)
        If True, the counts in each histogram bin will be returned as floats. If
        False, they will be returned as integers.

    Returns
    -------
    counts : :class:`~numpy.ndarray`
        The bin counts.
    bins : :class:`~numpy.ndarray`
        The bin edges.

    See Also
    --------
    vertex_hist : Vertex histograms.
    edge_hist : Edge histograms.
    vertex_average : Average of vertex degree, properties.
    distance_histogram : Shortest-distance histogram.

    Notes
    -----
    The algorithm runs in :math:`O(V^2)` time, or :math:`O(V^2\log V)` if
    `weight != None`. If `samples` is supplied, the complexities are
    :math:`O(\text{samples}\times V)`  and
    :math:`O(\text{samples}\times V\log V)`, respectively.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    .. testsetup::

       import numpy.random
       numpy.random.seed(42)
       gt.seed_rng(42)

    >>> g = gt.random_graph(100, lambda: (3, 3))
    >>> hist = gt.distance_histogram(g)
    >>> print(hist)
    [array([    0.,   300.,   865.,  2214.,  3857.,  2480.,   184.]), array([0, 1, 2, 3, 4, 5, 6, 7], dtype=uint64)]
    >>> hist = gt.distance_histogram(g, samples=10)
    >>> print(hist)
    [array([   0.,   30.,   88.,  226.,  391.,  240.,   15.]), array([0, 1, 2, 3, 4, 5, 6, 7], dtype=uint64)]
    """

    if samples != None:
        ret = libgraph_tool_stats.\
              sampled_distance_histogram(g._Graph__graph,
                                         _prop("e", g, weight),
                                         [float(x) for x in bins],
                                         samples, _get_rng())
    else:
        ret = libgraph_tool_stats.\
              distance_histogram(g._Graph__graph, _prop("e", g, weight), bins)
    return [array(ret[0], dtype="float64") if float_count else ret[0], ret[1]]
