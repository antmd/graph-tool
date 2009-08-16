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
``graph_tool.stats`` - Graph Statistics
---------------------------------------
"""

from .. dl_import import dl_import
dl_import("import libgraph_tool_stats")

from .. core import _degree, _prop
from numpy import *

__all__ = ["vertex_hist", "edge_hist", "vertex_average", "edge_average",
           "label_components", "label_biconnected_components",
           "label_parallel_edges", "remove_parallel_edges",
           "label_self_loops", "remove_self_loops", "remove_labeled_edges"]

def vertex_hist(g, deg, bins=[1], float_count=True):
    """
    Return the vertex histogram of the given degree type or property.

    Parameters
    ----------
    g : Graph
        Graph to be used.

    deg : string or PropertyMap
        Degree or property to be used for the histogram. It can be either "in",
        "out" or "total", for in-, out-, or total degree of the vertices. It can
        also be a vertex property map.

    bins : list of bins
        List of bins to be used for the histogram. The values given represent
        the edges of the bins (i,e, lower bounds). If the list contains only one
        value, this will be used to automatically create an appropriate bin
        range, with a constant lenght given by this value.

    float_count : bool (optional, default: True)
        If True, the counts in each histogram bin will be returned as floats. If
        False, they will be returned as integers.

    Returns
    -------
    counts : ndarray
        The bin counts.

    bins : ndarray
        The bin edges.

    See Also
    --------
    edge_hist: Edge histograms.
    vertex_average: Average of vertex properties, degrees.
    edge_average: Average of edge properties.

    Notes
    -----
    The algorithm runs in :math:`O(|V|)` time.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    >>> from numpy.random import poisson, seed
    >>> seed(42)
    >>> g = gt.random_graph(1000, lambda: (poisson(5), poisson(5)), seed=42)
    >>> print gt.vertex_hist(g, "out")
    [array([  10.,   35.,   92.,  140.,  162.,  160.,  150.,  109.,   66.,
             37.,   26.,   10.,    2.,    1.]), array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13], dtype=uint64)]

    """

    ret = libgraph_tool_stats.\
          get_vertex_histogram(g._Graph__graph, _degree(g, deg), bins)
    return [array(ret[0], dtype="float64") if float_count else ret[0], ret[1]]

def edge_hist(g, eprop, bins=[1], float_count=True):
    """
    Return the edge histogram of the given property.

    Parameters
    ----------
    g : Graph
        Graph to be used.

    eprop : PropertyMap
        Edge property to be used for the histogram.

    bins : list of bins
        List of bins to be used for the histogram. The values given represent
        the edges of the bins (i,e, lower bounds). If the list contains only one
        value, this will be used to automatically create an appropriate bin
        range, with a constant lenght given by this value.

    float_count : bool (optional, default: True)
        If True, the counts in each histogram bin will be returned as floats. If
        False, they will be returned as integers.

    Returns
    -------
    counts : ndarray
        The bin counts.

    bins : ndarray
        The bin edges.

    See Also
    --------
    vertex_hist : Vertex histograms.
    vertex_average : Average of vertex properties, degrees.
    edge_average : Average of edge properties.

    Notes
    -----
    The algorithm runs in :math:`O(|E|)` time.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    >>> from numpy import arange
    >>> from numpy.random import random
    >>> g = gt.random_graph(1000, lambda: (5, 5), seed=42)
    >>> eprop = g.new_edge_property("double")
    >>> eprop.get_array()[:] = random(g.num_edges())
    >>> print gt.edge_hist(g, eprop, arange(0, 1, 0.1))
    [array([ 500.,  440.,  487.,  494.,  507.,  496.,  524.,  526.,  486.,  540.]), array([ 0. ,  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9])]

    """

    ret = libgraph_tool_stats.\
          get_edge_histogram(g._Graph__graph, _prop("e", g, eprop), bins)
    return [array(ret[0], dtype="float64") if float_count else ret[0], ret[1]]

def vertex_average(g, deg):
    """
    Return the average of the given degree or vertex property.

    Parameters
    ----------
    g : Graph
        Graph to be used.

    deg : string or PropertyMap
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

    Notes
    -----
    The algorithm runs in :math:`O(|V|)` time.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    >>> from numpy.random import poisson, seed
    >>> seed(42)
    >>> g = gt.random_graph(1000, lambda: (poisson(5), poisson(5)), seed=42)
    >>> print gt.vertex_average(g, "in")
    (5.0179999999999998, 0.072661379012512559)
    """

    ret = libgraph_tool_stats.\
          get_vertex_average(g._Graph__graph, _degree(g, deg))
    return ret

def edge_average(g, eprop):
    """
    Return the average of the given degree or vertex property.

    Parameters
    ----------
    g : Graph
        Graph to be used.

    eprop : PropertyMap
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

    Notes
    -----
    The algorithm runs in :math:`O(|E|)` time.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    >>> from numpy import arange
    >>> from numpy.random import random
    >>> g = gt.random_graph(1000, lambda: (5, 5), seed=42)
    >>> eprop = g.new_edge_property("double")
    >>> eprop.get_array()[:] = random(g.num_edges())
    >>> print gt.edge_average(g, eprop)
    (0.50951471604395204, 0.0040790901147649975)
    """

    ret = libgraph_tool_stats.\
          get_edge_average(g._Graph__graph, _prop("e", g, eprop))
    return ret

def label_components(g, vprop=None, directed=None):
    """
    Labels the components to which each vertex in the graph belongs. If the
    graph is directed, it finds the strongly connected components.

    Parameters
    ----------
    g : Graph
        Graph to be used.

    vprop : PropertyMap (optional, default: None)
        Vertex property to store the component labels. If none is supplied, one
        is created.

    directed : bool (optional, default:None)
        Treat graph as directed or not, independently of its actual
        directionality.

    Returns
    -------
    comp : PropertyMap
        Vertex property map with component labels.

    Notes
    -----
    The components are arbitrarily labeled from 0 to N-1, where N is the total
    number of components.

    The algorithm runs in :math:`O(|V| + |E|)` time.

    Examples
    --------
    >>> g = gt.random_graph(100, lambda: (1, 1), seed=42)
    >>> comp = gt.label_components(g)
    >>> print comp.get_array()
    [0 1 2 3 4 0 3 3 4 4 2 3 4 0 3 3 3 3 0 3 2 1 3 0 0 2 2 3 3 3 0 1 2 3 2 3 0
     1 0 5 5 1 4 2 2 1 0 3 3 3 3 3 3 0 0 3 4 2 3 2 5 5 0 2 1 0 3 2 0 3 3 0 4 3
     2 6 2 2 1 3 1 1 0 3 0 1 3 0 3 0 2 0 2 2 0 6 1 1 0 2]
    """

    if directed != None:
        g.stash_filter(directed=True)
        g.set_directed(directed)

    if vprop == None:
        vprop = g.new_vertex_property("int32_t")
    libgraph_tool_stats.\
          label_components(g._Graph__graph, _prop("v", g, vprop))

    if directed != None:
        g.pop_filter(directed=True)
    return vprop

def label_biconnected_components(g, eprop=None, vprop=None):

    if vprop == None:
        vprop = g.new_vertex_property("bool")
    if eprop == None:
        eprop = g.new_edge_property("int32_t")

    g.stash_filter(directed=True)
    g.set_directed(False)
    nc = libgraph_tool_stats.\
          label_biconnected_components(g._Graph__graph, _prop("e", g, eprop),
                                       _prop("v", g, vprop))
    g.pop_filter(directed=True)
    return eprop, vprop, nc

def remove_labeled_edges(g, label):
    g.stash_filter(all=False, directed=True, reversed=True)
    libgraph_tool_stats.\
          remove_labeled_edges(g._Graph__graph, _prop("e", g, label))
    g.pop_filter(all=False, directed=True, reversed=True)

def label_parallel_edges(g, eprop=None):
    if eprop == None:
        eprop = g.new_edge_property("int32_t")
    libgraph_tool_stats.\
          label_parallel_edges(g._Graph__graph, _prop("e", g, eprop))
    return eprop

def remove_parallel_edges(g):
    eprop = label_parallel_edges(g)
    remove_labeled_edges(g, eprop)

def label_self_loops(g, eprop=None):
    if eprop == None:
        eprop = g.new_edge_property("int32_t")
    libgraph_tool_stats.\
          label_self_loops(g._Graph__graph, _prop("e", g, eprop))
    return eprop

def remove_self_loops(g):
    eprop = label_self_loops(g)
    remove_labeled_edges(g, eprop)

