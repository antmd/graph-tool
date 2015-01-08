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
``graph_tool.correlations`` - Correlations
------------------------------------------

Summary
+++++++

.. autosummary::
   :nosignatures:

   assortativity
   scalar_assortativity
   corr_hist
   combined_corr_hist
   avg_neighbour_corr
   avg_combined_corr

Contents
++++++++
"""

from __future__ import division, absolute_import, print_function

from .. dl_import import dl_import
dl_import("from . import libgraph_tool_correlations")

from .. import _degree, _prop
from numpy import *

__all__ = ["assortativity", "scalar_assortativity", "corr_hist",
           "combined_corr_hist", "avg_neighbour_corr", "avg_combined_corr"]


def assortativity(g, deg):
    r"""
    Obtain the assortativity coefficient for the given graph.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    deg : string or :class:`~graph_tool.PropertyMap`
        degree type ("in", "out" or "total") or vertex property map, which
        specifies the vertex types.

    Returns
    -------
    assortativity coefficient : tuple of two floats
        The assortativity coefficient, and its variance.

    See Also
    --------
    assortativity: assortativity coefficient
    scalar_assortativity: scalar assortativity coefficient
    corr_hist: vertex-vertex correlation histogram
    combined_corr_hist: combined single-vertex correlation histogram
    avg_neighbour_corr: average nearest-neighbour correlation
    avg_combined_corr: average combined single-vertex correlation

    Notes
    -----
    The assortativity coefficient [newman-mixing-2003]_ tells in a concise
    fashion how vertices of different types are preferentially connected amongst
    themselves, and is defined by

    .. math::
        r = \frac{\sum_i e_{ii} - \sum_i a_i b_i}{1-\sum_i a_i b_i}

    where :math:`a_i=\sum_je_{ij}` and :math:`b_j=\sum_ie_{ij}`, and
    :math:`e_{ij}` is the fraction of edges from a vertex of type
    i to a vertex of type j.


    The variance is obtained with the `jackknife method`_.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    .. testcode::
       :hide:

       np.random.seed(42)
       gt.seed_rng(42)
       from pylab import *

    >>> def sample_k(max):
    ...     accept = False
    ...     while not accept:
    ...         k = np.random.randint(1,max+1)
    ...         accept = random() < 1.0/k
    ...     return k
    ...
    >>> g = gt.random_graph(1000, lambda: sample_k(40), model="probabilistic",
    ...                     vertex_corr=lambda i,k: 1.0 / (1 + abs(i - k)), directed=False,
    ...                     n_iter=100)
    >>> gt.assortativity(g, "out")
    (0.1425854833491391, 0.005107908315851484)

    References
    ----------
    .. [newman-mixing-2003] M. E. J. Newman, "Mixing patterns in networks",
        Phys. Rev. E 67, 026126 (2003), :doi:`10.1103/PhysRevE.67.026126`
    .. _jackknife method: http://en.wikipedia.org/wiki/Resampling_%28statistics%29#Jackknife

    """
    return libgraph_tool_correlations.\
           assortativity_coefficient(g._Graph__graph, _degree(g, deg))


def scalar_assortativity(g, deg):
    r"""
    Obtain the scalar assortativity coefficient for the given graph.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    deg : string or :class:`~graph_tool.PropertyMap`
        degree type ("in", "out" or "total") or vertex property map, which
        specifies the vertex types.

    Returns
    -------
    scalar assortativity coefficient : tuple of two floats
        The scalar assortativity coefficient, and its variance.

    See Also
    --------
    assortativity: assortativity coefficient
    scalar_assortativity: scalar assortativity coefficient
    corr_hist: vertex-vertex correlation histogram
    combined_corr_hist: combined single-vertex correlation histogram
    avg_neighbour_corr: average nearest-neighbour correlation
    avg_combined_corr: average combined single-vertex correlation

    Notes
    -----
    The scalar assortativity coefficient [newman-mixing-2003]_ tells in a
    concise fashion how vertices of different types are preferentially connected
    amongst themselves, and is defined by

    .. math::
        r = \frac{\sum_{xy} xy(e_{xy} - a_x b_y)}{\sigma_a\sigma_b}

    where :math:`a_x=\sum_ye_{xy}` and :math:`b_y=\sum_xe_{xy}`, and
    :math:`e_{xy}` is the fraction of edges from a vertex of type
    x to a vertex of type y.

    The variance is obtained with the `jackknife method`_.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    .. testcode::
       :hide:

       np.random.seed(42)
       gt.seed_rng(42)
       from pylab import *

    >>> def sample_k(max):
    ...     accept = False
    ...     while not accept:
    ...         k = np.random.randint(1,max+1)
    ...         accept = random() < 1.0/k
    ...     return k
    ...
    >>> g = gt.random_graph(1000, lambda: sample_k(40), model="probabilistic",
    ...                     vertex_corr=lambda i,k: abs(i-k),
    ...                     directed=False, n_iter=100)
    >>> gt.scalar_assortativity(g, "out")
    (-0.44147938242732404, 0.01047064816533614)

    >>> g = gt.random_graph(1000, lambda: sample_k(40), model="probabilistic",
    ...                     vertex_corr=lambda i, k: 1.0 / (1 + abs(i - k)),
    ...                     directed=False, n_iter=100)
    >>> gt.scalar_assortativity(g, "out")
    (0.6216815584946155, 0.011152910837726114)

    References
    ----------
    .. [newman-mixing-2003] M. E. J. Newman, "Mixing patterns in networks",
        Phys. Rev. E 67, 026126 (2003), :doi:`10.1103/PhysRevE.67.026126`
    .. _jackknife method: http://en.wikipedia.org/wiki/Resampling_%28statistics%29#Jackknife
    """
    return libgraph_tool_correlations.\
           scalar_assortativity_coefficient(g._Graph__graph,
                                            _degree(g, deg))


def corr_hist(g, deg_source, deg_target, bins=[[0, 1], [0, 1]], weight=None,
              float_count=True):
    r"""
    Obtain the correlation histogram for the given graph.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    deg_source : string or :class:`~graph_tool.PropertyMap`
        degree type ("in", "out" or "total") or vertex property map for the
        source vertex.
    deg_target : string or :class:`~graph_tool.PropertyMap`
        degree type ("in", "out" or "total") or vertex property map for the
        target vertex.
    bins : list of lists (optional, default: [[0, 1], [0, 1]])
        A list of bin edges to be used for the source and target degrees. If any
        list has size 2, it is used to create an automatically generated bin
        range starting from the first value, and with constant bin width given
        by the second value.
    weight : edge property map (optional, default: None)
        Weight (multiplicative factor) to be used on each edge.
    float_count : bool (optional, default: True)
        If True, the bin counts are converted float variables, which is useful
        for normalization, and other processing. It False, the bin counts will
        be unsigned integers.

    Returns
    -------
    bin_counts : :class:`~numpy.ndarray`
        Two-dimensional array with the bin counts.
    source_bins : :class:`~numpy.ndarray`
        Source degree bins
    target_bins : :class:`~numpy.ndarray`
        Target degree bins


    See Also
    --------
    assortativity: assortativity coefficient
    scalar_assortativity: scalar assortativity coefficient
    corr_hist: vertex-vertex correlation histogram
    combined_corr_hist: combined single-vertex correlation histogram
    avg_neighbour_corr: average nearest-neighbour correlation
    avg_combined_corr: average combined single-vertex correlation

    Notes
    -----
    The correlation histogram counts, for every vertex with degree (or scalar
    property) 'source_deg', the number of out-neighbours with degree (or scalar
    property) 'target_deg'.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    .. testcode::
       :hide:

       np.random.seed(42)
       gt.seed_rng(42)
       from pylab import *

    >>> def sample_k(max):
    ...     accept = False
    ...     while not accept:
    ...         k = np.random.randint(1,max+1)
    ...         accept = random() < 1.0/k
    ...     return k
    ...
    >>> g = gt.random_graph(10000, lambda: sample_k(40), model="probabilistic",
    ...                     vertex_corr=lambda i, j: (sin(i / pi) * sin(j / pi) + 1) / 2,
    ...                     directed=False, n_iter=100)
    >>> h = gt.corr_hist(g, "out", "out")
    >>> clf()
    >>> xlabel("Source out-degree")
    <...>
    >>> ylabel("Target out-degree")
    <...>
    >>> imshow(h[0].T, interpolation="nearest", origin="lower")
    <...>
    >>> colorbar()
    <...>
    >>> savefig("corr.pdf")

    .. testcode::
       :hide:

       savefig("corr.png")

    .. figure:: corr.*
        :align: center

        Out/out-degree correlation histogram.
    """

    ret = libgraph_tool_correlations.\
          vertex_correlation_histogram(g._Graph__graph, _degree(g, deg_source),
                                       _degree(g, deg_target),
                                       _prop("e", g, weight),
                                       [float(x) for x in bins[0]],
                                       [float(x) for x in bins[1]])
    return [array(ret[0], dtype="float64") if float_count else ret[0],
            [ret[1][0], ret[1][1]]]


def combined_corr_hist(g, deg1, deg2, bins=[[0, 1], [0, 1]], float_count=True):
    r"""
    Obtain the single-vertex combined correlation histogram for the given graph.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    deg1 : string or :class:`~graph_tool.PropertyMap`
        first degree type ("in", "out" or "total") or vertex property map.
    deg2 : string or :class:`~graph_tool.PropertyMap`
        second degree type ("in", "out" or "total") or vertex property map.
    bins : list of lists (optional, default: [[0, 1], [0, 1]])
        A list of bin edges to be used for the first and second degrees. If any
        list has size 2, it is used to create an automatically generated bin
        range starting from the first value, and with constant bin width given
        by the second value.
    float_count : bool (optional, default: True)
        If True, the bin counts are converted float variables, which is useful
        for normalization, and other processing. It False, the bin counts will
        be unsigned integers.

    Returns
    -------
    bin_counts : :class:`~numpy.ndarray`
        Two-dimensional array with the bin counts.
    first_bins : :class:`~numpy.ndarray`
        First degree bins
    second_bins : :class:`~numpy.ndarray`
        Second degree bins

    Notes
    -----
    If enabled during compilation, this algorithm runs in parallel.

    See Also
    --------
    assortativity: assortativity coefficient
    scalar_assortativity: scalar assortativity coefficient
    corr_hist: vertex-vertex correlation histogram
    combined_corr_hist: combined single-vertex correlation histogram
    avg_neighbour_corr: average nearest-neighbour correlation
    avg_combined_corr: average combined single-vertex correlation

    Examples
    --------
    .. testcode::
       :hide:

       np.random.seed(42)
       gt.seed_rng(42)
       from pylab import *

    >>> def sample_k(max):
    ...     accept = False
    ...     while not accept:
    ...         i = randint(1, max + 1)
    ...         j = randint(1, max + 1)
    ...         accept = random() < (sin(i / pi) * sin(j / pi) + 1) / 2
    ...     return i,j
    ...
    >>> g = gt.random_graph(10000, lambda: sample_k(40))
    >>> h = gt.combined_corr_hist(g, "in", "out")
    >>> clf()
    >>> xlabel("In-degree")
    <...>
    >>> ylabel("Out-degree")
    <...>
    >>> imshow(h[0].T, interpolation="nearest", origin="lower")
    <...>
    >>> colorbar()
    <...>
    >>> savefig("combined_corr.pdf")

    .. testcode::
       :hide:

       savefig("combined_corr.png")

    .. figure:: combined_corr.*
        :align: center

        Combined in/out-degree correlation histogram.

    """
    ret = libgraph_tool_correlations.\
          vertex_combined_correlation_histogram(g._Graph__graph,
                                                _degree(g, deg1),
                                                _degree(g, deg2),
                                                [float(x) for x in bins[0]],
                                                [float(x) for x in bins[1]])
    return [array(ret[0], dtype="float64") if float_count else ret[0],
            [ret[1][0], ret[1][1]]]


def avg_neighbour_corr(g, deg_source, deg_target, bins=[0, 1], weight=None):
    r"""
    Obtain the average neighbour-neighbour correlation for the given graph.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    deg_source : string or :class:`~graph_tool.PropertyMap`
        degree type ("in", "out" or "total") or vertex property map for the
        source vertex.
    deg_target : string or :class:`~graph_tool.PropertyMap`
        degree type ("in", "out" or "total") or vertex property map for the
        target vertex.
    bins : list (optional, default: [0, 1])
        Bins to be used for the source degrees. If the list has size 2, it is
        used as the constant width of an automatically generated bin range,
        starting from the first value.
    weight : edge property map (optional, default: None)
        Weight (multiplicative factor) to be used on each edge.

    Returns
    -------
    bin_avg : :class:`~numpy.ndarray`
        Array with the deg_target average for the get_source bins.
    bin_dev : :class:`~numpy.ndarray`
        Array with the standard deviation of the deg_target average for the
        get_source bins.
    bins : :class:`~numpy.ndarray`
        Source degree bins,


    See Also
    --------
    assortativity: assortativity coefficient
    scalar_assortativity: scalar assortativity coefficient
    corr_hist: vertex-vertex correlation histogram
    combined_corr_hist: combined single-vertex correlation histogram
    avg_neighbour_corr: average nearest-neighbour correlation
    avg_combined_corr: average combined single-vertex correlation

    Notes
    -----
    The average correlation is the average, for every vertex with degree (or
    scalar property) 'source_deg', the of the 'target_deg' degree (or
    scalar property) of its neighbours.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    .. testcode::
       :hide:

       np.random.seed(42)
       gt.seed_rng(42)
       from pylab import *

    >>> def sample_k(max):
    ...     accept = False
    ...     while not accept:
    ...         k = randint(1,max+1)
    ...         accept = random() < 1.0 / k
    ...     return k
    ...
    >>> g = gt.random_graph(10000, lambda: sample_k(40), model="probabilistic",
    ...                     vertex_corr=lambda i, j: (sin(i / pi) * sin(j / pi) + 1) / 2,
    ...                     directed=False, n_iter=100)
    >>> h = gt.avg_neighbour_corr(g, "out", "out")
    >>> clf()
    >>> xlabel("Source out-degree")
    <...>
    >>> ylabel("Target out-degree")
    <...>
    >>> errorbar(h[2][:-1], h[0], yerr=h[1], fmt="o")
    <...>
    >>> savefig("avg_corr.pdf")

    .. testcode::
       :hide:

       savefig("avg_corr.png")

    .. figure:: avg_corr.*
        :align: center

        Average out/out degree correlation.
    """

    ret = libgraph_tool_correlations.\
          vertex_avg_correlation(g._Graph__graph, _degree(g, deg_source),
                                 _degree(g, deg_target), _prop("e", g, weight),
                                 [float(x) for x in bins])
    return [ret[0], ret[1], ret[2][0]]


def avg_combined_corr(g, deg1, deg2, bins=[0, 1]):
    r"""
    Obtain the single-vertex combined correlation histogram for the given graph.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    deg1 : string or :class:`~graph_tool.PropertyMap`
        first degree type ("in", "out" or "total") or vertex property map.
    deg2 : string or :class:`~graph_tool.PropertyMap`
        second degree type ("in", "out" or "total") or vertex property map.
    bins : list (optional, default: [0, 1])
        Bins to be used for the first degrees. If the list has size 2, it is
        used as the constant width of an automatically generated bin range,
        starting from the first value.

    Returns
    -------
    bin_avg : :class:`~numpy.ndarray`
        Array with the deg2 average for the deg1 bins.
    bin_dev : :class:`~numpy.ndarray`
        Array with the standard deviation of the deg2 average for the deg1 bins.
    bins : :class:`~numpy.ndarray`
        The deg1 bins.

    Notes
    -----
    If enabled during compilation, this algorithm runs in parallel.


    See Also
    --------
    assortativity: assortativity coefficient
    scalar_assortativity: scalar assortativity coefficient
    corr_hist: vertex-vertex correlation histogram
    combined_corr_hist: combined single-vertex correlation histogram
    avg_neighbour_corr: average nearest-neighbour correlation
    avg_combined_corr: average combined single-vertex correlation

    Examples
    --------
    .. testcode::
       :hide:

       np.random.seed(42)
       gt.seed_rng(42)
       from pylab import *

    >>> def sample_k(max):
    ...     accept = False
    ...     while not accept:
    ...         i = randint(1,max+1)
    ...         j = randint(1,max+1)
    ...         accept = random() < (sin(i/pi)*sin(j/pi)+1)/2
    ...     return i,j
    ...
    >>> g = gt.random_graph(10000, lambda: sample_k(40))
    >>> h = gt.avg_combined_corr(g, "in", "out")
    >>> clf()
    >>> xlabel("In-degree")
    <...>
    >>> ylabel("Out-degree")
    <...>
    >>> errorbar(h[2][:-1], h[0], yerr=h[1], fmt="o")
    <...>
    >>> savefig("combined_avg_corr.pdf")

    .. testcode::
       :hide:

       savefig("combined_avg_corr.png")

    .. figure:: combined_avg_corr.*
        :align: center

        Average combined in/out-degree correlation.
    """

    ret = libgraph_tool_correlations.\
          vertex_avg_combined_correlation(g._Graph__graph, _degree(g, deg1),
                                          _degree(g, deg2),
                                          [float(x) for x in bins])
    return [ret[0], ret[1], ret[2][0]]
