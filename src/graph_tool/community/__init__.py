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
``graph_tool.community`` - Community structure
----------------------------------------------

This module contains algorithms for the computation of community structure on
graphs.


Stochastic blockmodel inference
+++++++++++++++++++++++++++++++

Non-hierarchical models
=======================

Summary
^^^^^^^

.. autosummary::
   :nosignatures:

   minimize_blockmodel_dl
   BlockState
   OverlapBlockState
   mcmc_sweep
   MinimizeState
   multilevel_minimize
   collect_vertex_marginals
   collect_edge_marginals
   mf_entropy
   bethe_entropy
   model_entropy
   get_max_B
   get_akc
   condensation_graph
   get_block_edge_gradient

Hierarchical models
===================

Summary
^^^^^^^

.. autosummary::
   :nosignatures:

   minimize_nested_blockmodel_dl
   NestedBlockState
   NestedMinimizeState
   init_nested_state
   nested_mcmc_sweep
   nested_tree_sweep
   get_hierarchy_tree


Modularity-based community detection
++++++++++++++++++++++++++++++++++++

Summary
=======

.. autosummary::
   :nosignatures:

   community_structure
   modularity


Contents
++++++++
"""

from __future__ import division, absolute_import, print_function
import sys
if sys.version_info < (3,):
    range = xrange

from .. dl_import import dl_import
#dl_import("from . import libgraph_tool_community")

from .. import _degree, _prop, Graph, GraphView, libcore, _get_rng
import random
import sys

__all__ = ["minimize_blockmodel_dl",
           "BlockState",
           "mcmc_sweep",
           "MinimizeState",
           "multilevel_minimize",
           "collect_edge_marginals",
           "collect_vertex_marginals",
           "mf_entropy",
           "bethe_entropy",
           "model_entropy",
           "get_max_B",
           "get_akc",
           "condensation_graph",
           "OverlapBlockState",
           "minimize_nested_blockmodel_dl",
           "NestedBlockState",
           "NestedMinimizeState",
           "init_nested_state",
           "nested_mcmc_sweep",
           "nested_tree_sweep",
           "get_hierarchy_tree",
           "get_block_edge_gradient",
           "community_structure",
           "modularity"]

from . blockmodel import minimize_blockmodel_dl, BlockState, mcmc_sweep, \
    multilevel_minimize, model_entropy, get_max_B, get_akc, condensation_graph, \
    collect_edge_marginals, collect_vertex_marginals, bethe_entropy, mf_entropy, MinimizeState

from . overlap_blockmodel import OverlapBlockState, get_block_edge_gradient

from . nested_blockmodel import NestedBlockState, NestedMinimizeState, init_nested_state, \
    nested_mcmc_sweep, nested_tree_sweep, minimize_nested_blockmodel_dl, get_hierarchy_tree

def community_structure(g, n_iter, n_spins, gamma=1.0, corr="erdos",
                        spins=None, weight=None, t_range=(100.0, 0.01),
                        verbose=False, history_file=None):
    r"""Obtain the community structure for the given graph, using a Potts model
    approach, which is a generalization of modularity maximization).

    .. warning::

       **The use of this function is discouraged.** Although community detection
       based on modularity maximization is very common, it is very
       problematic. It will find high-scoring partitions where there is none
       [guimera-modularity-2004]_, and at the same time will not find actual
       structure in large graphs [fortunato-resolution-2007]_. Furthermore, in
       many empirical networks, the partitions found in this way are largely
       meaningless [good-performance-2010]_.

       One should use instead methods based on statistical inference
       (i.e. :func:`~graph_tool.community.minimize_blockmodel_dl` and
       :func:`~graph_tool.community.minimize_nested_blockmodel_dl`).

    Parameters
    ----------
    g :  :class:`~graph_tool.Graph`
        Graph to be used.
    n_iter : int
        Number of iterations.
    n_spins : int
        Number of maximum spins to be used.
    gamma : float (optional, default: 1.0)
        The :math:`\gamma` parameter of the hamiltonian.
    corr : string (optional, default: "erdos")
        Type of correlation to be assumed: Either "erdos", "uncorrelated" and
        "correlated".
    spins : :class:`~graph_tool.PropertyMap`
        Vertex property maps to store the spin variables. If this is specified,
        the values will not be initialized to a random value.
    weight : :class:`~graph_tool.PropertyMap` (optional, default: None)
        Edge property map with the optional edge weights.
    t_range : tuple of floats (optional, default: (100.0, 0.01))
        Temperature range.
    verbose : bool (optional, default: False)
        Display verbose information.
    history_file : string (optional, default: None)
        History file to keep information about the simulated annealing.

    Returns
    -------
    spins : :class:`~graph_tool.PropertyMap`
        Vertex property map with the spin values.

    See Also
    --------
    community_structure: Obtain the community structure
    modularity: Calculate the network modularity
    condensation_graph: Network of communities, or blocks

    Notes
    -----
    The method of community detection covered here is an implementation of what
    was proposed in [reichard-statistical-2006]_. It
    consists of a `simulated annealing`_ algorithm which tries to minimize the
    following hamiltonian:

    .. math::

        \mathcal{H}(\{\sigma\}) = - \sum_{i \neq j} \left(A_{ij} -
        \gamma p_{ij}\right) \delta(\sigma_i,\sigma_j)

    where :math:`p_{ij}` is the probability of vertices i and j being connected,
    which reduces the problem of community detection to finding the ground
    states of a Potts spin-glass model. It can be shown that minimizing this
    hamiltonan, with :math:`\gamma=1`, is equivalent to maximizing
    Newman's modularity ([newman-modularity-2006]_). By increasing the parameter
    :math:`\gamma`, it's possible also to find sub-communities.

    It is possible to select three policies for choosing :math:`p_{ij}` and thus
    choosing the null model: "erdos" selects a Erdos-Reyni random graph,
    "uncorrelated" selects an arbitrary random graph with no vertex-vertex
    correlations, and "correlated" selects a random graph with average
    correlation taken from the graph itself. Optionally a weight property
    can be given by the `weight` option.


    The most important parameters for the algorithm are the initial and final
    temperatures (`t_range`), and total number of iterations (`max_iter`). It
    normally takes some trial and error to determine the best values for a
    specific graph. To help with this, the `history` option can be used, which
    saves to a chosen file the temperature and number of spins per iteration,
    which can be used to determined whether or not the algorithm converged to
    the optimal solution. Also, the `verbose` option prints the computation
    status on the terminal.

    .. note::

        If the spin property already exists before the computation starts, it's
        not re-sampled at the beginning. This means that it's possible to
        continue a previous run, if you saved the graph, by properly setting
        `t_range` value, and using the same `spin` property.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------

    This example uses the network :download:`community.xml <community.xml>`.

    >>> from pylab import *
    >>> from numpy.random import seed
    >>> seed(42)
    >>> g = gt.load_graph("community.xml")
    >>> pos = g.vertex_properties["pos"]
    >>> spins = gt.community_structure(g, 10000, 20, t_range=(5, 0.1))
    >>> gt.graph_draw(g, pos=pos, vertex_fill_color=spins, output_size=(420, 420), output="comm1.pdf")
    <...>

    .. testcode::
       :hide:

       gt.graph_draw(g, pos=pos, vertex_fill_color=spins, output_size=(420, 420), output="comm1.png")

    >>> spins = gt.community_structure(g, 10000, 40, t_range=(5, 0.1), gamma=2.5)
    >>> gt.graph_draw(g, pos=pos, vertex_fill_color=spins, output_size=(420, 420), output="comm2.pdf")
    <...>

    .. testcode::
       :hide:

       gt.graph_draw(g, pos=pos, vertex_fill_color=spins, output_size=(420, 420), output="comm2.png")

    .. figure:: comm1.*
       :align: center

       The community structure with :math:`\gamma=1`.


    .. figure:: comm2.*
       :align: center

       The community structure with :math:`\gamma=2.5`.

    References
    ----------

    .. [guimera-modularity-2004] Roger Guimerà, Marta Sales-Pardo, and Luís
       A. Nunes Amaral, "Modularity from fluctuations in random graphs and
       complex networks", Phys. Rev. E 70, 025101(R) (2004),
       :doi:`10.1103/PhysRevE.70.025101`
    .. [fortunato-resolution-2007] Santo Fortunato and Marc Barthélemy,
       "Resolution limit in community detection", Proc. Natl. Acad. Sci. USA
       104(1): 36–41 (2007), :doi:`10.1073/pnas.0605965104`
    .. [good-performance-2010] Benjamin H. Good, Yves-Alexandre de Montjoye, and
       Aaron Clauset, "Performance of modularity maximization in practical
       contexts", Phys. Rev. E 81, 046106 (2010),
       :doi:`10.1103/PhysRevE.81.046106`
    .. [reichard-statistical-2006] Joerg Reichardt and Stefan Bornholdt,
       "Statistical Mechanics of Community Detection", Phys. Rev. E 74
       016110 (2006), :doi:`10.1103/PhysRevE.74.016110`, :arxiv:`cond-mat/0603718`
    .. [newman-modularity-2006] M. E. J. Newman, "Modularity and community
       structure in networks", Proc. Natl. Acad. Sci. USA 103, 8577-8582 (2006),
       :doi:`10.1073/pnas.0601602103`, :arxiv:`physics/0602124`
    .. _simulated annealing: http://en.wikipedia.org/wiki/Simulated_annealing
    """

    if spins is None:
        spins = g.new_vertex_property("int32_t")
    if history_file is None:
        history_file = ""
    ug = GraphView(g, directed=False)
    libgraph_tool_community.community_structure(ug._Graph__graph, gamma, corr,
                                                n_iter, t_range[1], t_range[0],
                                                n_spins, _get_rng(),
                                                verbose, history_file,
                                                _prop("e", ug, weight),
                                                _prop("v", ug, spins))
    return spins


def modularity(g, prop, weight=None):
    r"""
    Calculate Newman's modularity.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    prop : :class:`~graph_tool.PropertyMap`
        Vertex property map with the community partition.
    weight : :class:`~graph_tool.PropertyMap` (optional, default: None)
        Edge property map with the optional edge weights.

    Returns
    -------
    modularity : float
        Newman's modularity.

    See Also
    --------
    community_structure: obtain the community structure
    modularity: calculate the network modularity
    condensation_graph: Network of communities, or blocks

    Notes
    -----

    Given a specific graph partition specified by `prop`, Newman's modularity
    [newman-modularity-2006]_ is defined by:

    .. math::

          Q = \frac{1}{2E} \sum_r e_{rr}- \frac{e_r^2}{2E}

    where :math:`e_{rs}` is the number of edges which fall between
    vertices in communities s and r, or twice that number if :math:`r = s`, and
    :math:`e_r = \sum_s e_{rs}`.

    If weights are provided, the matrix :math:`e_{rs}` corresponds to the sum
    of edge weights instead of number of edges, and the value of :math:`E`
    becomes the total sum of edge weights.

    Examples
    --------
    >>> from pylab import *
    >>> from numpy.random import seed
    >>> seed(42)
    >>> g = gt.load_graph("community.xml")
    >>> b = gt.community_structure(g, 10000, 10)
    >>> gt.modularity(g, b)
    0.5353141885624041

    References
    ----------
    .. [newman-modularity-2006] M. E. J. Newman, "Modularity and community
       structure in networks", Proc. Natl. Acad. Sci. USA 103, 8577-8582 (2006),
       :doi:`10.1073/pnas.0601602103`, :arxiv:`physics/0602124`
    """

    ug = GraphView(g, directed=False)
    m = libgraph_tool_community.modularity(ug._Graph__graph,
                                           _prop("e", ug, weight),
                                           _prop("v", ug, prop))
    return m
