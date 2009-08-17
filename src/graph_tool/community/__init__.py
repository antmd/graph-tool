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
``graph_tool.community`` - Community structure
----------------------------------------------

This module contains algorithms for the computation of community structure on
graphs.
"""

from .. dl_import import dl_import
dl_import("import libgraph_tool_community")

from .. core import _degree, _prop, Graph, libcore
import random, sys

__all__ = ["community_structure", "modularity", "condensation_graph"]


def community_structure(g, n_iter, n_spins, gamma=1.0, corr= "erdos",
                        spins=None, weight=None, t_range=(100.0, 0.01),
                        verbose=False, history_file=None, seed=0):
    r"""
    Obtain the community structure for the given graph, used a Potts model
    approach.

    Parameters
    ----------
    g : Graph
        Graph to be used.
    n_iter : int
        Number of iterations.
    n_spins : int
        Number of maximum spins to be used.
    gamma : float (optional, default: 1.0)
        The :math:`\gamma` parameter of the hamiltonian.
    corr : string (optional, default: "erdos")
        Type of correlation to be assumed: Either "erdos", "random" and
        "correlated".
    spins : PropertyMap
        Vertex property maps to store the spin variables. If this is specified,
        the values will not be initialized to a random value.
    weight : PropertyMap (optional, default: None)
        Edge property map with the optional edge weights.
    t_range : tuple of floats (optional, default: (100.0, 0.01))
        Temperature range.
    verbose : bool (optional, default: False)
        Display verbose information.
    history_file : string (optional, default: None)
        History file to keep information about the simulated annealing.

    Returns
    -------
    spins : PropertyMap
        Vertex property map with the spin values.

    See Also
    --------
    community_structure: obtain the community structure
    modularity: calculate the network modularity
    condensation_graph: network of communities

    Notes
    -----
    The method of community detection covered here is an implementation of what
    was proposed in [reichard_statistical_2006]_. It
    consists of a `simulated annealing`_ algorithm which tries to minimize the
    following hamiltonian:

    .. math::

        \mathcal{H}(\{\sigma\}) = - \sum_{i \neq j} \left(A_{ij} -
        \gamma p_{ij}\right) \delta(\sigma_i,\sigma_j)

    where :math:`p_{ij}` is the probability of vertices i and j being connected,
    which reduces the problem of community detection to finding the ground
    states of a Potts spin-glass model. It can be shown that minimizing this
    hamiltonan, with :math:`\gamma=1`, is equivalent to maximizing
    Newman's modularity ([newman_modularity_2006]_). By increasing the parameter
    :math:`\gamma`, it's possible also to find sub-communities.

    It is possible to select three policies for choosing :math:`p_{ij}` and thus
    choosing the null model: "random" selects a Erdos-Reyni random graph,
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
    >>> from pylab import *
    >>> from numpy.random import seed
    >>> seed(42)
    >>> g = gt.load_graph("community.xml")
    >>> pos = (g.vertex_properties["pos_x"], g.vertex_properties["pos_y"])
    >>> spins = gt.community_structure(g, 10000, 20, t_range=(5, 0.1),
    ...                                history_file="community-history1")
    >>> gt.graph_draw(g, pos=pos, pin=True, vsize=0.3, vcolor=spins,
    ...               output="comm1.png")
    (...)
    >>> spins = gt.community_structure(g, 10000, 40, t_range=(5, 0.1),
    ...                                gamma=2.5,
    ...                                history_file="community-history2")
    >>> gt.graph_draw(g, pos=pos, pin=True, vsize=0.3, vcolor=spins,
    ...               output="comm2.png")
    (...)
    >>> clf()
    >>> xlabel("iterations")
    <...>
    >>> ylabel("number of communities")
    <...>
    >>> a = loadtxt("community-history1").transpose()
    >>> plot(a[0], a[2])
    [...]
    >>> savefig("comm1-hist.png")
    >>> clf()
    >>> xlabel("iterations")
    <...>
    >>> ylabel("number of communities")
    <...>
    >>> a = loadtxt("community-history2").transpose()
    >>> plot(a[0], a[2])
    [...]
    >>> savefig("comm2-hist.png")

    .. figure:: comm1.png
        :align: center

        Community structure with :math:`\gamma=1`.

    .. figure:: comm1-hist.png
        :align: center

        Algorithm evolution with :math:`\gamma=1`

    .. figure:: comm2.png
        :align: center

        Community structure with :math:`\gamma=2.5`.

    .. figure:: comm2-hist.png
        :align: center

        Algorithm evolution with :math:`\gamma=2.5`

    References
    ----------
    .. [reichard_statistical_2006] Joerg Reichardt and Stefan Bornholdt,
       "Statistical Mechanics of Community Detection", Phys. Rev. E 74 (2006)
       016110, arXiv:cond-mat/0603718
    .. [newman_modularity_2006] M. E. J. Newman, "Modularity and community
       structure in networks", Proc. Natl. Acad. Sci. USA 103, 8577-8582 (2006),
       arXiv:physics/0602124
    .. _simulated annealing: http://en.wikipedia.org/wiki/Simulated_annealing
    """

    if spins == None:
        spins = g.new_vertex_property("int32_t")
        new_spins = True
    else:
        new_spins = False
    if history_file == None:
        history_file = ""
    if seed != 0:
        seed = random.randint(0, sys.maxint)
    libgraph_tool_community.community_structure(g._Graph__graph, gamma, corr,
                                                n_iter, t_range[1], t_range[0],
                                                n_spins, new_spins, seed,
                                                verbose, history_file,
                                                _prop("e", g, weight),
                                                _prop("v", g, spins))
    return spins

def modularity(g, prop, weight=None):
    r"""
    Calculate Newman's modularity.

    Parameters
    ----------
    g : Graph
        Graph to be used.
    prop : PropertyMap
        Vertex property map with the community partition.
    weight : PropertyMap (optional, default: None)
        Edge property map with the optional edge weights.

    Returns
    -------
    modularity : float
        Newman's modularity.

    See Also
    --------
    community_structure: obtain the community structure
    modularity: calculate the network modularity
    condensation_graph: network of communities

    Notes
    -----

    Given a specific graph partition specified by `prop`, Newman's modularity
    ([newman_modularity_2006]_) is defined by:

    .. math::

          Q = \sum_s e_{ss}-\left(\sum_r e_{rs}\right)^2

    where :math:`e_{rs}` is the fraction of edges which fall between
    vertices with spin s and r.

    If enabled during compilation, this algorithm runs in parallel.

    Examples
    --------
    >>> from pylab import *
    >>> from numpy.random import seed
    >>> seed(42)
    >>> g = gt.load_graph("community.dot")
    >>> spins = gt.community_structure(g, 10000, 10)
    >>> gt.modularity(g, spins)
    0.53531418856240398

    References
    ----------
    .. [newman_modularity_2006] M. E. J. Newman, "Modularity and community
       structure in networks", Proc. Natl. Acad. Sci. USA 103, 8577-8582 (2006),
       arXiv:physics/0602124
    """

    m = libgraph_tool_community.modularity(g._Graph__graph,
                                           _prop("e", g, weight),
                                           _prop("v", g, prop))
    return m

def condensation_graph(g, prop, weight=None):
    r"""
    Obtain the condensation graph, where each vertex with the same 'prop' value
    is condensed in one vertex.

    Parameters
    ----------
    g : Graph
        Graph to be used.
    prop : PropertyMap
        Vertex property map with the community partition.
    weight : PropertyMap (optional, default: None)
        Edge property map with the optional edge weights.

    Returns
    -------
    condensation_graph : Graph
        The community network
    vcount : PropertyMap
        A vertex property map with the vertex count for each community.
    ecount : PropertyMap
        An edge property map with the inter-community edge count for each edge.

    See Also
    --------
    community_structure: obtain the community structure
    modularity: calculate the network modularity
    condensation_graph:  network of communities

    Notes
    -----
    Each vertex in the condensation graph represents one community in the
    original graph (vertices with the same 'prop' value'), and the edges
    represent existent edges between vertices of the respective communities in
    the original graph.

    Examples
    --------
    >>> from pylab import *
    >>> from numpy.random import poisson, seed
    >>> seed(42)
    >>> g = gt.random_graph(1000, lambda: poisson(3), directed=False)
    >>> spins = gt.community_structure(g, 10000, 100)
    >>> ng = gt.condensation_graph(g, spins)
    >>> size = ng[0].new_vertex_property("double")
    >>> size.get_array()[:] = log(ng[1].get_array()+1)
    >>> gt.graph_draw(ng[0], vsize=size, vcolor=size, splines=True,
    ...               eprops={"len":20, "penwidth":10}, vprops={"penwidth":10},
    ...               output="comm-network.png")
    (...)

    .. figure:: comm-network.png
        :align: center

        Community network of a random graph. The color and sizes of the nodes
        indicate the size of the corresponding community.
    """
    gp = Graph()
    vcount = gp.new_vertex_property("int32_t")
    if weight != None:
        ecount = gp.new_edge_property("double")
    else:
        ecount = gp.new_edge_property("int32_t")
    libgraph_tool_community.community_network(g._Graph__graph,
                                              gp._Graph__graph,
                                              _prop("v", g, prop),
                                              _prop("v", g, vcount),
                                              _prop("e", g, ecount),
                                              _prop("e", g, weight))
    return gp, vcount, ecount

