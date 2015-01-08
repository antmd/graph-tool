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

from __future__ import division, absolute_import, print_function
import sys
if sys.version_info < (3,):
    range = xrange

from .. import _degree, _prop, Graph, GraphView, libcore, _get_rng, PropertyMap
import random
from numpy import *
import numpy
from collections import defaultdict

from .. import group_vector_property, ungroup_vector_property

from .. dl_import import dl_import
dl_import("from . import libgraph_tool_community as libcommunity")

from . blockmodel import *
from . blockmodel import __test__

class OverlapBlockState(BlockState):
    r"""This class encapsulates the overlapping block state of a given graph.

    This must be instantiated and used by functions such as :func:`mcmc_sweep`.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be modelled.
    b : :class:`~graph_tool.PropertyMap` or :class:`numpy.ndarray` (optional, default: ``None``)
        Initial block labels on the vertices or half-edges. If not supplied, it
        will be randomly sampled.
        If the value passed is a vertex property map, it will be assumed to be a
        non-overlapping partition of the vertices. If it is an edge property
        map, it should contain a vector for each edge, with the block labels at
        each end point (sorted according to their vertex index, in the case of
        undirected graphs, otherwise from source to target). If the value is an
        :class:`numpy.ndarray`, it will be assumed to correspond directly to a
        partition of the list of half-edges.
    B : ``int`` (optional, default: ``None``)
        Number of blocks. If not supplied it will be either obtained from the
        parameter ``b``, or set to the maximum possible value according to the
        minimum description length.
    clabel : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Constraint labels on the vertices. If supplied, vertices with different
        label values will not be clustered in the same group.
    deg_corr : ``bool`` (optional, default: ``True``)
        If ``True``, the degree-corrected version of the blockmodel ensemble will
        be assumed, otherwise the traditional variant will be used.
    max_BE : ``int`` (optional, default: ``1000``)
        If the number of blocks exceeds this number, a sparse representation of
        the block graph is used, which is slightly less efficient, but uses less
        memory,
    """

    def __init__(self, g, b=None, B=None, clabel=None, deg_corr=True,
                 max_BE=1000, **kwargs):

        BlockState._state_ref_count += 1

        # determine if there is a base graph, and overlapping structure
        self.base_g = kwargs.get("base_g", None)

        # overlapping information
        node_index = kwargs.get("node_index", None)
        node_in_degs = kwargs.get("node_in_degs", None)
        node_out_degs = kwargs.get("node_out_degs", None)
        half_edges = kwargs.get("half_edges", None)

        if node_index is not None and self.base_g is None:
            raise ValueError("Must specify base graph if node_index is specified...")

        # create overlapping structure
        if node_index is None:
            # keep base graph
            self.base_g = g

            # substitute provided graph by its half-edge graph
            g, b, node_index, half_edges = half_edge_graph(g, b, B)

        # create half edges set if absent
        if node_index is not None and half_edges is None:
            half_edges = self.base_g.new_vertex_property("vector<int>")
            libcommunity.get_nodeset_overlap(g._Graph__graph,
                                             _prop("v", g, node_index),
                                             _prop("v", self.base_g, half_edges))

        self.overlap = True
        self.node_index = node_index
        self.half_edges = half_edges
        if self.node_index is None:
            self.node_index = g.new_vertex_property("int")
        if self.half_edges is None:
            self.half_edges = g.new_vertex_property("vector<int>")

        # configure the main graph and block model parameters
        self.g = g

        self.E = self.g.num_edges()
        self.N = self.base_g.num_vertices()

        self.deg_corr = deg_corr

        # initialize weights to unity
        eweight = g.new_edge_property("int")
        eweight.fa = 1
        vweight = g.new_vertex_property("int")
        vweight.fa = 1
        self.eweight = eweight
        self.vweight = vweight
        self.is_weighted = False

        # ensure we have at most as many blocks as nodes
        if B is not None:
            B = min(B, self.g.num_vertices())

        if b is None:
            # create a random partition into B blocks.
            if B is None:
                B = get_max_B(self.N, self.E, directed=g.is_directed())
            B = min(B, self.g.num_vertices())
            ba = random.randint(0, B, self.g.num_vertices())
            ba[:B] = arange(B)        # avoid empty blocks
            if B < self.g.num_vertices():
                random.shuffle(ba)
            b = g.new_vertex_property("int")
            b.fa = ba
            self.b = b
        else:
            # if a partition is available, we will incorporate it.
            # in the overlapping case
            # at this point, *b* must correspond to the partition of
            # *half-edges*
            if isinstance(b, numpy.ndarray):
                self.b = g.new_vertex_property("int")
                self.b.fa = b
            else:
                b = b.copy(value_type="int")
                b = g.own_property(b)
                self.b = b
            if B is None:
                B = int(self.b.fa.max()) + 1

        # if B > self.N:
        #     raise ValueError("B > N!")

        if self.b.fa.max() >= B:
            raise ValueError("Maximum value of b is larger or equal to B!")

        # Construct block-graph
        self.bg = get_block_graph(g, B, self.b, self.vweight, self.eweight)
        self.bg.set_fast_edge_removal()

        self.mrs = self.bg.ep["count"]
        self.wr = self.bg.vp["count"]

        del self.bg.ep["count"]
        del self.bg.vp["count"]

        self.mrp = self.bg.degree_property_map("out", weight=self.mrs)

        if g.is_directed():
            self.mrm = self.bg.degree_property_map("in", weight=self.mrs)
        else:
            self.mrm = self.mrp

        self.vertices = libcommunity.get_vector(B)
        self.vertices.a = arange(B)
        self.B = B

        # if overlapping, the node counts must correspond to actual nodes,
        # not half-edges.
        self.wr.a = 0
        bv = self.get_overlap_blocks()[0]
        libcommunity.get_wr_overlap(self.base_g._Graph__graph,
                                    _prop("v", self.base_g, bv),
                                    _prop("v", self.bg, self.wr))

        self.overlap_stats = libcommunity.init_overlap_stats(self.g._Graph__graph,
                                                             _prop("v", self.g, self.b),
                                                             _prop("v", self.g, self.half_edges),
                                                             _prop("v", self.g, self.node_index),
                                                             self.B)

        if clabel is not None:
            if isinstance(clabel, PropertyMap):
                # if clabel is a property map, we will assume it constraints the *nodes*
                if __test__:
                    assert len(clabel.a) < self.g.num_vertices()
                self.clabel = self.node_index.copy()
                pmap(self.clabel, clabel)
            else:
                # if clabel is an array, we will assume it constraints the *half-edges*
                self.clabel = self.g.new_vertex_property("int")
                self.clabel.a = clabel
        else:
            self.clabel = self.g.new_vertex_property("int")

        self.emat = None
        if max_BE is None:
            max_BE = 1000
        self.max_BE = max_BE

        # used by mcmc_sweep()
        self.egroups = None
        self.nsampler = None
        self.sweep_vertices = None
        self.partition_stats = libcommunity.overlap_partition_stats()

        # computation cache
        libcommunity.init_safelog(int(5 * max(self.E, self.N)))
        libcommunity.init_xlogx(int(5 * max(self.E, self.N)))
        libcommunity.init_lgamma(int(3 * max(self.E, self.N)))

    def __del__(self):
        try:
            BlockState.__del__(self)
        except (TypeError, AttributeError):
            pass

    def __repr__(self):
        return "<OverlapBlockState object with %d blocks,%s for graph %s, at 0x%x>" % \
            (self.B, " degree corrected," if self.deg_corr else "",
             str(self.base_g), id(self))

    def __init_partition_stats(self, empty=True):
        if not empty:
            self.partition_stats = libcommunity.init_overlap_partition_stats(self.g._Graph__graph,
                                                                             _prop("v", self.g, self.b),
                                                                             _prop("e", self.g, self.eweight),
                                                                             self.N, self.B,
                                                                             self.overlap_stats)
        else:
            self.partition_stats = libcommunity.overlap_partition_stats()


    def copy(self, b=None, B=None, deg_corr=None, clabel=None, overlap=True):
        r"""Copies the block state. The parameters override the state properties, and
         have the same meaning as in the constructor. If ``overlap=False`` an
         instance of :class:`~graph_tool.community.BlockState` is returned."""

        if overlap:
            state = OverlapBlockState(self.g,
                                      eweight=self.eweight,
                                      vweight=self.vweight,
                                      b=self.b if b is None else b,
                                      B=(self.B if b is None else None) if B is None else B,
                                      clabel=self.clabel.a if clabel is None else clabel,
                                      deg_corr=self.deg_corr if deg_corr is None else deg_corr,
                                      base_g=self.base_g,
                                      half_edges=self.half_edges,
                                      node_index=self.node_index,
                                      max_BE=self.max_BE)
        else:
            state = BlockState(self.base_g,
                               b=b if b is not None else self.get_nonoverlap_blocks(),
                               B=B,
                               clabel=clabel if clabel is not None else None,
                               deg_corr=self.deg_corr if deg_corr is None else deg_corr,
                               max_BE=self.max_BE)

        if not state._BlockState__check_clabel():
            b = state.b.a + state.clabel.a * state.B
            continuous_map(b)
            state = state.copy(b=b)

            if __test__:
                assert state._BlockState__check_clabel()

        return state

    def __getstate__(self):
        state = dict(g=self.g,
                     eweight=self.eweight,
                     vweight=self.vweight,
                     b=self.b,
                     B=self.B,
                     clabel=array(self.clabel.a),
                     deg_corr=self.deg_corr,
                     base_g=self.base_g,
                     half_edges=self.half_edges,
                     node_index=self.node_index,
                     max_BE=self.max_BE)
        return state

    def __setstate__(self, state):
        self.__init__(**state)
        return state

    def get_block_state(self, b=None, vweight=False, overlap=False,
                        deg_corr=False):
        r"""Returns a :class:`~graph_tool.community.BlockState`` (or
        :class:`~graph_tool.community.OverlapBlockState`` if ``overlap==True``)
        corresponding to the block graph. The parameters have the same meaning
        as the in the constructor."""

        if not overlap:
            bg = self.bg.copy()
            mrs = bg.own_property(self.mrs.copy())
            wr = bg.own_property(self.wr.copy())
            state = BlockState(bg, eweight=mrs,
                               vweight=wr if vweight else None,
                               b=bg.vertex_index.copy("int") if b is None else b,
                               deg_corr=deg_corr,
                               clabel=self.get_bclabel(),
                               max_BE=self.max_BE)
            n_map = self.b.copy()
        else:
            ## FIXME: Move this to C++
            bg = Graph(directed=self.g.is_directed())
            bg.add_vertex(self.B)
            for e in self.g.edges():
                r, s = tuple(e)
                r = self.b[r]
                s = self.b[s]
                bg.add_edge(bg.vertex(r), bg.vertex(s))

            state = OverlapBlockState(self.g,
                                      b=self.b.a if b is None else b.copy(),
                                      deg_corr=deg_corr,
                                      clabel=self.clabel,
                                      base_g=bg,
                                      node_index=self.b,
                                      max_BE=self.max_BE)
            n_map = self.g.vertex_index.copy("int")
        return state, n_map

    def get_edge_blocks(self):
        r"""Returns an edge property map which contains the block labels pairs for each
        edge."""
        be = self.base_g.new_edge_property("vector<int>")
        libcommunity.get_be_overlap(self.base_g._Graph__graph,
                                    self.g._Graph__graph,
                                    _prop("e", self.base_g, be),
                                    _prop("v", self.g, self.b),
                                    _prop("v", self.g, self.node_index))
        return be

    def get_overlap_blocks(self):
        r"""Returns the mixed membership of each vertex.

        Returns
        -------
        bv : :class:`PropertyMap`
           A vector-valued vertex property map containing the block memberships
           of each node.
        bc_in : :class:`PropertyMap`
           The labelled in-degrees of each node, i.e. how many in-edges belong
           to each group, in the same order as the ``bv`` property above.
        bc_out : :class:`PropertyMap`
           The labelled out-degrees of each node, i.e. how many out-edges belong
           to each group, in the same order as the ``bv`` property above.
        bc_total : :class:`PropertyMap`
           The labelled total degrees of each node, i.e. how many incident edges
           belong to each group, in the same order as the ``bv`` property above.

        """
        bv = self.base_g.new_vertex_property("vector<int>")
        bc_in = self.base_g.new_vertex_property("vector<int>")
        bc_out = self.base_g.new_vertex_property("vector<int>")
        bc_total = self.base_g.new_vertex_property("vector<int>")
        libcommunity.get_bv_overlap(self.g._Graph__graph,
                                    _prop("v", self.g, self.b),
                                    _prop("v", self.g, self.node_index),
                                    _prop("v", self.base_g, bv),
                                    _prop("v", self.base_g, bc_in),
                                    _prop("v", self.base_g, bc_out),
                                    _prop("v", self.base_g, bc_total))
        return bv, bc_in, bc_out, bc_total

    def get_nonoverlap_blocks(self):
        r"""Returns a scalar-valued vertex property map with the block mixture
        represented as a single number."""

        bv = self.get_overlap_blocks()[0]
        b = self.base_g.new_vertex_property("int")
        libcommunity.get_overlap_split(self.base_g._Graph__graph,
                                       _prop("v", self.base_g, bv),
                                       _prop("v", self.base_g, b))
        return b

    # def get_overlap_projection(self):
    #     b = self.b.copy()
    #     libcommunity.get_overlap_proj(self.g._Graph__graph,
    #                                   self.bg._Graph__graph,
    #                                   _prop("e", self.bg, self.mrs),
    #                                   _prop("v", self.g, b))
    #     return b

    def entropy(self, complete=True, dl=False, partition_dl=False,
                multigraph=True, norm=True, dl_ent=False, **kwargs):
        r"""Calculate the entropy associated with the current overlapping partition.

        Parameters
        ----------
        complete : ``bool`` (optional, default: ``False``)
            If ``True``, the complete entropy will be returned, including constant
            terms not relevant to the block partition.
        dl : ``bool`` (optional, default: ``False``)
            If ``True``, the full description length will be returned.
        partition_dl : ``bool`` (optional, default: ``False``)
            If ``True``, and ``dl == True`` only the partition description
            length will be considered.
        multigraph : ``bool`` (optional, default: ``False``)
            If ``True``, the multigraph entropy will be used. Only has an effect
            if ``dense == True``.
        norm : ``bool`` (optional, default: ``True``)
            If ``True``, the entropy will be "normalized" by dividing by the
            number of edges.
        dl_ent : ``bool`` (optional, default: ``False``)
            If ``True``, the description length of the degree sequence will be
            approximated by its entropy.

        Notes
        -----
        For the traditional blockmodel (``deg_corr == False``), the entropy is
        given by

        .. math::

          \mathcal{S}_t &\cong E - \frac{1}{2} \sum_{rs}e_{rs}\ln\left(\frac{e_{rs}}{n_rn_s}\right), \\
          \mathcal{S}^d_t &\cong E - \sum_{rs}e_{rs}\ln\left(\frac{e_{rs}}{n_rn_s}\right),

        for undirected and directed graphs, respectively, where :math:`e_{rs}`
        is the number of edges from block :math:`r` to :math:`s` (or the number
        of half-edges for the undirected case when :math:`r=s`), and :math:`n_r
        = \sum_ib_i^r` is the number of vertices in block :math:`r`, with
        :math:`\vec{b}_i` being a binary vector with :math:`B` entries
        describing the mixed membership of node :math:`i`. Note that because of
        the possible overlaps, we have :math:`\sum_rn_r \ge N`.

        For the degree-corrected variant with "hard" degree constraints the
        equivalent expressions are

        .. math::

            \mathcal{S}_c &\cong -E -\sum_{ir}\ln k_i^r! - \frac{1}{2} \sum_{rs}e_{rs}\ln\left(\frac{e_{rs}}{e_re_s}\right), \\
            \mathcal{S}^d_c &\cong -E -\sum_{ir}\left(\ln {k^+}_i^r! + \ln {k^-}_i^r! \right) - \sum_{rs}e_{rs}\ln\left(\frac{e_{rs}}{e^+_re^-_s}\right),

        where :math:`e_r = \sum_se_{rs}` is the number of half-edges incident on
        block :math:`r`, :math:`e^+_r = \sum_se_{rs}` and :math:`e^-_r =
        \sum_se_{sr}` are the numbers of out- and in-edges adjacent to block
        :math:`r`, respectively, and :math:`k_i^r` is the degree of node
        :math:`i` of type :math:`r`, whereas :math:`{k^+}_i^r` and
        :math:`{k^-}_i^r` are the labeled out- and in-degrees, respectively.

        If ``multigraph == True``, the entropy used will be of the "Poisson"
        model, with the additional term:

        .. math::

            {\mathcal{S}_{cm}^{(d)}} = \mathcal{S}_c^{(d)} + \sum_{i>j}\sum_{rs} \ln A^{rs}_{ij}! + \sum_i\sum_{rs} \ln A^{rs}_{ii}!!


        If ``complete == False`` only the last term of the equations above will
        be returned.

        If ``dl == True``, the description length :math:`\mathcal{L}_t` of the
        model will be returned as well:

        .. math::

           \mathcal{L}_t = \ln\Omega_m + \ln\left(\!\!{D \choose N}\!\!\right) + \sum_d \ln {\left(\!\!{{B\choose d}\choose n_d}\!\!\right)} + \ln N! - \sum_{\vec{b}}\ln n_{\vec{b}}!,

        where :math:`d \equiv |\vec{b}|_1 = \sum_rb_r` is the mixture
        size, :math:`n_d` is the number of nodes in a mixture of size :math:`d`,
        :math:`D` is the maximum value of :math:`d`, :math:`n_{\vec{b}}` is the
        number of nodes in mixture :math:`\vec{b}`, and

        .. math::

            \Omega_m = \left(\!\!{\left(\!{B \choose 2}\!\right) \choose E}\!\!\right)

        is the number of different :math:`e_{rs}` matrices for undirected
        graphs, whereas for directed graphs we have

        .. math::

            \Omega_m = \left(\!\!{B^2 \choose E}\!\!\right).


        Note that for the degree-corrected version the description length is

        .. math::

            \mathcal{L}_c = \mathcal{L}_t + \sum_r\ln\left(\!\!{m_r \choose e_r}\!\!\right) + \sum_{\vec{b}}\min\left(\mathcal{L}^{(1)}_{\vec{b}}, \mathcal{L}^{(2)}_{\vec{b}}\right),

        where :math:`m_r` is the number of non-empty mixtures which contain type :math:`r`, and

        .. math::

            \mathcal{L}^{(1)}_{\vec{b}} &= \sum_r\ln{\left(\!\!{n_{\vec{b}}\choose e^r_{\vec{b}}}\!\!\right)}, \\
            \mathcal{L}^{(2)}_{\vec{b}} &= \sum_rb_r\ln\Xi_{\vec{b}}^r + \ln n_{\vec{b}}! - \sum_{\vec{k}} \ln n^{\vec{b}}_{\vec{k}}!,

        and :math:`\ln\Xi_{\vec{b}}^r \simeq 2\sqrt{\zeta(2)e_{\vec{b}}^r}`,
        where :math:`\zeta(x)` is the `Riemann zeta function
        <https://en.wikipedia.org/wiki/Riemann_zeta_function>`_, and
        :math:`n^{\vec{b}}_{\vec{k}}` is the number of nodes in mixture
        :math:`\vec{b}` with labelled degree :math:`\vec{k}`. For directed
        graphs we have instead :math:`\vec{k} \to (\vec{k}^-, \vec{k}^+)`, and
        :math:`\ln\Xi_{\vec{b}}^r \to \ln{\Xi^+}_{\vec{b}}^r +
        \ln{\Xi^-}_{\vec{b}}^r`.

        If ``dl_ent=True`` is passed, this will be approximated instead by

        .. math::

            \mathcal{L}_c \simeq \mathcal{L}_t - \sum_{\vec{b}}n_{\vec{b}}\sum_{\vec{k}}p^{\vec{b}}_{\vec{k}}\ln p^{\vec{b}}_{\vec{k}},

        where :math:`p^{\vec{b}}_{\vec{k}} = n^{\vec{b}}_{\vec{k}} / n_{\vec{b}}`.

        If ``complete == False`` constants in the above equations which do not
        depend on the partition of the nodes will be omitted.

        Note that in all cases if ``norm == True`` the value returned
        corresponds to the entropy `per edge`, i.e. :math:`(\mathcal{S}_{t/c}\;[\,+\,\mathcal{L}_{t/c}])/ E`.
        """


        xi_fast = kwargs.get("xi_fast", False)
        dl_deg_alt = kwargs.get("dl_deg_alt", True)  # compare the two deg encodings
        dense = kwargs.get("dense", False)

        if dense:
            raise NotImplementedError("Dense entropy for overlapping model not yet implemented")

            S = libcommunity.entropy_dense(self.bg._Graph__graph,
                                           _prop("e", self.bg, self.mrs),
                                           _prop("v", self.bg, self.wr),
                                           multigraph)
        else:
            S = libcommunity.entropy(self.bg._Graph__graph,
                                     _prop("e", self.bg, self.mrs),
                                     _prop("v", self.bg, self.mrp),
                                     _prop("v", self.bg, self.mrm),
                                     _prop("v", self.bg, self.wr),
                                     self.deg_corr)

            if complete:
                if self.deg_corr:
                    S -= self.E
                else:
                    S += self.E

            if multigraph:
                S += libcommunity.overlap_parallel_entropy(self.g._Graph__graph,
                                                           _prop("v", self.g, self.b),
                                                           self.overlap_stats)

            if self.deg_corr:
                S += libcommunity.deg_entropy_term(self.g._Graph__graph,
                                                   _prop("v", self.g, self.b),
                                                   self.overlap_stats,
                                                   self.N)

        if __test__:
            assert not isnan(S) and not isinf(S), "invalid entropy %g (%s) " % (S, str(dict(complete=complete,
                                                                                            random=random, dl=dl,
                                                                                            partition_dl=partition_dl,
                                                                                            dense=dense, multigraph=multigraph,
                                                                                            norm=norm)))
        if dl:
            N = self.base_g.num_vertices()
            E = self.E
            if not partition_dl:
                S += model_entropy(self.B, N, E, directed=self.g.is_directed(),
                                   nr=self.wr.a) * E - partition_entropy(self.B, N, self.wr.a)
                if __test__:
                    assert not isnan(S) and not isinf(S), "invalid entropy %g (%s) " % (S, str(dict(complete=complete,
                                                                                                    random=random, dl=dl,
                                                                                                    partition_dl=partition_dl,
                                                                                                    dense=dense, multigraph=multigraph,
                                                                                                    norm=norm)))

            if self.deg_corr:
                if self.partition_stats.is_enabled():
                    S += self.partition_stats.get_deg_dl(dl_ent, dl_deg_alt, xi_fast)
                else:
                    self.__init_partition_stats(empty=False)
                    S += self.partition_stats.get_deg_dl(dl_ent, dl_deg_alt, xi_fast)
                    self.__init_partition_stats(empty=True)

                if __test__:
                    assert not isnan(S) and not isinf(S), "invalid entropy %g (%s) " % (S, str(dict(complete=complete,
                                                                                                    random=random, dl=dl,
                                                                                                    partition_dl=partition_dl,
                                                                                                    dense=dense, multigraph=multigraph,
                                                                                                    norm=norm)))
            if self.partition_stats.is_enabled():
                S += self.partition_stats.get_partition_dl()
            else:
                self.__init_partition_stats(empty=False)
                S += self.partition_stats.get_partition_dl()
                self.__init_partition_stats(empty=True)

        if __test__:
            assert not isnan(S) and not isinf(S), "invalid entropy %g (%s) " % (S, str(dict(complete=complete,
                                                                                            random=random, dl=dl,
                                                                                            partition_dl=partition_dl,
                                                                                            dense=dense, multigraph=multigraph,
                                                                                            norm=norm)))
        if norm:
            return S / self.E
        else:
            return S

def half_edge_graph(g, b=None, B=None):
    r"""Generate a half-edge graph, where each half-edge is represented by a node,
    and an edge connects the half-edges like in the original graph."""

    E = g.num_edges()

    b_array = None
    if b is None:
        # if no partition is given, obtain a random one.
        ba = random.randint(0, B, 2 * E)
        ba[:B] = arange(B)        # avoid empty blocks
        if B < len(ba):
            random.shuffle(ba)
        b = ba

    if isinstance(b, numpy.ndarray):
        # if given an array, assume it corresponds to the *final* half-edge
        # partitions
        b_array = b
        b = g.new_vertex_property("int")

    if b.key_type() == "v":
        # If a vertex partition is given, we convert it into a
        # non-overlapping edge partition
        be = g.new_edge_property("vector<int>")
        libcommunity.get_be_from_b_overlap(g._Graph__graph,
                                           _prop("e", g, be),
                                           _prop("v", g, b))
        b = be
    else:
        # If an half-edge partition is provided, we incorporate it
        b = b.copy(value_type="vector<int32_t>")

    if B is None:
        if b_array is None:
            bs, bt = ungroup_vector_property(b, [0, 1])
            B = int(max(bs.fa.max(), bt.fa.max())) + 1
        else:
            B = b_array.max() + 1

    bs, bt = ungroup_vector_property(b, [0, 1])

    if bs.fa.max() >= B or bt.fa.max() >= B or (b_array is not None and b_array.max() >= B):
        raise ValueError("Maximum value of b is larger or equal to B!")

    eg = Graph(directed=g.is_directed())
    node_index = eg.new_vertex_property("int")
    half_edges = g.new_vertex_property("vector<int>")
    be = eg.new_vertex_property("int")

    # create half-edge graph
    libcommunity.get_eg_overlap(g._Graph__graph,
                                eg._Graph__graph,
                                _prop("e", g, b),
                                _prop("v", eg, be),
                                _prop("v", eg, node_index),
                                _prop("v", g, half_edges))

    if b_array is not None:
        be.a = b_array

    return eg, be, node_index, half_edges

def augmented_graph(g, b, node_index, eweight=None):
    r"""Generates an augmented graph from the half-edge graph ``g`` partitioned
    according to ``b``, where each half-edge belonging to a different group
    inside each node forms a new node."""

    node_map = g.new_vertex_property("int")
    br_b = libcore.Vector_int32_t()
    br_ni = libcore.Vector_int32_t()
    libcommunity.get_augmented_overlap(g._Graph__graph,
                                       _prop("v", g, b),
                                       _prop("v", g, node_index),
                                       _prop("v", g, node_map),
                                       br_b, br_ni)


    au, idx, vcount, ecount = condensation_graph(g, node_map,
                                                 eweight=eweight,
                                                 self_loops=True)[:4]
    anidx = idx.copy("int")
    libcommunity.vector_map(anidx.a, br_ni.a)

    ab = idx.copy("int")
    libcommunity.vector_map(ab.a, br_b.a)

    return au, ab, anidx, ecount, node_map

def get_block_edge_gradient(g, be, cmap=None):
    r"""Get edge gradients corresponding to the block membership at the endpoints of
    the edges given by the ``be`` edge property map.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        The graph.
    be : :class:`~graph_tool.PropertyMap`
        Vector-valued edge property map with the block membership at each
        endpoint.
    cmap : :class:`matplotlib.colors.Colormap` (optional, default: ``default_cm``)
        Color map used to construct the gradient.

    Returns
    -------
    cp : :class:`~graph_tool.PropertyMap`
       A vector-valued edge property map containing a color gradient.
    """

    if cmap is None:
        from .. draw import default_cm
        cmap = default_cm

    cp = g.new_edge_property("vector<double>")
    rg = [float("inf"), -float("inf")]
    for e in g.edges():
        s, t = be[e]
        rg[0] = min(s, rg[0])
        rg[0] = min(t, rg[0])
        rg[1] = max(s, rg[1])
        rg[1] = max(t, rg[1])

    for e in g.edges():
        if int(e.source()) < int(e.target()) or g.is_directed():
            s, t = be[e]
        else:
            t, s = be[e]
        cs = cmap((s - rg[0]) / max(rg[1] - rg[0], 1))
        ct = cmap((t - rg[0]) / max(rg[1] - rg[0], 1))
        cp[e] = [0] + list(cs) + [1] + list(ct)
    return cp
