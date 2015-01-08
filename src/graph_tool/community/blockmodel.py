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
from .. stats import label_self_loops
from .. spectral import adjacency
import random
from numpy import *
import numpy
from scipy.optimize import fsolve, fminbound
import scipy.special
from collections import defaultdict
import copy
import heapq

from .. dl_import import dl_import
dl_import("from . import libgraph_tool_community as libcommunity")

__test__ = False

def get_block_graph(g, B, b, vcount, ecount):
    cg, br, vcount, ecount = condensation_graph(g, b,
                                                vweight=vcount,
                                                eweight=ecount,
                                                self_loops=True)[:4]
    cg.vp["count"] = vcount
    cg.ep["count"] = ecount
    cg = Graph(cg, vorder=br)

    cg.add_vertex(B - cg.num_vertices())
    return cg

class BlockState(object):
    r"""This class encapsulates the block state of a given graph.

    This must be instantiated and used by functions such as :func:`mcmc_sweep`.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be modelled.
    eweight : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Edge multiplicities (for multigraphs or block graphs).
    vweight : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Vertex multiplicities (for block graphs).
    b : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Initial block labels on the vertices. If not supplied, it will be
        randomly sampled.
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

    _state_ref_count = 0

    def __init__(self, g, eweight=None, vweight=None, b=None,
                 B=None, clabel=None, deg_corr=True,
                 max_BE=1000, **kwargs):

        BlockState._state_ref_count += 1

        # initialize weights to unity, if necessary
        if eweight is None:
            eweight = g.new_edge_property("int")
            eweight.fa = 1
        elif eweight.value_type() != "int32_t":
            eweight = eweight.copy(value_type="int32_t")
        if vweight is None:
            vweight = g.new_vertex_property("int")
            vweight.fa = 1
        elif vweight.value_type() != "int32_t":
            vweight = vweight.copy(value_type="int32_t")
        self.eweight = g.own_property(eweight)
        self.vweight = g.own_property(vweight)

        self.is_weighted = True if self.eweight.fa.max() > 1 else False

        # configure the main graph and block model parameters
        self.g = g

        self.E = int(self.eweight.fa.sum())
        self.N = int(self.vweight.fa.sum())

        self.deg_corr = deg_corr

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
            if isinstance(b, numpy.ndarray):
                self.b = g.new_vertex_property("int")
                self.b.fa = b
            else:
                self.b = b = g.own_property(b.copy(value_type="int"))
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

        if clabel is not None:
            if isinstance(clabel, PropertyMap):
                self.clabel = clabel
            else:
                self.clabel = self.g.new_vertex_property("int")
                self.clabel.a = clabel
        else:
            self.clabel = self.g.new_vertex_property("int")

        self.emat = None
        if max_BE is None:
            max_BE = 1000
        self.max_BE = max_BE

        self.overlap = False

        # used by mcmc_sweep()
        self.egroups = None
        self.nsampler = None
        self.sweep_vertices = None
        self.overlap_stats = libcommunity.overlap_stats()
        self.partition_stats = libcommunity.partition_stats()

        # computation cache
        libcommunity.init_safelog(int(5 * max(self.E, self.N)))
        libcommunity.init_xlogx(int(5 * max(self.E, self.N)))
        libcommunity.init_lgamma(int(3 * max(self.E, self.N)))

    def __del__(self):
        try:
            BlockState._state_ref_count -= 1
            if BlockState._state_ref_count == 0:
                libcommunity.clear_safelog()
                libcommunity.clear_xlogx()
                libcommunity.clear_lgamma()
        except (ValueError, AttributeError, TypeError):
            pass

    def __repr__(self):
        return "<BlockState object with %d blocks,%s for graph %s, at 0x%x>" % \
            (self.B, " degree corrected," if self.deg_corr else "", str(self.g),
             id(self))


    def __init_partition_stats(self, empty=True):
        if not empty:
            self.partition_stats = libcommunity.init_partition_stats(self.g._Graph__graph,
                                                                     _prop("v", self.g, self.b),
                                                                     _prop("e", self.g, self.eweight),
                                                                     self.N, self.B)
        else:
            self.partition_stats = libcommunity.partition_stats()



    def copy(self, b=None, B=None, deg_corr=None, clabel=None, overlap=False):
        r"""Copies the block state. The parameters override the state properties, and
         have the same meaning as in the constructor. If ``overlap=True`` an
         instance of :class:`~graph_tool.community.OverlapBlockState` is
         returned."""

        if not overlap:
            state = BlockState(self.g,
                               eweight=self.eweight,
                               vweight=self.vweight,
                               b=self.b.copy() if b is None else b,
                               B=(self.B if b is None else None) if B is None else B,
                               clabel=self.clabel if clabel is None else clabel,
                               deg_corr=self.deg_corr if deg_corr is None else deg_corr,
                               max_BE=self.max_BE)
        else:
            state = OverlapBlockState(self.g,
                                      b=b if b is not None else self.b,
                                      B=(self.B if b is None else None) if B is None else B,
                                      clabel=self.clabel if clabel is None else clabel,
                                      deg_corr=self.deg_corr if deg_corr is None else deg_corr,
                                      max_BE=self.max_BE)

        if not state.__check_clabel():
            b = state.b.a + state.clabel.a * state.B
            continuous_map(b)
            state = state.copy(b=b)

            if __test__:
                assert state.__check_clabel()

        return state


    def __getstate__(self):
        state = dict(g=self.g,
                     eweight=self.eweight,
                     vweight=self.vweight,
                     b=self.b,
                     B=self.B,
                     clabel=self.clabel,
                     deg_corr=self.deg_corr,
                     max_BE=self.max_BE)
        return state

    def __setstate__(self, state):
        self.__init__(**state)
        return state

    def get_block_state(self, b=None, vweight=False, deg_corr=False, overlap=False):
        r"""Returns a :class:`~graph_tool.community.BlockState`` corresponding to the
        block graph. The parameters have the same meaning as the in the constructor."""


        state = BlockState(self.bg, eweight=self.mrs,
                           vweight=self.wr if vweight else None,
                           b=self.bg.vertex_index.copy("int") if b is None else b,
                           clabel=self.get_bclabel(),
                           deg_corr=deg_corr,
                           max_BE=self.max_BE)
        if overlap:
            state = state.copy(overlap=True)
        n_map = self.b.copy()
        return state, n_map

    def get_bclabel(self):
        r"""Returns a :class:`~graph_tool.PropertyMap`` corresponding to constraint
        labels for the block graph."""

        bclabel = self.bg.new_vertex_property("int")
        reverse_map(self.b, bclabel)
        pmap(bclabel, self.clabel)
        return bclabel

    def __check_clabel(self):
        b = self.b.a + self.clabel.a * self.B
        continuous_map(b)
        b2 = self.b.copy()
        continuous_map(b2.a)
        return (b == b2.a).all()

    def __get_emat(self):
        if self.emat is None:
            self.__regen_emat()
        return self.emat

    def __regen_emat(self):
        if self.B <= self.max_BE:
            self.emat = libcommunity.create_emat(self.bg._Graph__graph)
        else:
            self.emat = libcommunity.create_ehash(self.bg._Graph__graph)

    def __build_egroups(self, empty=False):
        self.esrcpos = self.g.new_edge_property("int")
        self.etgtpos = self.g.new_edge_property("int")

        self.egroups = libcommunity.build_egroups(self.g._Graph__graph,
                                                  self.bg._Graph__graph,
                                                  _prop("v", self.g, self.b),
                                                  _prop("e", self.g, self.eweight),
                                                  _prop("e", self.g, self.esrcpos),
                                                  _prop("e", self.g, self.etgtpos),
                                                  self.is_weighted, empty)

    def __build_nsampler(self, empty=False):
        self.nsampler = libcommunity.init_neighbour_sampler(self.g._Graph__graph,
                                                            _prop("e", self.g, self.eweight),
                                                            True, empty)

    def __cleanup_bg(self):
        emask = self.bg.new_edge_property("bool")
        emask.a = self.mrs.a[:len(emask.a)] > 0
        self.bg.set_edge_filter(emask)
        self.bg.purge_edges()
        self.emat = None

    def get_blocks(self):
        r"""Returns the property map which contains the block labels for each vertex."""
        return self.b

    def get_bg(self):
        r"""Returns the block graph."""
        return self.bg

    def get_ers(self):
        r"""Returns the edge property map of the block graph which contains the :math:`e_{rs}` matrix entries.
        For undirected graphs, the diagonal values (self-loops) contain :math:`e_{rr}/2`."""
        return self.mrs

    def get_er(self):
        r"""Returns the vertex property map of the block graph which contains the number
        :math:`e_r` of half-edges incident on block :math:`r`. If the graph is
        directed, a pair of property maps is returned, with the number of
        out-edges :math:`e^+_r` and in-edges :math:`e^-_r`, respectively."""
        if self.bg.is_directed():
            return self.mrp, self.mrm
        else:
            return self.mrp

    def get_nr(self):
        r"""Returns the vertex property map of the block graph which contains the block sizes :math:`n_r`."""
        return self.wr

    def entropy(self, complete=True, dl=False, partition_dl=False, dense=False,
                multigraph=True, norm=True, dl_ent=False, **kwargs):
        r"""Calculate the entropy associated with the current block partition.

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
        dense : ``bool`` (optional, default: ``False``)
            If ``True``, the "dense" variant of the entropy will be computed.
        multigraph : ``bool`` (optional, default: ``False``)
            If ``True``, the multigraph entropy will be used.
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
        of half-edges for the undirected case when :math:`r=s`), and :math:`n_r`
        is the number of vertices in block :math:`r` .

        For the degree-corrected variant with "hard" degree constraints the
        equivalent expressions are

        .. math::

            \mathcal{S}_c &\cong -E -\sum_kN_k\ln k! - \frac{1}{2} \sum_{rs}e_{rs}\ln\left(\frac{e_{rs}}{e_re_s}\right), \\
            \mathcal{S}^d_c &\cong -E -\sum_{k^+}N_{k^+}\ln k^+!  -\sum_{k^-}N_{k^-}\ln k^-! - \sum_{rs}e_{rs}\ln\left(\frac{e_{rs}}{e^+_re^-_s}\right),

        where :math:`e_r = \sum_se_{rs}` is the number of half-edges incident on
        block :math:`r`, and :math:`e^+_r = \sum_se_{rs}` and :math:`e^-_r =
        \sum_se_{sr}` are the numbers of out- and in-edges adjacent to block
        :math:`r`, respectively.

        If ``dense == False`` and ``multigraph == True``, the entropy used will
        be of the "Poisson" model, with the additional term:

        .. math::

            {\mathcal{S}_{cm}^{(d)}} = \mathcal{S}_c^{(d)} + \sum_{i>j} \ln A_{ij}! + \sum_i \ln A_{ii}!!


        If ``dl == True``, the description length :math:`\mathcal{L}_t` of the
        model will be returned as well, as described in
        :func:`model_entropy`. Note that for the degree-corrected version the
        description length is

        .. math::

            \mathcal{L}_c = \mathcal{L}_t + \sum_r\min\left(\mathcal{L}^{(1)}_r, \mathcal{L}^{(2)}_r\right),

        with

        .. math::

              \mathcal{L}^{(1)}_r &= \ln{\left(\!\!{n_r \choose e_r}\!\!\right)}, \\
              \mathcal{L}^{(2)}_r &= \ln\Xi_r + \ln n_r! - \sum_k \ln n^r_k!,

        and :math:`\ln\Xi_r \simeq 2\sqrt{\zeta(2)e_r}`, where :math:`\zeta(x)`
        is the `Riemann zeta function
        <https://en.wikipedia.org/wiki/Riemann_zeta_function>`_, and
        :math:`n^r_k` is the number of nodes in block :math:`r` with degree
        :math:`k`. For directed graphs we have instead :math:`k \to (k^-, k^+)`,
        and :math:`\ln\Xi_r \to \ln\Xi^+_r + \ln\Xi^-_r` with :math:`\ln\Xi_r^+
        \simeq 2\sqrt{\zeta(2)e^+_r}` and :math:`\ln\Xi_r^- \simeq
        2\sqrt{\zeta(2)e^-_r}`.

        If ``dl_ent=True`` is passed, this will be approximated instead by

        .. math::

            \mathcal{L}_c \simeq \mathcal{L}_t - \sum_rn_r\sum_kp^r_k\ln p^r_k,

        where :math:`p^r_k = n^r_k / n_r`.

        If the "dense" entropies are requested (``dense=True``), they will be
        computed as

        .. math::

            \mathcal{S}_t  &= \sum_{r>s} \ln{\textstyle {n_rn_s \choose e_{rs}}} + \sum_r \ln{\textstyle {{n_r\choose 2}\choose e_{rr}/2}}\\
            \mathcal{S}^d_t  &= \sum_{rs} \ln{\textstyle {n_rn_s \choose e_{rs}}},

        for simple graphs, and

        .. math::

            \mathcal{S}_m  &= \sum_{r>s} \ln{\textstyle \left(\!\!{n_rn_s \choose e_{rs}}\!\!\right)} + \sum_r \ln{\textstyle \left(\!\!{\left(\!{n_r\choose 2}\!\right)\choose e_{rr}/2}\!\!\right)}\\
            \mathcal{S}^d_m  &= \sum_{rs} \ln{\textstyle \left(\!\!{n_rn_s \choose e_{rs}}\!\!\right)},

        for multigraphs (i.e. ``multigraph == True``). A dense entropy for the
        degree-corrected model is not available, and if requested will return a
        :exc:`NotImplementedError`.

        If ``complete == False`` constants in the above equations which do not
        depend on the partition of the nodes will be omitted.

        Note that in all cases if ``norm==True`` the value returned corresponds
        to the entropy `per edge`, i.e. :math:`(\mathcal{S}_{t/c}\; [\,+\,\mathcal{L}_{t/c}])/ E`.
        """

        xi_fast = kwargs.get("xi_fast", False)
        dl_deg_alt = kwargs.get("dl_deg_alt", True)

        E = self.E
        N = self.N

        if dense:
            if self.deg_corr:
                raise NotImplementedError('A degree-corrected "dense" entropy is not yet implemented')

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

            if __test__:
                assert not isnan(S) and not isinf(S), "invalid entropy %g (%s) " % (S, str(dict(complete=complete,
                                                                                                random=random, dl=dl,
                                                                                                partition_dl=partition_dl,
                                                                                                dense=dense, multigraph=multigraph,
                                                                                                norm=norm)))
            if complete:
                if self.deg_corr:
                    S -= E
                    S += libcommunity.deg_entropy_term(self.g._Graph__graph,
                                                       libcore.any(),
                                                       self.overlap_stats,
                                                       self.N)
                else:
                    S += E

                if multigraph:
                    S += libcommunity.entropy_parallel(self.g._Graph__graph,
                                                       _prop("e", self.g, self.eweight))

                if __test__:
                    assert not isnan(S) and not isinf(S), "invalid entropy %g (%s) " % (S, str(dict(complete=complete,
                                                                                                    random=random, dl=dl,
                                                                                                    partition_dl=partition_dl,
                                                                                                    dense=dense, multigraph=multigraph,
                                                                                                    norm=norm)))
        if dl:
            if partition_dl:
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
            else:
                S += model_entropy(self.B, N, E, directed=self.g.is_directed(), nr=self.wr.a) * E

            if self.deg_corr:
                if self.partition_stats.is_enabled():
                    S_seq = self.partition_stats.get_deg_dl(dl_ent, dl_deg_alt, xi_fast)
                else:
                    self.__init_partition_stats(empty=False)
                    S_seq = self.partition_stats.get_deg_dl(dl_ent, dl_deg_alt, xi_fast)
                    self.__init_partition_stats(empty=True)

                S += S_seq

                if __test__:
                    assert not isnan(S_seq) and not isinf(S_seq), "invalid entropy %g (%s) " % (S_seq, str(dict(complete=complete,
                                                                                                                random=random, dl=dl,
                                                                                                                partition_dl=partition_dl,
                                                                                                                dense=dense, multigraph=multigraph,
                                                                                                                norm=norm)))

        if __test__:
            assert not isnan(S) and not isinf(S), "invalid entropy %g (%s) " % (S, str(dict(complete=complete,
                                                                                            random=random, dl=dl,
                                                                                            partition_dl=partition_dl,
                                                                                            dense=dense, multigraph=multigraph,
                                                                                            norm=norm)))

        if norm:
            return S / E
        else:
            return S

    def get_matrix(self):
        r"""Returns the block matrix (as a sparse :class:`~scipy.sparse.csr_matrix`),
        which contains the number of edges between each block pair.

        Examples
        --------

        .. testsetup:: get_matrix

           gt.seed_rng(42)
           np.random.seed(42)
           from pylab import *

        .. doctest:: get_matrix

           >>> g = gt.collection.data["polbooks"]
           >>> state = gt.BlockState(g, B=5, deg_corr=True)
           >>> for i in range(1000):
           ...     ds, nmoves = gt.mcmc_sweep(state)
           >>> m = state.get_matrix()
           >>> figure()
           <...>
           >>> matshow(m.todense())
           <...>
           >>> savefig("bloc_mat.pdf")

        .. testcleanup:: get_matrix

           savefig("bloc_mat.png")

        .. figure:: bloc_mat.*
           :align: center

           A  5x5 block matrix.

       """

        return adjacency(self.bg, weight=self.mrs)


def model_entropy(B, N, E, directed=False, nr=None):
    r"""Computes the amount of information necessary for the parameters of the traditional blockmodel ensemble, for ``B`` blocks, ``N`` vertices, ``E`` edges, and either a directed or undirected graph.

    A traditional blockmodel is defined as a set of :math:`N` vertices which can
    belong to one of :math:`B` blocks, and the matrix :math:`e_{rs}` describes
    the number of edges from block :math:`r` to :math:`s` (or twice that number
    if :math:`r=s` and the graph is undirected).

    For an undirected graph, the number of distinct :math:`e_{rs}` matrices is given by,

    .. math::

       \Omega_m = \left(\!\!{\left(\!{B \choose 2}\!\right) \choose E}\!\!\right)

    and for a directed graph,

    .. math::
       \Omega_m = \left(\!\!{B^2 \choose E}\!\!\right)


    where :math:`\left(\!{n \choose k}\!\right) = {n+k-1\choose k}` is the
    number of :math:`k` combinations with repetitions from a set of size :math:`n`.

    The total information necessary to describe the model is then,

    .. math::

       \mathcal{L}_t = \ln\Omega_m + \ln\left(\!\!{B \choose N}\!\!\right) + \ln N! - \sum_r \ln n_r!,


    where the remaining term is the information necessary to describe the
    block partition, where :math:`n_r` is the number of nodes in block :math:`r`.

    If ``nr`` is ``None``, it is assumed :math:`n_r=N/B`.

    The value returned corresponds to the information per edge, i.e.
    :math:`\mathcal{L}_t/E`.

    References
    ----------

    .. [peixoto-parsimonious-2013] Tiago P. Peixoto, "Parsimonious module inference in large networks",
       Phys. Rev. Lett. 110, 148701 (2013), :doi:`10.1103/PhysRevLett.110.148701`, :arxiv:`1212.4794`.
    .. [peixoto-hierarchical-2014] Tiago P. Peixoto, "Hierarchical block structures and high-resolution
       model selection in large networks ", Phys. Rev. X 4, 011047 (2014), :doi:`10.1103/PhysRevX.4.011047`,
       :arxiv:`1310.4377`.

    """

    if directed:
        x = (B * B);
    else:
        x = (B * (B + 1)) / 2;
    L = lbinom(x + E - 1, E) + partition_entropy(B, N, nr)
    return L / E

def Sdl(B, S, N, E, directed=False):
    return S + model_entropy(B, N, E, directed)

def lbinom(n, k):
    return scipy.special.gammaln(float(n + 1)) - scipy.special.gammaln(float(n - k + 1)) - scipy.special.gammaln(float(k + 1))

def partition_entropy(B, N, nr=None):
    if nr is None:
        S = N * log(B) + log1p(-(1 - 1./B) ** N)
    else:
        S = lbinom(B + N - 1, N) + scipy.special.gammaln(N + 1) - scipy.special.gammaln(nr + 1).sum()
    return S

def get_max_B(N, E, directed=False):
    r"""Return the maximum detectable number of blocks, obtained by minimizing:

    .. math::

        \mathcal{L}_t(B, N, E) - E\ln B

    where :math:`\mathcal{L}_t(B, N, E)` is the information necessary to
    describe a traditional blockmodel with `B` blocks, `N` nodes and `E`
    edges (see :func:`model_entropy`).

    Examples
    --------

    >>> gt.get_max_B(N=1e6, E=5e6)
    1572

    References
    ----------
    .. [peixoto-parsimonious-2013] Tiago P. Peixoto, "Parsimonious module inference in large networks",
       Phys. Rev. Lett. 110, 148701 (2013), :doi:`10.1103/PhysRevLett.110.148701`, :arxiv:`1212.4794`.


    """

    B = fminbound(lambda B: Sdl(B, -log(B), N, E, directed), 1, E,
                  xtol=1e-6, maxfun=1500, disp=0)
    if isnan(B):
        B = 1
    return max(int(ceil(B)), 2)

def get_akc(B, I, N=float("inf"), directed=False):
    r"""Return the minimum value of the average degree of the network, so that some block structure with :math:`B` blocks can be detected, according to the minimum description length criterion.

    This is obtained by solving

    .. math::

       \Sigma_b = \mathcal{L}_t(B, N, E) - E\mathcal{I}_{t/c} = 0,

    where :math:`\mathcal{L}_{t}` is the necessary information to describe the
    blockmodel parameters (see :func:`model_entropy`), and
    :math:`\mathcal{I}_{t/c}` characterizes the planted block structure, and is
    given by

    .. math::

        \mathcal{I}_t &= \sum_{rs}m_{rs}\ln\left(\frac{m_{rs}}{w_rw_s}\right),\\
        \mathcal{I}_c &= \sum_{rs}m_{rs}\ln\left(\frac{m_{rs}}{m_rm_s}\right),

    where :math:`m_{rs} = e_{rs}/2E` (or :math:`m_{rs} = e_{rs}/E` for directed
    graphs) and :math:`w_r=n_r/N`. We note that :math:`\mathcal{I}_{t/c}\in[0,
    \ln B]`. In the case where :math:`E \gg B^2`, this simplifies to

    .. math::

       \left<k\right>_c &= \frac{2\ln B}{\mathcal{I}_{t/c}},\\
       \left<k^{-/+}\right>_c &= \frac{\ln B}{\mathcal{I}_{t/c}},

    for undirected and directed graphs, respectively. This limit is assumed if
    ``N == inf``.

    Examples
    --------

    >>> gt.get_akc(10, log(10) / 100, N=100)
    2.414413200430159

    References
    ----------
    .. [peixoto-parsimonious-2013] Tiago P. Peixoto, "Parsimonious module inference in large networks",
       Phys. Rev. Lett. 110, 148701 (2013), :doi:`10.1103/PhysRevLett.110.148701`, :arxiv:`1212.4794`.

    """
    if N != float("inf"):
        if directed:
            get_dl = lambda ak: model_entropy(B, N, N * ak, directed) - N * ak * I
        else:
            get_dl = lambda ak: model_entropy(B, N, N * ak / 2., directed) - N * ak * I / 2.
        ak = fsolve(lambda ak: get_dl(ak), 10)
        ak = float(ak)
    else:
        ak = 2 * log(B) / S
        if directed:
            ak /= 2
    return ak

def mcmc_sweep(state, beta=1., c=1., dl=False, dense=False, multigraph=False,
               node_coherent=False, nmerges=0, nmerge_sweeps=1, merge_map=None,
               coherent_merge=False, sequential=True, parallel=False,
               vertices=None, verbose=False, **kwargs):
    r"""Performs a Markov chain Monte Carlo sweep on the network, to sample the block partition according to a probability :math:`\propto e^{-\beta \mathcal{S}_{t/c}}`, where :math:`\mathcal{S}_{t/c}` is the blockmodel entropy.

    Parameters
    ----------
    state : :class:`~graph_tool.community.BlockState` or :class:`~graph_tool.community.OverlapBlockState`
        The block state.
    beta : ``float`` (optional, default: `1.0`)
        The inverse temperature parameter :math:`\beta`.
    c : ``float`` (optional, default: ``1.0``)
        This parameter specifies how often fully random moves are attempted,
        instead of more likely moves based on the inferred block partition.
        For ``c == 0``, no fully random moves are attempted, and for ``c == inf``
        they are always attempted.
    dl : ``bool`` (optional, default: ``False``)
        If ``True``, the change in the whole description length will be
        considered after each vertex move, not only the entropy.
    dense : ``bool`` (optional, default: ``False``)
        If ``True``, the "dense" variant of the entropy will be computed.
    multigraph : ``bool`` (optional, default: ``False``)
        If ``True``, the multigraph entropy will be used. Only has an effect
        if ``dense == True``.
    node_coherent : ``bool`` (optional, default: ``False``)
        If ``True``, and if the ``state`` is an instance of
        :class:`~graph_tool.community.OverlapBlockState`, then all half-edges
        incident on the same node are moved simultaneously.
    sequential : ``bool`` (optional, default: ``True``)
        If ``True``, the move attempts on the vertices are done in sequential
        random order. Otherwise a total of `N` moves attempts are made, where
        `N` is the number of vertices, where each vertex can be selected with
        equal probability.
    parallel : ``bool`` (optional, default: ``False``)
        If ``True``, the updates are performed in parallel (multiple
        threads).

        .. warning::

            If this is used, the Markov Chain is not guaranteed to be sampled with
            the correct probabilities. This is better used in conjunction with
            ``beta=float('inf')``, where this is not an issue.

    vertices : ``list of ints`` (optional, default: ``None``)
        A list of vertices which will be attempted to be moved. If ``None``, all
        vertices will be attempted.
    verbose : ``bool`` (optional, default: ``False``)
        If ``True``, verbose information is displayed.

    Returns
    -------

    dS : ``float``
       The entropy difference (per edge) after a full sweep.
    nmoves : ``int``
       The number of accepted block membership moves.


    Notes
    -----

    This algorithm performs a Markov chain Monte Carlo sweep on the network,
    where the block memberships are randomly moved, and either accepted or
    rejected, so that after sufficiently many sweeps the partitions are sampled
    with probability proportional to :math:`e^{-\beta\mathcal{S}_{t/c}}`, where
    :math:`\mathcal{S}_{t/c}` is the blockmodel entropy, given by

    .. math::

      \mathcal{S}_t &\cong - \frac{1}{2} \sum_{rs}e_{rs}\ln\left(\frac{e_{rs}}{n_rn_s}\right), \\
      \mathcal{S}^d_t &\cong - \sum_{rs}e_{rs}\ln\left(\frac{e_{rs}}{n_rn_s}\right),

    for undirected and directed traditional blockmodels (``deg_corr == False``),
    respectively, where :math:`e_{rs}` is the number of edges from block
    :math:`r` to :math:`s` (or the number of half-edges for the undirected case
    when :math:`r=s`), and :math:`n_r` is the number of vertices in block
    :math:`r`, and constant terms which are independent of the block partition
    were dropped (see :meth:`BlockState.entropy` for the complete entropy). For
    the degree-corrected variant with "hard" degree constraints the equivalent
    expressions are

    .. math::

       \mathcal{S}_c &\cong  - \frac{1}{2} \sum_{rs}e_{rs}\ln\left(\frac{e_{rs}}{e_re_s}\right), \\
       \mathcal{S}^d_c &\cong - \sum_{rs}e_{rs}\ln\left(\frac{e_{rs}}{e^+_re^-_s}\right),

    where :math:`e_r = \sum_se_{rs}` is the number of half-edges incident on
    block :math:`r`, and :math:`e^+_r = \sum_se_{rs}` and :math:`e^-_r =
    \sum_se_{sr}` are the number of out- and in-edges adjacent to block
    :math:`r`, respectively.

    The Monte Carlo algorithm employed attempts to improve the mixing time of
    the Markov chain by proposing membership moves :math:`r\to s` with
    probability :math:`p(r\to s|t) \propto e_{ts} + c`, where :math:`t` is the
    block label of a random neighbour of the vertex being moved. See
    [peixoto-efficient-2014]_ for more details.

    This algorithm has a complexity of :math:`O(E)`, where :math:`E` is the
    number of edges in the network.

    Examples
    --------
    .. testsetup:: mcmc

       gt.seed_rng(42)
       np.random.seed(42)

    .. doctest:: mcmc

       >>> g = gt.collection.data["polbooks"]
       >>> state = gt.BlockState(g, B=3, deg_corr=True)
       >>> pv = None
       >>> for i in range(1000):        # remove part of the transient
       ...     ds, nmoves = gt.mcmc_sweep(state)
       >>> for i in range(1000):
       ...     ds, nmoves = gt.mcmc_sweep(state)
       ...     pv = gt.collect_vertex_marginals(state, pv)
       >>> gt.graph_draw(g, pos=g.vp["pos"], vertex_shape="pie", vertex_pie_fractions=pv, output="polbooks_blocks_soft.pdf")
       <...>

    .. testcleanup:: mcmc

       gt.graph_draw(g, pos=g.vp["pos"], vertex_shape="pie", vertex_pie_fractions=pv, output="polbooks_blocks_soft.png")

    .. figure:: polbooks_blocks_soft.*
       :align: center

       "Soft" block partition of a political books network with :math:`B=3`.

     References
    ----------

    .. [holland-stochastic-1983] Paul W. Holland, Kathryn Blackmond Laskey,
       Samuel Leinhardt, "Stochastic blockmodels: First steps",
       Carnegie-Mellon University, Pittsburgh, PA 15213, U.S.A., :doi:`10.1016/0378-8733(83)90021-7`
    .. [faust-blockmodels-1992] Katherine Faust, and Stanley
       Wasserman. "Blockmodels: Interpretation and Evaluation." Social Networks
       14, no. 1-2 (1992): 5-61. :doi:`10.1016/0378-8733(92)90013-W`
    .. [karrer-stochastic-2011] Brian Karrer, and M. E. J. Newman. "Stochastic
       Blockmodels and Community Structure in Networks." Physical Review E 83,
       no. 1 (2011): 016107. :doi:`10.1103/PhysRevE.83.016107`.
    .. [peixoto-entropy-2012] Tiago P. Peixoto "Entropy of Stochastic Blockmodel
       Ensembles." Physical Review E 85, no. 5 (2012): 056122. :doi:`10.1103/PhysRevE.85.056122`,
       :arxiv:`1112.6028`.
    .. [peixoto-parsimonious-2013] Tiago P. Peixoto, "Parsimonious module inference in large networks",
       Phys. Rev. Lett. 110, 148701 (2013), :doi:`10.1103/PhysRevLett.110.148701`, :arxiv:`1212.4794`.
    .. [peixoto-efficient-2014] Tiago P. Peixoto, "Efficient Monte Carlo and greedy
       heuristic for the inference of stochastic block models", Phys. Rev. E 89, 012804 (2014),
       :doi:`10.1103/PhysRevE.89.012804`, :arxiv:`1310.4378`.
    .. [peixoto-model-2014] Tiago P. Peixoto, "Model selection and hypothesis
       testing for large-scale network models with overlapping groups",
       :arxiv:`1409.3059`.
    """

    if state.B == 1:
        return 0., 0

    if vertices is not None:
        vlist = libcommunity.get_vector(len(vertices))
        vlist.a = vertices
        vertices = vlist
        #state.sweep_vertices = vertices

    if state.sweep_vertices is None:
        vertices = libcommunity.get_vector(state.g.num_vertices())
        vertices.a = state.g.vertex_index.copy("int").fa
        state.sweep_vertices = vertices

    random_move = c == float("inf")

    if random_move or nmerges > 0:
        state._BlockState__build_egroups(empty=True)
    elif state.egroups is None:
        state._BlockState__build_egroups(empty=False)

    bclabel = state.get_bclabel()

    if nmerges == 0:
        nmerge_sweeps = 1
        if state.nsampler is None:
            state._BlockState__build_nsampler(empty=state.overlap)
        nsampler = state.nsampler
        ncavity_sampler = state.nsampler
        merge_map = state.g.new_vertex_property("int")
    else:
        if kwargs.get("unweighted_merge", False):
            emask = state.mrs
        else:
            emask = state.mrs.copy()
            emask.a = emask.a > 0

        nsampler = libcommunity.init_neighbour_sampler(state.bg._Graph__graph,
                                                       _prop("e", state.bg, emask),
                                                       True, False)
        ncavity_sampler = libcommunity.init_neighbour_sampler(state.bg._Graph__graph,
                                                              _prop("e", state.bg, emask),
                                                              False, False)
        beta = float("inf")

        if merge_map is None:
            merge_map = state.g.vertex_index.copy("int")

    if dl and not state.partition_stats.is_enabled():
        if state.overlap:
            state._OverlapBlockState__init_partition_stats(empty=False)
        else:
            state._BlockState__init_partition_stats(empty=False)

    if not dl and state.partition_stats.is_enabled():
        if state.overlap:
            state._OverlapBlockState__init_partition_stats(empty=True)
        else:
            state._BlockState__init_partition_stats(empty=True)

    if __test__:
        assert state._BlockState__check_clabel(), "clabel already invalid!"
        S = state.entropy(dense=dense, multigraph=multigraph, complete=False, dl=dl, dl_deg_alt=False, xi_fast=True)
        assert not (isinf(S) or isnan(S)), "invalid entropy before sweep: %g" % S

    try:
        if not state.overlap:
            dS, nmoves = libcommunity.move_sweep(state.g._Graph__graph,
                                                 state.bg._Graph__graph,
                                                 state._BlockState__get_emat(),
                                                 nsampler, ncavity_sampler,
                                                 _prop("e", state.bg, state.mrs),
                                                 _prop("v", state.bg, state.mrp),
                                                 _prop("v", state.bg, state.mrm),
                                                 _prop("v", state.bg, state.wr),
                                                 _prop("v", state.g, state.b),
                                                 _prop("v", state.bg, bclabel),
                                                 state.sweep_vertices,
                                                 state.deg_corr, dense, multigraph,
                                                 _prop("e", state.g, state.eweight),
                                                 _prop("v", state.g, state.vweight),
                                                 state.egroups,
                                                 _prop("e", state.g, state.esrcpos),
                                                 _prop("e", state.g, state.etgtpos),
                                                 float(beta), sequential,
                                                 parallel, random_move,
                                                 c, state.is_weighted,
                                                 nmerges, nmerge_sweeps,
                                                 _prop("v", state.g, merge_map),
                                                 state.partition_stats,
                                                 verbose, _get_rng())
        else:
            dS, nmoves = libcommunity.move_sweep_overlap(state.g._Graph__graph,
                                                         state.bg._Graph__graph,
                                                         state._BlockState__get_emat(),
                                                         nsampler,
                                                         ncavity_sampler,
                                                         _prop("e", state.bg, state.mrs),
                                                         _prop("v", state.bg, state.mrp),
                                                         _prop("v", state.bg, state.mrm),
                                                         _prop("v", state.bg, state.wr),
                                                         _prop("v", state.g, state.b),
                                                         _prop("v", state.bg, bclabel),
                                                         state.sweep_vertices,
                                                         state.deg_corr, dense, multigraph,
                                                         multigraph,
                                                         _prop("e", state.g, state.eweight),
                                                         _prop("v", state.g, state.vweight),
                                                         state.egroups,
                                                         _prop("e", state.g, state.esrcpos),
                                                         _prop("e", state.g, state.etgtpos),
                                                         float(beta),
                                                         sequential, parallel,
                                                         random_move, float(c),
                                                         ((nmerges == 0 and node_coherent) or
                                                          (nmerges > 0 and coherent_merge)),
                                                         state.is_weighted,
                                                         nmerges, nmerge_sweeps,
                                                         _prop("v", state.g, merge_map),
                                                         state.overlap_stats,
                                                         state.partition_stats,
                                                         verbose, _get_rng())

    finally:
        if random_move:
            state.egroups = None
        if nmerges > 0:
            state.nsampler = None
            state.egroups = None

    if __test__:
        assert state._BlockState__check_clabel(), "clabel invalidated!"
        assert not (isinf(dS) or isnan(dS)), "invalid after sweep: %g" % dS
        S2 = state.entropy(dense=dense, multigraph=multigraph, complete=False, dl=dl, dl_deg_alt=False, xi_fast=True)
        c_dS = S2 - S
        if not abs(dS / state.E - c_dS) < 1e-6:
            print(dS / state.E, c_dS, nmoves, state.overlap, dense, multigraph, state.deg_corr, state.is_weighted, node_coherent)
        assert abs(dS / state.E - c_dS) < 1e-6, "invalid delta S (%g, %g)" % (dS / state.E, c_dS)

    return dS / state.E, nmoves


def pmap(prop, value_map):
    """Maps all the values of `prop` to the values given by `value_map`, which
    is indexed by the values of `prop`."""
    if isinstance(prop, PropertyMap):
        prop = prop.a
    if isinstance(value_map, PropertyMap):
        value_map = value_map.a
    if prop.max() >= len(value_map):
        raise ValueError("value map is not large enough!")
    libcommunity.vector_map(prop, value_map)

def reverse_map(prop, value_map):
    """Modify `value_map` such that the positions indexed by the values in `prop`
    correspond to their index in `prop`."""
    if isinstance(prop, PropertyMap):
        prop = prop.a
    if isinstance(value_map, PropertyMap):
        value_map = value_map.a
    if prop.max() >= len(value_map):
        raise ValueError("value map is not large enough!")
    libcommunity.vector_rmap(prop, value_map)

def continuous_map(prop):
    """Remap the values of ``prop`` in the continuous range :math:`[0, N-1]`."""
    if isinstance(prop, PropertyMap):
        prop = prop.a
    if prop.max() < len(prop):
        rmap = -ones(len(prop), dtype=prop.dtype)
        libcommunity.vector_map(prop, rmap)
    else:
        libcommunity.vector_continuous_map(prop)

def greedy_shrink(state, B, **kwargs):
    if B > state.B:
        raise ValueError("Cannot shrink to a larger size!")

    kwargs = kwargs.copy()
    if kwargs.get("nmerge_sweeps", None) is None:
        kwargs["nmerge_sweeps"] = max((2 * state.g.num_edges()) // state.g.num_vertices(), 1)

    verbose = kwargs.get("verbose", False)

    orig_state = state
    state = state.copy(B=state.B)

    # merge according to indirect neighbourhood; we put all group-nodes in their
    # own groups, and merge/move them until the desired size is reached
    curr_B = (state.wr.a > 0).sum()
    assert curr_B >= B, "shrinking to a larger size ?! (%d, %d)" % (curr_B, B)

    random = kwargs.get("random_move", False)
    old_state = state
    if not state.overlap:
        state, n_map = state.get_block_state(vweight=True, deg_corr=state.deg_corr)
        unilevel_minimize(state, **kwargs)
    merge_map = state.g.vertex_index.copy("int")

    unweighted = False
    kwargs["c"] = 0 if not random else float("inf")
    kwargs["dl"] = False
    while curr_B > B:
        dS, nmoves = mcmc_sweep(state, beta=float("inf"),
                                nmerges=curr_B - B,
                                merge_map=merge_map,
                                unweighted_merge=unweighted,
                                **kwargs)

        #assert nmoves == curr_B - (state.wr.a > 0).sum()
        curr_B = (state.wr.a > 0).sum()

        if verbose:
            print("merging, B=%d" % curr_B, "left:", curr_B - B,
                  "(%g, %d%s%s)" % (dS, nmoves, ", random" if random else "",
                                    ", unweighted" if unweighted else ""))

        if nmoves == 0:
            if not unweighted:
                unweighted = True
            else:
                kwargs["c"] = float("inf")
                random = True

    #assert curr_B == (state.wr.a > 0).sum()

    if not state.overlap:
        unilevel_minimize(state, **kwargs)  # block level moves
        pmap(merge_map, state.b)
        pmap(n_map, merge_map)
        continuous_map(n_map)
        state = orig_state.copy(b=n_map, B=B)
    else:
        pmap(merge_map, state.b)
        continuous_map(merge_map)
        state = orig_state.copy(b=merge_map, B=B)


    if __test__:
        assert state._BlockState__check_clabel(), "clabel already invalidated!"
        assert curr_B == (state.wr.a > 0).sum()
        curr_B = (state.wr.a > 0).sum()
        assert state.B == curr_B
        assert state.B == B

    return state


class MinimizeState(object):
    r"""This object stores information regarding the current entropy minimization
    state, so that the algorithms can resume previously started runs.
    This object can be saved to disk via the :mod:`pickle` interface."""

    def __init__(self):
        self.b_cache = {}
        self.checkpoint_state = defaultdict(dict)
        self.init = True

    def clear(self):
        r"""Clear state."""
        self.b_cache.clear()
        self.checkpoint_state.clear()

def unilevel_minimize(state, nsweeps=10, adaptive_sweeps=True, epsilon=0,
                      anneal=(1., 1.), greedy=True, c=0., dl=False, dense=False,
                      multigraph=True, sequential=True, parallel=False,
                      verbose=False, **kwargs):
    kwargs = kwargs.copy()
    kwargs.update(dict(c=c, dl=dl, dense=dense, multigraph=multigraph,
                       sequential=sequential, parallel=parallel))

    t_dS, t_nmoves = 0, 0

    S = state.entropy()

    if not adaptive_sweeps:
        ntotal = nsweeps if greedy else 2 * nsweeps
        if verbose:
            print("Performing %d sweeps for B=%d..." % (ntotal, state.B))

        for i in range(ntotal):
            if i < niter:
                continue
            if i < nsweeps and not greedy:
                beta = anneal[0]
            else:
                beta = float("inf")

            delta, nmoves = mcmc_sweep(state, beta=beta, **kwargs)

            S += delta
            t_dS += delta
            t_nmoves += nmoves
            niter += 1
    else:
        # adaptive mode
        min_dl = S
        max_dl = S
        count = 0
        bump = False
        beta =  anneal[0]
        last_min = min_dl
        greedy_step = greedy
        total_nmoves = 0

        if verbose and not greedy:
            print("Performing sweeps for beta = %g, B=%d (N=%d)..." % \
                   (beta, state.B, state.g.num_vertices()))

        eps = 1e-8
        niter = 0
        while True:
            if greedy_step:
                break
            if count > nsweeps:
                if not bump:
                    min_dl = max_dl = S
                    bump = True
                    count = 0
                else:
                    if anneal[1] <= 1 or min_dl == last_min:
                        break
                    else:
                        beta *= anneal[1]
                        count = 0
                        last_min = min_dl
                        if verbose:
                            print("Performing sweeps for beta = %g, B=%d (N=%d)..." % \
                                   (beta, state.B, state.g.num_vertices()))

            delta, nmoves = mcmc_sweep(state, beta=beta, **kwargs)

            if state.overlap and beta == float("inf"):
                ds, nm = mcmc_sweep(state, beta=beta, node_coherent=True, **kwargs)
                delta += ds
                nmoves += nm

            S += delta
            niter += 1
            total_nmoves += nmoves

            t_dS += delta
            t_nmoves += nmoves

            if S > max_dl + eps:
                max_dl = S
                count = 0
            elif S < min_dl - eps:
                min_dl = S
                count = 0
            else:
                count += 1

        if verbose:
            if not greedy_step:
                print("... performed %d sweeps with %d vertex moves" % (niter, total_nmoves))
            print("Performing sweeps for beta = ∞, B=%d (N=%d)..." % \
                  (state.B, state.g.num_vertices()))

        if not greedy_step:
            min_dl = S
            count = 0

        niter = 0
        total_nmoves = 0
        while count <= nsweeps:
            delta, nmoves = mcmc_sweep(state, beta=float("inf"), **kwargs)

            if state.overlap:
                ds, nm = mcmc_sweep(state, beta=float("inf"), node_coherent=True, **kwargs)
                delta += ds
                nmoves += nm

            S += delta
            niter += 1
            total_nmoves += nmoves

            t_dS += delta
            t_nmoves += nmoves

            if abs(delta) > eps and nmoves / state.g.num_vertices() > epsilon:
                min_dl = S
                count = 0
            else:
                count += 1

        if verbose:
            print("... performed %d sweeps with %d vertex moves" % (niter, total_nmoves))

        bi = state.b

    return t_dS, t_nmoves


def multilevel_minimize(state, B, nsweeps=10, adaptive_sweeps=True, epsilon=0,
                        anneal=(1., 1.), r=2., nmerge_sweeps=10, greedy=True,
                        c=0., dl=False, dense=False, multigraph=True,
                        sequential=True, parallel=False, checkpoint=None,
                        minimize_state=None, verbose=False, **kwargs):
    r"""Performs an agglomerative heuristic, which progressively merges blocks together (while allowing individual node moves) to achieve a good partition in ``B`` blocks.

    Parameters
    ----------
    state : :class:`~graph_tool.community.BlockState` or :class:`~graph_tool.community.OverlapBlockState`
        The block state.
    B : ``int``
        The desired number of blocks.
    nsweeps : ``int`` (optional, default: ``10``)
        The number of sweeps done after each merge step to reach the local
        minimum.
    adaptive_sweeps : ``bool`` (optional, default: ``True``)
        If ``True``, the number of sweeps necessary for the local minimum will
        be estimated to be enough so that no more than ``epsilon * N`` nodes
        changes their states in the last ``nsweeps`` sweeps.
    epsilon : ``float`` (optional, default: ``0``)
        Converge criterion for ``adaptive_sweeps``.
    anneal : pair of ``floats`` (optional, default: ``(1., 1.)``)
        The first value specifies the starting value for  ``beta`` of the MCMC
        steps, and the second value is the factor which is multiplied to ``beta``
        after each estimated equilibration (according to ``nsweeps`` and
        ``adaptive_sweeps``).
    r : ``float`` (optional, default: ``2.``)
        Agglomeration ratio for the merging steps. Each merge step will attempt
        to find the best partition into :math:`B_{i-1} / r` blocks, where
        :math:`B_{i-1}` is the number of blocks in the last step.
    nmerge_sweeps : `int` (optional, default: `10`)
        The number of merge sweeps done, where in each sweep a better merge
        candidate is searched for every block.
    greedy : ``bool`` (optional, default: ``True``)
        If ``True``, the value of ``beta`` of the MCMC steps are kept at
        infinity for all steps. Otherwise they change according to the ``anneal``
        parameter.
    c : ``float`` (optional, default: ``0.0``)
        This parameter specifies how often fully random moves are attempted,
        instead of more likely moves based on the inferred block partition.
        For ``c == 0``, no fully random moves are attempted, and for ``c == inf``
        they are always attempted.
    dl : ``bool`` (optional, default: ``False``)
        If ``True``, the change in the whole description length will be
        considered after each vertex move, not only the entropy.
    dense : ``bool`` (optional, default: ``False``)
        If ``True``, the "dense" variant of the entropy will be computed.
    multigraph : ``bool`` (optional, default: ``False``)
        If ``True``, the multigraph entropy will be used. Only has an effect
        if ``dense == True``.
    sequential : ``bool`` (optional, default: ``True``)
        If ``True``, the move attempts on the vertices are done in sequential
        random order. Otherwise a total of `N` moves attempts are made, where
        `N` is the number of vertices, where each vertex can be selected with
        equal probability.
    parallel : ``bool`` (optional, default: ``False``)
        If ``True``, the updates are performed in parallel (multiple
        threads).
    vertices: ``list of ints`` (optional, default: ``None``)
        A list of vertices which will be attempted to be moved. If ``None``, all
        vertices will be attempted.
    checkpoint : function (optional, default: ``None``)
        If provided, this function will be called after each call to
        :func:`mcmc_sweep`. This can be used to store the current state, so it
        can be continued later. The function must have the following signature:

        .. code-block:: python

            def checkpoint(state, S, delta, nmoves, minimize_state):
                ...

        where `state` is either a :class:`~graph_tool.community.BlockState`
        instance or ``None``, `S` is the current entropy value, `delta` is
        the entropy difference in the last MCMC sweep, and `nmoves` is the
        number of accepted block membership moves. The ``minimize_state``
        argument is a :class:`MinimizeState` instance which specifies the current
        state of the algorithm, which can be stored via :mod:`pickle`, and
        supplied via the ``minimize_state`` option below to continue from an
        interrupted run.

        This function will also be called when the MCMC has finished for the
        current value of :math:`B`, in which case ``state == None``, and the
        remaining parameters will be zero, except the last.
    minimize_state : :class:`MinimizeState` (optional, default: ``None``)
        If provided, this will specify an exact point of execution from which
        the algorithm will continue. The expected object is a :class:`MinimizeState`
        instance which will be passed to the callback of the ``checkpoint``
        option above, and  can be stored by :mod:`pickle`.
    verbose : ``bool`` (optional, default: ``False``)
        If ``True``, verbose information is displayed.

    Returns
    -------

    state : :class:`~graph_tool.community.BlockState`
        The new :class:`~graph_tool.community.BlockState` with ``B`` blocks.

    Notes
    -----

    This algorithm performs an agglomerative heuristic on the current block state,
    where blocks are progressively merged together, using repeated applications of
    the :func:`mcmc_sweep` moves, at different scales. See [peixoto-efficient-2014]_
    for more details.

    This algorithm has a complexity of :math:`O(N\ln^2 N)`, where :math:`N` is the
    number of nodes in the network.

    Examples
    --------
    .. testsetup:: multilevel_minimize

       gt.seed_rng(42)
       np.random.seed(42)

    .. doctest:: multilevel_minimize

       >>> g = gt.collection.data["polblogs"]
       >>> g = gt.GraphView(g, vfilt=gt.label_largest_component(gt.GraphView(g, directed=False)))
       >>> state = gt.BlockState(g, B=g.num_vertices(), deg_corr=True)
       >>> state = gt.multilevel_minimize(state, B=2)
       >>> gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=state.get_blocks(), output="polblogs_agg.pdf")
       <...>

    .. testcleanup:: multilevel_minimize

       gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=state.get_blocks(), output="polblogs_agg.png")

    .. figure:: polblogs_agg.*
       :align: center

       Block partition of a political blogs network with :math:`B=2`.

     References
    ----------

    .. [peixoto-efficient-2014] Tiago P. Peixoto, "Efficient Monte Carlo and greedy
       heuristic for the inference of stochastic block models", Phys. Rev. E 89, 012804 (2014),
       :doi:`10.1103/PhysRevE.89.012804`, :arxiv:`1310.4378`.
    """

    if minimize_state is None:
        minimize_state = MinimizeState()
    b_cache = minimize_state.b_cache
    checkpoint_state = minimize_state.checkpoint_state

    nkwargs = dict(nsweeps=nsweeps, epsilon=epsilon, c=c,
                   dl=dl, dense=dense, multigraph=multigraph,
                   nmerge_sweeps=nmerge_sweeps,
                   sequential=sequential, parallel=parallel)
    kwargs = copy.copy(kwargs)
    kwargs.update(nkwargs)

    nonoverlap_compare = kwargs.get("nonoverlap_compare", False)
    if "nonoverlap_compare" in kwargs:
        del kwargs["nonoverlap_compare"]
    assert not nonoverlap_compare, "don't do this"

    orig_state = state

    if __test__:
        assert state._BlockState__check_clabel(), "orig clabel already invalidated!"

    # some simple boundary conditions
    if B == 1:
        if state.clabel.fa.max() > 0:
            raise ValueError("Cannot shrink to B = 1 without invalidating constraints")
        bi = state.g.new_vertex_property("int")
        state = state.copy(b=bi)
        return state
    if B == state.g.num_vertices():
        bi = state.g.vertex_index.copy("int")
        state = state.copy(b=bi)
        return state

    Bi = state.B
    while True:
        # keep reducing B by a factor of "r", until desired size is reached
        Bi = max(int(round(Bi / r)), B)
        if Bi == state.B and Bi > B:
            Bi -= 1

        # check cache for previous results
        if b_cache is not None and Bi in b_cache:
            if __test__:
                assert (state.clabel.a == b_cache[Bi][1].clabel.a).all(), "wrong clabel in cache"
                assert state._BlockState__check_clabel(), "clabel already invalidated before cache"
                assert b_cache[Bi][1]._BlockState__check_clabel(), "clabel already invalidated after cache"
            state = b_cache[Bi][1].copy()

        if __test__:
            assert state._BlockState__check_clabel(), "clabel already invalidated!"

        # if necessary, shrink state
        if Bi < state.B:
            if verbose:
                print("Shrinking:", state.B, "->", Bi)

            state = greedy_shrink(state, B=Bi, verbose=verbose, **kwargs)

            if __test__:
                assert state._BlockState__check_clabel(), "clabel invalidated after shrink"

        dS, nmoves = unilevel_minimize(state, checkpoint=checkpoint, verbose=verbose, **kwargs)

        if __test__:
            assert state._BlockState__check_clabel(), "clabel invalidated after unilevel minimize!"

        if state.overlap and state.deg_corr and nonoverlap_compare:
            if verbose:
                print("Attempting nonoverlapping minimize...")
            nstate = state.copy(b=state.get_nonoverlap_blocks(), overlap=False)
            assert nstate.B <= nstate.N
            nstate = multilevel_minimize(nstate, B=Bi, verbose=verbose, **kwargs)
            nstate = nstate.copy(overlap=True, clabel=state.clabel.a)
            unilevel_minimize(nstate, **kwargs)

            if nstate.B > Bi:
                nstate = multilevel_minimize(nstate, B=Bi, verbose=verbose,
                                             nonoverlap_compare=False, **kwargs)

            if nstate.entropy(dense=dense, multigraph=multigraph) < state.entropy(dense=dense, multigraph=multigraph):
                if verbose:
                    print("Nonoverlapping minimize improved.")
                state = nstate

                if __test__:
                    assert state._BlockState__check_clabel(), "clabel invalidated after nonoverlap compare!"

        if Bi == B:
            break

    return state

def get_b_dl(state, dense, multigraph, nested_dl, complete=False,
             nested_overlap=False, dl_ent=False):
    if not nested_dl:
        dl = state.entropy(dense=dense, multigraph=multigraph, dl=True,
                           complete=complete, dl_ent=dl_ent)
    else:
        dl = state.entropy(dense=dense, multigraph=multigraph, dl=True,
                           partition_dl=True, complete=complete,
                           dl_ent=dl_ent)

        bclabel = state.get_bclabel()

        bstate = state.get_block_state(b=bclabel, overlap=nested_overlap)[0]

        dl += bstate.entropy(dl=True, partition_dl=True, dense=True,
                             multigraph=True, dl_ent=dl_ent)
    return dl


def get_state_dl(B, minimize_state, checkpoint, sparse_heuristic, **kwargs):
    bs = minimize_state.b_cache
    checkpoint_state = minimize_state.checkpoint_state
    previous = None
    verbose = kwargs.get("verbose", False)

    if B in bs and checkpoint_state[B].get("done", False):
        # A previous finished result is available. Use that and keep going.
        if verbose:
            print("(using previous finished result for B=%d)" % B)
        if __test__:
            dl = get_b_dl(bs[B][1],
                          kwargs.get("dense", False),
                          kwargs.get("multigraph", False),
                          kwargs.get("nested_dl", False),
                          kwargs.get("complete", False),
                          kwargs.get("nested_overlap", False),
                          kwargs.get("dl_ent", False))
            assert abs(dl - bs[B][0]) < 1e-8, "inconsistent DL values! (%g, %g, overlap: %s)" % (dl, bs[B][0], str(bs[B][1].overlap))
        return bs[B][0]
    elif B in bs:
        # A previous unfinished result is available. Use that as the starting point.
        if verbose:
            print("(starting from previous result for B=%d)" % B)
        previous = bs[B]
        state = previous[1].copy()

        if __test__:
            assert state._BlockState__check_clabel(), "previous clabel already invalidated!"
            dl = get_b_dl(state, kwargs.get("dense", False),
                          kwargs.get("multigraph", False),
                          kwargs.get("nested_dl", False),
                          kwargs.get("complete", False),
                          kwargs.get("nested_overlap", False),
                          kwargs.get("dl_ent", False))
            assert abs(dl - bs[B][0]) < 1e-8, "inconsistent DL values! (%g, %g)" % (dl, bs[B][0])
            dl = get_b_dl(previous[1], kwargs.get("dense", False),
                          kwargs.get("multigraph", False),
                          kwargs.get("nested_dl", False),
                          kwargs.get("complete", False),
                          kwargs.get("nested_overlap", False),
                          kwargs.get("dl_ent", False))
            assert abs(dl - previous[0]) < 1e-8, "inconsistent DL values! (%g, %g) (!?)" % (dl, previous[0])
    else:
        # No previous result is available.
        bs_keys = [k for k in bs.keys() if type(k) != str]
        B_sup = max(max(bs_keys), B) if len(bs_keys) > 0 else B
        for Bi in bs_keys:
            if Bi > B and Bi < B_sup:
                B_sup = Bi
        if B_sup == B or not kwargs["shrink"]:
            # Start from scratch.
            raise RuntimeError("should not happen!")
        else:
            # Start from result with B_sup > B, and shrink it.
            if kwargs.get("verbose", False):
                print("(shrinking from B=%d to B=%d)" % (B_sup, B))
            state = bs[B_sup][1].copy()
        if __test__:
            assert state._BlockState__check_clabel(), "larger B clabel already invalidated!"

    # perform the actual minimization
    args = kwargs.copy()
    args["minimize_state"] = minimize_state
    if sparse_heuristic:
        args["dense"] = False
        args["multigraph"] = False
    #args["verbose"] = False

    state = multilevel_minimize(state, B, checkpoint=checkpoint, **args)

    if __test__:
        assert state._BlockState__check_clabel(), "clabel invalidated after minimize"
        assert state.B == B

    dl = get_b_dl(state, kwargs.get("dense", False),
                  kwargs.get("multigraph", False),
                  kwargs.get("nested_dl", False),
                  kwargs.get("complete", False),
                  kwargs.get("nested_overlap", False),
                  kwargs.get("dl_ent", False))


    if __test__:
        assert state._BlockState__check_clabel(), "clabel invalidated after minimize (?!)"

    if previous is None or dl < previous[0]:
        # the current result improved the previous one
        bs[B] = [dl, state]
        if kwargs.get("verbose", False):
            print("(using new result for B=%d with L=%g)" % (B, dl))

    else:
        # the previous result is better than the current one
        if kwargs.get("verbose", False):
            print("(kept old result for B=%d with L=%g [vs L=%g])" % (B, previous[0], dl))
        dl = previous[0]

        if __test__:
            tdl = get_b_dl(previous[1], kwargs.get("dense", False),
                           kwargs.get("multigraph", False),
                           kwargs.get("nested_dl", False),
                           kwargs.get("complete", False),
                           kwargs.get("nested_overlap", False),
                           kwargs.get("dl_ent", False))
            assert abs(dl - tdl) < 1e-8, "inconsistent DL values! (%g, %g)" % (dl, tdl)


    checkpoint_state[B]["done"] = True
    if __test__:
        assert not isinf(dl)
    return dl

def fibo(n):
    phi = (1 + sqrt(5)) / 2
    return int(round(phi ** n / sqrt(5)))

def fibo_n_floor(x):
    phi = (1 + sqrt(5)) / 2
    n = floor(log(x * sqrt(5) + 0.5) / log(phi))
    return int(n)

def get_mid(a, b):
    n = fibo_n_floor(b - a)
    return b - fibo(n - 1)

def is_fibo(x):
    return fibo(fibo_n_floor(x)) == x


def minimize_blockmodel_dl(g, deg_corr=True, overlap=False,
                           nonoverlap_init=False, dl=True, multigraph=True,
                           dense=False, sparse_heuristic=False, eweight=None,
                           vweight=None, clabel=None, c=0, nsweeps=100,
                           adaptive_sweeps=True, epsilon=1e-3, anneal=(1., 1.),
                           greedy_cooling=True, sequential=True, parallel=False,
                           r=2, nmerge_sweeps=10, max_B=None, min_B=None,
                           mid_B=None, checkpoint=None, minimize_state=None,
                           exhaustive=False, init_states=None, max_BE=None,
                           verbose=False, **kwargs):
    r"""Find the block partition of an unspecified size which minimizes the description
    length of the network, according to the stochastic blockmodel ensemble which
    best describes it.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be modelled.
    deg_corr : ``bool`` (optional, default: ``True``)
        If ``True``, the degree-corrected version of the blockmodel ensemble will
        be assumed, otherwise the traditional variant will be used.
    overlap : ``bool`` (optional, default: ``False``)
        If ``True``, the mixed-membership version of the blockmodel will be used.
    nonoverlap_init : ``bool`` (optional, default: ``False``)
        If ``True``, and `´overlap==True``, the minimization starts by first
        fitting the non-overlapping model, and using that as a starting state.
    dl : ``bool`` (optional, default: ``True``)
        If ``True``, the change in the whole description length will be
        considered after each vertex move, not only the entropy.
    multigraph : ``bool`` (optional, default: ``False``)
        If ``True``, the multigraph entropy will be used.
    dense : ``bool`` (optional, default: ``False``)
        If ``True``, the "dense" variant of the entropy will be computed.
    sparse_heuristic : ``bool`` (optional, default: ``False``)
        If ``True``, the sparse entropy will be used to find the best partition,
        but the dense entropy will be used to compare different partitions. This
        has an effect only if ``dense == True``.
    eweight : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Edge multiplicities (for multigraphs or block graphs).
    vweight : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Vertex multiplicities (for block graphs).
    clabel : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Constraint labels on the vertices. If supplied, vertices with different
        label values will not be clustered in the same group.
    c : ``float`` (optional, default: ``1.0``)
        This parameter specifies how often fully random moves are attempted,
        instead of more likely moves based on the inferred block partition.
        For ``c == 0``, no fully random moves are attempted, and for ``c == inf``
        they are always attempted.
    nsweeps : ``int`` (optional, default: ``10``)
        The number of sweeps done after each merge step to reach the local
        minimum.
    adaptive_sweeps : ``bool`` (optional, default: ``True``)
        If ``True``, the number of sweeps necessary for the local minimum will
        be estimated to be enough so that no more than ``epsilon * N`` nodes
        changes their states in the last ``nsweeps`` sweeps.
    epsilon : ``float`` (optional, default: ``1e-3``)
        Converge criterion for ``adaptive_sweeps``.
    anneal : pair of ``floats`` (optional, default: ``(1., 1.)``)
        The first value specifies the starting value for  ``beta`` of the MCMC
        steps, and the second value is the factor which is multiplied to ``beta``
        after each estimated equilibration (according to ``nsweeps`` and
        ``adaptive_sweeps``).
    greedy_cooling : ``bool`` (optional, default: ``True``)
        If ``True``, the value of ``beta`` of the MCMC steps are kept at
        infinity for all steps. Otherwise they change according to the ``anneal``
        parameter.
    sequential : ``bool`` (optional, default: ``True``)
        If ``True``, the move attempts on the vertices are done in sequential
        random order. Otherwise a total of `N` moves attempts are made, where
        `N` is the number of vertices, where each vertex can be selected with
        equal probability.
    parallel : ``bool`` (optional, default: ``False``)
        If ``True``, the updates are performed in parallel (multiple
        threads).
    r : ``float`` (optional, default: ``2.``)
        Agglomeration ratio for the merging steps. Each merge step will attempt
        to find the best partition into :math:`B_{i-1} / r` blocks, where
        :math:`B_{i-1}` is the number of blocks in the last step.
    nmerge_sweeps : `int` (optional, default: `10`)
        The number of merge sweeps done, where in each sweep a better merge
        candidate is searched for every block.
    max_B : ``int`` (optional, default: ``None``)
        Maximum number of blocks tried. If not supplied, it will be
        automatically determined.
    min_B : ``int`` (optional, default: `1`)
        Minimum number of blocks tried.
    mid_B : ``int`` (optional, default: ``None``)
        Middle of the range which brackets the minimum. If not supplied, will be
        automatically determined.
    checkpoint : function (optional, default: ``None``)
        If provided, this function will be called after each call to
        :func:`mcmc_sweep`. This can be used to store the current state, so it
        can be continued later. The function must have the following signature:

        .. code-block:: python

            def checkpoint(state, L, delta, nmoves, minimize_state):
                ...

        where `state` is either a :class:`~graph_tool.community.BlockState`
        instance or ``None``, `L` is the current description length, `delta` is
        the entropy difference in the last MCMC sweep, and `nmoves` is the
        number of accepted block membership moves. The ``minimize_state``
        argument is a :class:`~graph_tool.community.MinimizeState` instance
        which specifies the current state of the algorithm, which can be stored
        via :mod:`pickle`, and supplied via the ``minimize_state`` option below
        to continue from an interrupted run.

        This function will also be called when the MCMC has finished for the
        current value of :math:`B`, in which case ``state == None``, and the
        remaining parameters will be zero, except the last.
    minimize_state : :class:`~graph_tool.community.MinimizeState` (optional, default: ``None``)
        If provided, this will specify an exact point of execution from which
        the algorithm will continue. The expected object is a
        :class:`~graph_tool.community.MinimizeState`
        instance which will be passed to the callback of the ``checkpoint``
        option above, and  can be stored by :mod:`pickle`.
    exhaustive : ``bool`` (optional, default: ``False``)
        If ``True``, the best value of ``B`` will be found by testing all possible
        values, instead of performing a bisection search.
    init_states : ``list`` of :class:`~graph_tool.community.BlockState` or :class:`~graph_tool.community.OverlapBlockState` (optional, default: ``None``)
        If provided, this list of block states will be used when performing the
        minimization.
    max_BE : ``int`` (optional, default: ``1000``)
        If the number of blocks exceeds this number, a sparse representation of
        the block graph is used, which is slightly less efficient, but uses less
        memory,
    verbose : ``bool`` (optional, default: ``False``)
        If ``True``, verbose information is displayed.

    Returns
    -------
    state : :class:`~graph_tool.community.BlockState` or :class:`~graph_tool.community.OverlapBlockState`
       The block state containing the best model fit.

    Notes
    -----

    This algorithm attempts to find a block partition of an unspecified size
    which minimizes the description length of the network,

    .. math::

       \Sigma_{t/c} = \mathcal{S}_{t/c} + \mathcal{L}_{t/c},

    where :math:`\mathcal{S}_{t/c}` is the blockmodel entropy (as described in
    the docstring of :func:`mcmc_sweep`, :meth:`BlockState.entropy`
    :meth:`OverlapBlockState.entropy`) and :math:`\mathcal{L}_{t/c}` is the
    information necessary to describe the model (as described in the docstring
    of :func:`model_entropy`, :meth:`BlockState.entropy` and
    :meth:`OverlapBlockState.entropy`).

    The algorithm works by minimizing the entropy :math:`\mathcal{S}_{t/c}` for
    specific values of :math:`B` via :func:`mcmc_sweep` (with :math:`\beta = 1`
    and :math:`\beta\to\infty`), and minimizing :math:`\Sigma_{t/c}` via an
    one-dimensional Fibonacci search on :math:`B`. See
    [peixoto-parsimonious-2013]_ and [peixoto-efficient-2014]_ for more details.

    This algorithm has a complexity of :math:`O(\tau N\ln^2 B_{\text{max}})`,
    where :math:`N` is the number of nodes in the network, :math:`\tau` is the
    mixing time of the MCMC, and :math:`B_{\text{max}}` is the maximum number of
    blocks considered. If :math:`B_{\text{max}}` is not supplied, it is computed
    as :math:`\sim\sqrt{E}` via :func:`get_max_B`, in which case the complexity
    becomes :math:`O(\tau E\ln E)`.


    Examples
    --------
    .. testsetup:: mdl

       gt.seed_rng(42)
       np.random.seed(42)

    .. doctest:: mdl

       >>> g = gt.collection.data["polbooks"]
       >>> state = gt.minimize_blockmodel_dl(g)
       >>> b = state.b
       >>> gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=b, vertex_shape=b, output="polbooks_blocks_mdl.pdf")
       <...>

    .. testcleanup:: mdl

       gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=b, vertex_shape=b, output="polbooks_blocks_mdl.png")

    .. figure:: polbooks_blocks_mdl.*
       :align: center

       Block partition of a political books network, which minimizes the description
       length of the network according to the degree-corrected stochastic blockmodel.


    .. testsetup:: mdl_overlap

       gt.seed_rng(42)
       np.random.seed(42)

    .. doctest:: mdl_overlap

       >>> g = gt.collection.data["polbooks"]
       >>> state = gt.minimize_blockmodel_dl(g, overlap=True)
       >>> bv, *rest, bc = state.get_overlap_blocks()
       >>> eg = gt.get_block_edge_gradient(g, state.get_edge_blocks())
       >>> gt.graph_draw(g, g.vp["pos"], vertex_pie_fractions=bc,
       ...               vertex_pie_colors=bv, vertex_shape="pie",
       ...               edge_gradient=eg,
       ...               output="polbooks_overlap_blocks_mdl.pdf")
       <...>

    .. testcleanup:: mdl_overlap

       gt.graph_draw(g, g.vp["pos"], vertex_pie_fractions=bc,
                     vertex_pie_colors=bv, vertex_shape="pie",
                     edge_gradient=eg,
                     output="polbooks_overlap_blocks_mdl.png")

    .. figure:: polbooks_overlap_blocks_mdl.*
       :align: center

       Overlapping partition of a political books network, which minimizes the
       description length of the network according to the overlapping
       degree-corrected stochastic blockmodel.

    References
    ----------

    .. [holland-stochastic-1983] Paul W. Holland, Kathryn Blackmond Laskey,
       Samuel Leinhardt, "Stochastic blockmodels: First steps",
       Carnegie-Mellon University, Pittsburgh, PA 15213, U.S.A., :doi:`10.1016/0378-8733(83)90021-7`
    .. [faust-blockmodels-1992] Katherine Faust, and Stanley
       Wasserman. "Blockmodels: Interpretation and Evaluation." Social Networks
       14, no. 1-2 (1992): 5-61. :doi:`10.1016/0378-8733(92)90013-W`
    .. [karrer-stochastic-2011] Brian Karrer, and M. E. J. Newman. "Stochastic
       Blockmodels and Community Structure in Networks." Physical Review E 83,
       no. 1 (2011): 016107. :doi:`10.1103/PhysRevE.83.016107`.
    .. [peixoto-entropy-2012] Tiago P. Peixoto "Entropy of Stochastic Blockmodel
       Ensembles." Physical Review E 85, no. 5 (2012): 056122. :doi:`10.1103/PhysRevE.85.056122`,
       :arxiv:`1112.6028`.
    .. [peixoto-parsimonious-2013] Tiago P. Peixoto, "Parsimonious module inference in large networks",
       Phys. Rev. Lett. 110, 148701 (2013), :doi:`10.1103/PhysRevLett.110.148701`, :arxiv:`1212.4794`.
    .. [peixoto-efficient-2014] Tiago P. Peixoto, "Efficient Monte Carlo and greedy
       heuristic for the inference of stochastic block models", Phys. Rev. E 89, 012804 (2014),
       :doi:`10.1103/PhysRevE.89.012804`, :arxiv:`1310.4378`.
    .. [peixoto-model-2014] Tiago P. Peixoto, "Model selection and hypothesis
       testing for large-scale network models with overlapping groups",
       :arxiv:`1409.3059`.
    """

    nested_dl = kwargs.get("nested_dl", False)
    nested_overlap = kwargs.get("nested_overlap", False)
    nonoverlap_compare = kwargs.get("nonoverlap_compare", False)
    dl_ent = kwargs.get("dl_ent", False)

    if minimize_state is None:
        minimize_state = MinimizeState()

    if overlap and nonoverlap_init and minimize_state.init:
        if verbose:
            print("Non-overlapping initialization...")
        state = minimize_blockmodel_dl(g=g, eweight=eweight, vweight=vweight,
                                       deg_corr=deg_corr, dl=dl, dense=dense,
                                       multigraph=multigraph,
                                       sparse_heuristic=sparse_heuristic, c=c,
                                       nsweeps=nsweeps,
                                       adaptive_sweeps=adaptive_sweeps,
                                       epsilon=epsilon, anneal=anneal,
                                       greedy_cooling=greedy_cooling,
                                       sequential=sequential, parallel=parallel,
                                       r=r, nmerge_sweeps=nmerge_sweeps,
                                       max_B=max_B, min_B=min_B, mid_B=mid_B,
                                       clabel=clabel if isinstance(clabel, PropertyMap) else None,
                                       checkpoint=checkpoint,
                                       minimize_state=minimize_state,
                                       exhaustive=exhaustive, max_BE=max_BE,
                                       nested_dl=nested_dl, overlap=False,
                                       init_states=None, dl_ent=dl_ent,
                                       verbose=verbose)
        state = state.copy(overlap=True, clabel=clabel)
        unilevel_minimize(state, nsweeps=nsweeps, epsilon=epsilon, c=c, dl=dl,
                          nmerge_sweeps=nmerge_sweeps, sequential=sequential)
        max_B = state.B
        init_states = [state]

        minimize_state.clear()
        minimize_state.init = False

        if min_B is None:
            min_B = state.clabel.a.max() + 1

        if verbose:
            print("Overlapping minimization starting from B=", max_B)

    if min_B is None:
        if clabel is None:
            min_B = 1
        elif isinstance(clabel, PropertyMap):
            min_B = clabel.fa.max() + 1
        else:
            min_B = clabel.max() + 1
    elif clabel is not None:
        C = clabel.fa.max() + 1 if isinstance(clabel, PropertyMap) else clabel.max() + 1
        if C > min_B:
            raise ValueError("value of min_B=%d is not consistent with the enforced constraints of size %d" % (min_B, C))

    if max_B is None:
        if dense:
            max_B = max(g.num_vertices(), 1)
        else:
            max_B = get_max_B(g.num_vertices(), g.num_edges(), g.is_directed())
        if verbose:
            print("max_B:", max_B)

    if mid_B is None:
        mid_B = get_mid(min_B, max_B)

    B_lims = (min_B, max_B)

    greedy = greedy_cooling
    shrink = True

    b_cache = minimize_state.b_cache
    checkpoint_state = minimize_state.checkpoint_state

    kwargs = dict(nsweeps=nsweeps, adaptive_sweeps=adaptive_sweeps, c=c,
                  sequential=sequential, parallel=parallel, shrink=shrink, r=r,
                  anneal=anneal, greedy=greedy, epsilon=epsilon,
                  nmerge_sweeps=nmerge_sweeps, deg_corr=deg_corr, dense=dense,
                  multigraph=multigraph, dl=dl,
                  sparse_heuristic=sparse_heuristic, checkpoint=checkpoint,
                  minimize_state=minimize_state, nested_dl=nested_dl,
                  nested_overlap=nested_overlap,
                  nonoverlap_compare=nonoverlap_compare, dl_ent=dl_ent,
                  verbose=verbose)

    if init_states is not None:
        for state in init_states:
            dl = get_b_dl(state, kwargs.get("dense", False),
                          kwargs.get("multigraph", False),
                          kwargs.get("nested_dl", False),
                          kwargs.get("complete", False),
                          kwargs.get("nested_overlap", False),
                          kwargs.get("dl_ent", False))
            b_cache[state.B] = [dl, state]

    B_init = True
    for Bi, bstate in b_cache.items():
        if Bi >= max_B:
            B_init = False

    if B_init:
        if overlap:
            state = OverlapBlockState(g, B=2 * g.num_edges(), deg_corr=deg_corr,
                                      vweight=vweight, eweight=eweight,
                                      clabel=clabel, max_BE=max_BE)
        else:
            state = BlockState(g, B=g.num_vertices(), deg_corr=deg_corr,
                               vweight=vweight, eweight=eweight, clabel=clabel,
                               max_BE=max_BE)

            if __test__:
                assert state._BlockState__check_clabel(), "clabel invalid at copying!"

        if __test__:
            assert state._BlockState__check_clabel(), "clabel invalid at creation!"

        dl = get_b_dl(state, kwargs.get("dense", False),
                      kwargs.get("multigraph", False),
                      kwargs.get("nested_dl", False),
                      kwargs.get("complete", False),
                      kwargs.get("nested_overlap", False),
                      kwargs.get("dl_ent", False))
        b_cache[state.B] = [dl, state]

    if exhaustive:
        if max_B not in b_cache:
            Bi = max(b_cache.keys())
            state = b_cache[Bi][1]
            state = multilevel_minimize(state, B=max_B, **kwargs)

        for B in reversed(range(min_B, max_B + 1)):
            if B in b_cache:
                state = b_cache[B][1]
                if checkpoint_state[B].get("done", False):
                    continue

            args = kwargs.copy()
            if sparse_heuristic:
                args["dense"] = False

            state = multilevel_minimize(state, B, **args)

            dl = get_b_dl(state, kwargs.get("dense", False),
                          kwargs.get("multigraph", False),
                          kwargs.get("nested_dl", False),
                          kwargs.get("complete", False),
                          kwargs.get("nested_overlap", False),
                          kwargs.get("dl_ent", False))

            b_cache[B] = [dl, state]

            if verbose:
                print("Result for B=%d: L=%g" % (B, dl))

        min_dl = float(inf)
        best_B = None
        for Bi in b_cache.keys():
            if b_cache[Bi][0] <= min_dl:
                min_dl = b_cache[Bi][0]
                best_B = Bi
        if verbose:
            print("Best result: B=%d, L=%g" % (best_B, min_dl))

        return b_cache[best_B][1]


    def cleanup_cache(b_cache, B_min, B_max):
        best_B = None
        min_dl = float("inf")
        for Bi in b_cache.keys():
            if b_cache[Bi][0] <= min_dl:
                min_dl = b_cache[Bi][0]
                best_B = Bi

        del_Bs = []

        for Bi in b_cache.keys():
            if (Bi < B_min or Bi > B_max) and Bi != best_B:
                del_Bs.append(Bi)

        for Bi in del_Bs:
            del b_cache[Bi]

    # Initial bracketing
    while True:
        f_max = get_state_dl(B=max_B, **kwargs)
        f_mid = get_state_dl(B=mid_B, **kwargs)
        f_min = get_state_dl(B=min_B, **kwargs)

        if verbose:
            print("Current bracket:", (min_B, mid_B, max_B), (f_min, f_mid, f_max))

        if checkpoint is not None:
             checkpoint(None, 0, 0, 0, minimize_state)

        cleanup_cache(b_cache, min_B, max_B)

        if f_max > f_mid > f_min:
            max_B = mid_B
            mid_B = get_mid(min_B, mid_B)
        elif f_max < f_mid < f_min:
            min_B = mid_B
            mid_B = get_mid(mid_B, max_B)
        else:
            break

    # Fibonacci search
    while True:
        if max_B - mid_B > mid_B - min_B:
            x = get_mid(mid_B, max_B)
        else:
            x = get_mid(min_B, mid_B)

        f_x = get_state_dl(B=x, **kwargs)
        f_mid = get_state_dl(B=mid_B, **kwargs)

        if verbose:
            print("Current bracket:",
                  (min_B, mid_B, max_B), (get_state_dl(B=min_B, **kwargs), f_mid,
                                          get_state_dl(B=max_B, **kwargs)))
            print("Bisect at", x, "with L=%g" % f_x)

        if max_B - mid_B <= 1:
            min_dl = float(inf)
            best_B = None
            for Bi in b_cache.keys():
                if Bi < B_lims[0] or Bi > B_lims[1]:
                    continue
                if b_cache[Bi][0] <= min_dl:
                    min_dl = b_cache[Bi][0]
                    best_B = Bi
            if verbose:
                print("Best result: B=%d, L=%g" % (best_B, min_dl))

            return b_cache[best_B][1]

        if checkpoint is not None:
            checkpoint(None, 0, 0, 0, minimize_state)

        if f_x < f_mid:
            if max_B - mid_B > mid_B - min_B:
                min_B = mid_B
                mid_B = x
            else:
                max_B = mid_B
                mid_B = x
        else:
            if max_B - mid_B > mid_B - min_B:
                max_B = x
            else:
                min_B = x

        cleanup_cache(b_cache, min_B, max_B)



def collect_edge_marginals(state, p=None):
    r"""Collect the edge marginal histogram, which counts the number of times
    the endpoints of each node have been assigned to a given block pair.

    This should be called multiple times, after repeated runs of the
    :func:`mcmc_sweep` function.

    Parameters
    ----------
    state : :class:`~graph_tool.community.BlockState`
        The block state.
    p : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Edge property map with vector-type values, storing the previous block
        membership counts.  Each vector entry corresponds to ``b[i] + B *
        b[j]``, where ``b`` is the block membership and ``i = min(source(e),
        target(e))`` and ``j = max(source(e), target(e))``. If not provided, an
        empty histogram will be created.

    Returns
    -------
    p : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Vertex property map with vector-type values, storing the accumulated
        block membership counts.


    Examples
    --------
    .. testsetup:: collect_edge_marginals

       gt.seed_rng(42)
       np.random.seed(42)

    .. doctest:: collect_edge_marginals

       >>> g = gt.collection.data["polbooks"]
       >>> state = gt.BlockState(g, B=4, deg_corr=True)
       >>> pe = None
       >>> for i in range(1000):        # remove part of the transient
       ...     ds, nmoves = gt.mcmc_sweep(state)
       >>> for i in range(1000):
       ...     ds, nmoves = gt.mcmc_sweep(state)
       ...     pe = gt.collect_edge_marginals(state, pe)
       >>> gt.bethe_entropy(state, pe)[0]
       17.609773262509986
    """

    if p is None:
        p = state.g.new_edge_property("vector<int>")

    libcommunity.edge_marginals(state.g._Graph__graph,
                                state.bg._Graph__graph,
                                state.B,
                                _prop("v", state.g, state.b),
                                _prop("e", state.g, p))
    return p

def collect_vertex_marginals(state, p=None):
    r"""Collect the vertex marginal histogram, which counts the number of times a
    node was assigned to a given block.

    This should be called multiple times, after repeated runs of the
    :func:`mcmc_sweep` function.

    Parameters
    ----------
    state : :class:`~graph_tool.community.BlockState`
        The block state.
    p : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Vertex property map with vector-type values, storing the previous block
        membership counts. If not provided, an empty histogram will be created.

    Returns
    -------
    p : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Vertex property map with vector-type values, storing the accumulated
        block membership counts.

    Examples
    --------
    .. testsetup:: cvm

       gt.seed_rng(42)
       np.random.seed(42)

    .. doctest:: cvm

       >>> g = gt.collection.data["polbooks"]
       >>> state = gt.BlockState(g, B=4, deg_corr=True)
       >>> pv = None
       >>> for i in range(1000):        # remove part of the transient
       ...     ds, nmoves = gt.mcmc_sweep(state)
       >>> for i in range(1000):
       ...     ds, nmoves = gt.mcmc_sweep(state)
       ...     pv = gt.collect_vertex_marginals(state, pv)
       >>> gt.mf_entropy(state, pv)
       20.117550557730116
       >>> gt.graph_draw(g, pos=g.vp["pos"], vertex_shape="pie", vertex_pie_fractions=pv, output="polbooks_blocks_soft_B4.pdf")
       <...>

    .. testcleanup:: cvm

       gt.graph_draw(g, pos=g.vp["pos"], vertex_shape="pie", vertex_pie_fractions=pv, output="polbooks_blocks_soft_B4.png")

    .. figure:: polbooks_blocks_soft_B4.*
       :align: center

       "Soft" block partition of a political books network with :math:`B=4`.

    """
    B = state.B

    if p is None:
        p = state.g.new_vertex_property("vector<int>")

    libcommunity.vertex_marginals(state.g._Graph__graph,
                                  _prop("v", state.g, state.b),
                                  _prop("v", state.g, p))
    return p

def bethe_entropy(state, p):
    r"""Compute the Bethe entropy given the edge block membership marginals.

    Parameters
    ----------
    state : :class:`~graph_tool.community.BlockState`
        The block state.
    p : :class:`~graph_tool.PropertyMap`
        Edge property map with vector-type values, storing the previous block
        membership counts.  Each vector entry corresponds to ``b[i] + B *
        b[j]``, where ``b`` is the block membership and ``i = min(source(e),
        target(e))`` and ``j = max(source(e), target(e))``.

    Returns
    -------
    H : ``float``
        The Bethe entropy value (in `nats <http://en.wikipedia.org/wiki/Nat_%28information%29>`_)
    Hmf : ``float``
        The "mean field" entropy value (in `nats <http://en.wikipedia.org/wiki/Nat_%28information%29>`_),
        as would be returned by the :func:`mf_entropy` function.
    pv : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Vertex property map with vector-type values, storing the accumulated
        block membership counts. These are the node marginals, as would be
        returned by the :func:`collect_vertex_marginals` function.

    Notes
    -----

    The Bethe entropy is defined as,

    .. math::

        H = -\sum_{e,(r,s)}\pi_{(r,s)}^e\ln\pi_{(r,s)}^e - \sum_{v,r}(1-k_i)\pi_r^v\ln\pi_r^v,

    where :math:`\pi_{(r,s)}^e` is the marginal probability that the endpoints
    of the edge :math:`e` belong to blocks :math:`(r,s)`, and :math:`\pi_r^v` is
    the marginal probability that vertex :math:`v` belongs to block :math:`r`,
    and :math:`k_i` is the degree of vertex :math:`v` (or total degree for
    directed graphs).

    References
    ----------
    .. [mezard-information-2009] Marc Mézard, Andrea Montanari, "Information,
       Physics, and Computation", Oxford Univ Press, 2009.
    """
    B = state.B
    H = 0
    pv =  state.g.new_vertex_property("vector<double>")

    H, sH, Hmf, sHmf  = libcommunity.bethe_entropy(state.g._Graph__graph,
                                                   state.B,
                                                   _prop("e", state.g, p),
                                                   _prop("v", state.g, pv))
    return H, Hmf, pv


def mf_entropy(state, p):
    r"""Compute the "mean field" entropy given the vertex block membership marginals.

    Parameters
    ----------
    state : :class:`~graph_tool.community.BlockState`
        The block state.
    p : :class:`~graph_tool.PropertyMap`
        Vertex property map with vector-type values, storing the accumulated block
        membership counts.

    Returns
    -------
    Hmf : ``float``
        The "mean field" entropy value (in `nats <http://en.wikipedia.org/wiki/Nat_%28information%29>`_).

    Notes
    -----

    The "mean field" entropy is defined as,

    .. math::

        H = - \sum_{v,r}\pi_r^v\ln\pi_r^v,

    where :math:`\pi_r^v` is the marginal probability that vertex :math:`v`
    belongs to block :math:`r`.

    References
    ----------
    .. [mezard-information-2009] Marc Mézard, Andrea Montanari, "Information,
       Physics, and Computation", Oxford Univ Press, 2009.

    """
    H = 0
    for v in state.g.vertices():
        N = p[v].a.sum()
        if N == 0:
            continue
        pvi = asarray(p[v].a, dtype="float") /  N
        pvi = pvi[pvi > 0]
        H -= (pvi * log(pvi)).sum()
    return H


def condensation_graph(g, prop, vweight=None, eweight=None, avprops=None,
                       aeprops=None, self_loops=False):
    r"""
    Obtain the condensation graph, where each vertex with the same 'prop' value is condensed in one vertex.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be modelled.
    prop : :class:`~graph_tool.PropertyMap`
        Vertex property map with the community partition.
    vweight : :class:`~graph_tool.PropertyMap` (optional, default: None)
        Vertex property map with the optional vertex weights.
    eweight : :class:`~graph_tool.PropertyMap` (optional, default: None)
        Edge property map with the optional edge weights.
    avprops : list of :class:`~graph_tool.PropertyMap` (optional, default: None)
        If provided, the average value of each property map in this list for
        each vertex in the condensed graph will be computed and returned.
    aeprops : list of :class:`~graph_tool.PropertyMap` (optional, default: None)
        If provided, the average value of each property map in this list for
        each edge in the condensed graph will be computed and returned.
    self_loops : ``bool`` (optional, default: ``False``)
        If ``True``, self-loops due to intra-block edges are also included in
        the condensation graph.

    Returns
    -------
    condensation_graph : :class:`~graph_tool.Graph`
        The community network
    prop : :class:`~graph_tool.PropertyMap`
        The community values.
    vcount : :class:`~graph_tool.PropertyMap`
        A vertex property map with the vertex count for each community.
    ecount : :class:`~graph_tool.PropertyMap`
        An edge property map with the inter-community edge count for each edge.
    va : list of :class:`~graph_tool.PropertyMap`
        A list of vertex property maps with average values of the properties
        passed via the ``avprops`` parameter.
    ea : list of :class:`~graph_tool.PropertyMap`
        A list of edge property maps with average values of the properties
        passed via the ``avprops`` parameter.

    See Also
    --------
    community_structure: Obtain the community structure
    modularity: Calculate the network modularity
    condensation_graph: Network of communities, or blocks

    Notes
    -----
    Each vertex in the condensation graph represents one community in the
    original graph (vertices with the same 'prop' value), and the edges
    represent existent edges between vertices of the respective communities in
    the original graph.

    Examples
    --------

    .. testsetup:: condensation_graph

       gt.seed_rng(43)
       np.random.seed(42)

    Let's first obtain the best block partition with ``B=5``.

    .. doctest:: condensation_graph

       >>> g = gt.collection.data["polbooks"]
       >>> state = gt.BlockState(g, B=5, deg_corr=True)
       >>> for i in range(1000):        # remove part of the transient
       ...     ds, nmoves = gt.mcmc_sweep(state)
       >>> for i in range(1000):
       ...     ds, nmoves = gt.mcmc_sweep(state, beta=float("inf"))
       >>> b = state.get_blocks()
       >>> gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=b, vertex_shape=b, output="polbooks_blocks_B5.pdf")
       <...>

    Now we get the condensation graph:

    .. doctest:: condensation_graph

       >>> bg, bb, vcount, ecount, avp, aep = gt.condensation_graph(g, b, avprops=[g.vp["pos"]], self_loops=True)
       >>> gt.graph_draw(bg, pos=avp[0], vertex_fill_color=bb, vertex_shape=bb,
       ...               vertex_size=gt.prop_to_size(vcount, mi=40, ma=100),
       ...               edge_pen_width=gt.prop_to_size(ecount, mi=2, ma=10),
       ...               output="polbooks_blocks_B5_cond.pdf")
       <...>

    .. testcleanup:: condensation_graph

       gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=b, vertex_shape=b, output="polbooks_blocks_B5.png")
       gt.graph_draw(bg, pos=avp[0], vertex_fill_color=bb, vertex_shape=bb,
                     vertex_size=gt.prop_to_size(vcount, mi=40, ma=100),
                     edge_pen_width=gt.prop_to_size(ecount, mi=2, ma=10),
                     output="polbooks_blocks_B5_cond.png")

    .. figure:: polbooks_blocks_B5.*
       :align: center

       Block partition of a political books network with :math:`B=5`.

    .. figure:: polbooks_blocks_B5_cond.*
       :align: center

       Condensation graph of the obtained block partition.

    """
    gp = Graph(directed=g.is_directed())
    if vweight is None:
        vcount = gp.new_vertex_property("int32_t")
    else:
        vcount = gp.new_vertex_property(vweight.value_type())
    if eweight is None:
        ecount = gp.new_edge_property("int32_t")
    else:
        ecount = gp.new_edge_property(eweight.value_type())

    if prop is g.vertex_index:
        prop = prop.copy(value_type="int32_t")
    cprop = gp.new_vertex_property(prop.value_type())

    if avprops is None:
        avprops = []
    avp = []
    r_avp = []
    for p in avprops:
        if p is g.vertex_index:
            p = p.copy(value_type="int")
        if "string" in p.value_type():
            raise ValueError("Cannot compute average of string properties!")
        temp = g.new_vertex_property(p.value_type())
        cp = gp.new_vertex_property(p.value_type())
        avp.append((_prop("v", g, p), _prop("v", g, temp), _prop("v", g, cp)))
        r_avp.append(cp)

    if aeprops is None:
        aeprops = []
    aep = []
    r_aep = []
    for p in aeprops:
        if p is g.edge_index:
            p = p.copy(value_type="int")
        if "string" in p.value_type():
            raise ValueError("Cannot compute average of string properties!")
        temp = g.new_edge_property(p.value_type())
        cp = gp.new_edge_property(p.value_type())
        aep.append((_prop("e", g, p), _prop("e", g, temp), _prop("e", g, cp)))
        r_aep.append(cp)

    libcommunity.community_network(g._Graph__graph,
                                   gp._Graph__graph,
                                   _prop("v", g, prop),
                                   _prop("v", gp, cprop),
                                   _prop("v", gp, vcount),
                                   _prop("e", gp, ecount),
                                   _prop("v", g, vweight),
                                   _prop("e", g, eweight),
                                   self_loops)
    u = GraphView(g, directed=True, reversed=False)
    libcommunity.community_network_vavg(u._Graph__graph,
                                        gp._Graph__graph,
                                        _prop("v", g, prop),
                                        _prop("v", gp, cprop),
                                        _prop("v", gp, vcount),
                                        _prop("v", g, vweight),
                                        avp)
    u = GraphView(g, directed=True)
    libcommunity.community_network_eavg(u._Graph__graph,
                                        gp._Graph__graph,
                                        _prop("v", g, prop),
                                        _prop("v", gp, cprop),
                                        _prop("e", gp, ecount),
                                        _prop("e", g, eweight),
                                        aep,
                                        self_loops)
    return gp, cprop, vcount, ecount, r_avp, r_aep

from . overlap_blockmodel import *
