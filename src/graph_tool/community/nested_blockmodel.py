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
import random
from numpy import *
import numpy
from scipy.optimize import fsolve, fminbound
import scipy.special
from collections import defaultdict
import copy
import heapq

from . blockmodel import *
from . blockmodel import __test__

class NestedBlockState(object):
    r"""This class encapsulates the nested block state of a given graph.

    This must be instantiated and used by functions such as :func:`nested_mcmc_sweep`.

    The instances of this class contain a data member called ``levels``, which
    is a list of :class:`~graph_tool.community.BlockState` (or
    :class:`~graph_tool.community.OverlapBlockState`) instances, containing the
    entire nested hierarchy.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    eweight : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Edge weights (i.e. multiplicity).
    vweight : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Vertex weights (i.e. multiplicity).
    bs : list of :class:`~graph_tool.PropertyMap` or :class:`~numpy.ndarray` instances (optional, default: ``None``)
        Initial block labels on the vertices, for each hierarchy level.
    Bs : list of ``int`` (optional, default: ``None``)
        Number of blocks for each hierarchy level.
    deg_corr : ``bool`` (optional, default: ``True``)
        If ``True``, the degree-corrected version of the blockmodel ensemble will
        be used in the bottom level, otherwise the traditional variant will be used.
    overlap : ``bool`` (optional, default: ``False``)
        If ``True``, the mixed-membership version of the blockmodel will be used
        at the lowest level.
    clabel : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Constraint labels on the vertices. If supplied, vertices with different
        label values will not be clustered in the same group.
    max_BE : ``int`` (optional, default: ``1000``)
        If the number of blocks exceeds this number, a sparse representation of
        the block graph is used, which is slightly less efficient, but uses less
        memory,
    """

    def __init__(self, g, eweight=None, vweight=None, bs=None,
                 Bs=None, deg_corr=True, overlap=False, clabel=None,
                 max_BE=1000):
        L = len(Bs) if Bs is not None else len(bs)
        self.g = cg = g
        self.vweight = vcount = vweight
        self.eweight = ecount = eweight

        self.levels = []
        self.overlap = overlap
        self.deg_corr = deg_corr
        self.clabel = clabel if clabel is not None else g.new_vertex_property("int")

        for l in range(L):
            Bl = Bs[l] if Bs is not None else None
            bl = None
            if bs is not None:
                if isinstance(bs[l], PropertyMap):
                    bl = cg.own_property(bs[l])
                else:
                    bl = bs[l]

            if l == 0:
                if overlap:
                    state = OverlapBlockState(g, B=Bl, b=bl,
                                              eweight=ecount,
                                              vweight=vcount,
                                              deg_corr=deg_corr != False,
                                              clabel=self.clabel,
                                              max_BE=max_BE)
                    self.clabel = state.clabel.copy()
                    state.clabel.a = 0
                else:
                    state = BlockState(g, B=Bl, b=bl,
                                       eweight=ecount,
                                       vweight=vcount,
                                       deg_corr=deg_corr != False,
                                       #clabel=self.clabel,
                                       max_BE=max_BE)
            else:
                state = self.levels[-1].get_block_state(b=bl,
                                                        overlap=self.overlap == "full",
                                                        deg_corr=self.deg_corr == "full")[0]
                if __test__:
                    assert not state.deg_corr, "upper levels must be non-deg-corr"

            self.levels.append(state)

        # for l in range(len(self.levels) - 1):
        #     clabel = self.__project_partition(l, l+1)
        #     self.levels[l].clabel = clabel


    def copy(self, deg_corr=None, overlap=None, clabel=None):
        r"""Copies the block state. The parameters override the state properties, and
         have the same meaning as in the constructor.."""

        bs = [s.b.a for s in self.levels]
        if overlap is None:
            overlap = self.overlap
        elif self.overlap and not overlap:
            raise ValueError("Cannot automatically convert overlapping nested state to nonoverlapping")
        elif not self.overlap and overlap:
            s = self.level[0].copy(overlap=True)
            bs[0] = s.b.a
        if deg_corr is None:
            deg_corr = self.deg_corr
        return NestedBlockState(self.g, self.eweight, self.vweight,
                                bs, deg_corr=deg_corr, overlap=overlap,
                                clabel=clabel, max_BE=self.levels[0].max_BE)

    def __getstate__(self):
        state = dict(g=self.g,
                     eweight=self.eweight,
                     vweight=self.vweight,
                     overlap=self.overlap,
                     bs=[array(s.b.a) for s in self.levels],
                     clabel=self.clabel,
                     deg_corr=self.deg_corr,
                     max_BE=self.levels[0].max_BE)
        return state

    def __setstate__(self, state):
        self.__init__(**state)
        return state

    def __project_partition(self, l, j):
        """Project partition of level 'j' onto level 'l'"""
        if self.overlap != "full":
            b = self.levels[l].b.copy()
            for i in range(l + 1, j + 1):
                clabel = self.levels[i].b.copy()
                pmap(b, clabel)
        else:
            b = self.levels[j].b.copy()
        return b

    def __propagate_clabel(self, l):
        clabel = self.clabel.copy()
        for j in range(l):
            bclabel = self.levels[j].bg.new_vertex_property("int")
            reverse_map(self.levels[j].b, bclabel)
            pmap(bclabel, clabel)
            clabel = bclabel
        return clabel

    def __consistency_check(self, op, l):

        print("consistency check after", op, "at level", l)

        for j in range(len(self.levels)):

            c_state = self.levels[j].copy()
            S1 = self.levels[j].entropy()
            S2 = c_state.entropy()

            assert abs(S1 - S2) < 1e-8 and not isnan(S2) and not isnan(S1), "inconsistency at level %d after %s of level %d, different entropies of copies! (%g, %g)" % (j, op, l, S1, S2)


            if self.levels[j].wr.a.min() == 0:
                print("WARNING: empty blocks at level", j)
            if self.levels[j].b.a.max() + 1 != self.levels[j].B:
                print("WARNING: b.max() + 1 != B at level", j, self.levels[j].b.max() + 1, self.levels[j].B)

        for j in range(len(self.levels) - 1):

            B = self.levels[j].b.a.max() + 1
            bg_state = self.levels[j].get_block_state(b=self.levels[j+1].b.copy(),
                                                      overlap=self.overlap == "full",
                                                      deg_corr=self.deg_corr == "full")[0]

            S1 = bg_state.entropy(dense=True and self.deg_corr != "full", multigraph=True)
            S2 = self.levels[j+1].entropy(dense=True and self.deg_corr != "full", multigraph=True)

            if self.levels[j].B != self.levels[j+1].N or S1 != S2:
                self.print_summary()

                from graph_tool.topology import similarity
                print(bg_state)
                print(self.levels[j+1])
                print("N, B:", bg_state.N, bg_state.B)
                print("N, B:", self.levels[j + 1].N, self.levels[j + 1].B)
                print("similarity:", similarity(bg_state.g, self.levels[j+1].g))
                print("b:", bg_state.b.a)
                print("b:", self.levels[j+1].b.a)

                print("wr:", bg_state.wr.a)
                print("wr:", self.levels[j+1].wr.a)

                print("mrs:", bg_state.mrs.a)
                print("mrs:", self.levels[j+1].mrs.a)

                print("eweight:", bg_state.eweight.a)
                print("eweight:", self.levels[j+1].eweight.a)

                print("vweight:", bg_state.vweight.a)
                print("vweight:", self.levels[j+1].vweight.a)


            assert abs(S1 - S2) < 1e-6 and not isnan(S2) and not isnan(S1), "inconsistency at level %d after %s of level %d, different entropies (%g, %g)" % (j, op, l, S1, S2)
            assert self.levels[j].B == self.levels[j+1].N, "inconsistency  at level %d after %s of level %d, different sizes" % (j + 1, op, l)


            # verify hierarchy / clabel consistency
            clabel = self.__project_partition(0, j)
            self.levels[0].clabel.a = self.clabel.a
            assert self.levels[0]._BlockState__check_clabel(),  "inconsistency at level %d after %s of level %d, clabel invalidated" % (j + 1, op, l)
            self.levels[0].clabel.a = 0

            # verify hierarchy consistency
            clabel = self.__project_partition(j, j + 1)
            self.levels[0].clabel.a = self.clabel.a
            assert self.levels[0]._BlockState__check_clabel(),  "inconsistency at level %d after %s of level %d, partition not compatible with upper level" % (j + 1, op, l)
            self.levels[0].clabel.a = 0


    def __rebuild_level(self, l, b, clabel=None):
        r"""Replace level ``l`` given the new partition ``b``, and the
        projected upper level partition clabel."""

        if __test__:
            assert clabel is not None or l == len(self.levels) - 1, "clabel not given for intermediary level"

        if clabel is None:
            clabel = self.levels[l].g.new_vertex_property("int")

        old_b = b.copy()

        state = self.levels[l].copy(b=b, clabel=clabel.a)
        self.levels[l] = state

        # upper level
        bclabel = state.get_bclabel()
        bstate = self.levels[l].get_block_state(b=bclabel,
                                                overlap=self.overlap == "full",
                                                deg_corr=self.deg_corr == "full")[0]
        if l == len(self.levels) - 1:
            self.levels.append(None)
        self.levels[l + 1] = bstate

        self.levels[l].clabel.a = 0
        self.levels[l + 1].clabel.a = 0

        # if l + 1 < len(self.levels) - 1:
        #     self.levels[l + 1].clabel = self.__project_partition(l + 1, l + 2)


        if l + 1 < len(self.levels) - 1:
            bstate = self.levels[l + 1].get_block_state(b=self.levels[l + 2].b,
                                                        overlap=self.overlap == "full",
                                                        deg_corr=self.deg_corr == "full")[0]

            if __test__:
                from graph_tool.topology import similarity
                print("- similarity:", similarity(bstate.g, self.levels[l + 2].g))
                # if similarity(bstate.g, self.levels[l + 2].g) < 1:
                #     from graph_tool.draw import graph_draw
                #     graph_draw(bstate.g)
                #     graph_draw(self.levels[l + 2].g)


                if abs(bstate.entropy() - self.levels[l + 2].entropy()) > 1e-6:
                    print("********************** inconsistent rebuild! **************************")

                    print(bstate.b.a)
                    print(self.levels[l + 2].b.a)

                    print(bstate.wr.a)
                    print(self.levels[l + 2].wr.a)

                    print(bstate.eweight.a)
                    print(self.levels[l + 2].eweight.a)

                    # from graph_tool.draw import graph_draw, prop_to_size
                    # graph_draw(bstate.g, vertex_fill_color=bstate.b)
                    # graph_draw(self.levels[l + 2].g, vertex_fill_color=self.levels[l + 2].b)
                    nclabel = self.__project_partition(l, l + 1)
                    print(nclabel.a)
                    print(clabel.a)
                    print(self.levels[l].b.a)
                    print(self.levels[l+1].b.a)
                    print(self.levels[l+2].b.a)
                    print(bstate.b.a)
                print ("DS", l, l + 1, bstate.entropy(), self.levels[l + 2].entropy())

        B = self.levels[l].B

        if __test__:
            self.__consistency_check("rebuild", l)

    def __delete_level(self, l):
        if l == 0:
            raise ValueError("cannot delete level l=0")
        b = self.__project_partition(l - 1, l)
        if l < len(self.levels) - 1:
            clabel = self.__project_partition(l - 1, l + 1)
        else:
            clabel = None

        del self.levels[l]

        self.__rebuild_level(l - 1, b=b, clabel=clabel)

        if __test__:
            self.__consistency_check("delete", l)

    def __duplicate_level(self, l):
        assert l > 0, "attempted to duplicate level 0"
        if not self.levels[l].overlap:
            bstate = self.levels[l].copy(b=self.levels[l].g.vertex_index.copy("int"))
        else:
            bstate = self.levels[l].copy(b=arange(self.levels[l].g.num_vertices()))
        self.levels.insert(l, bstate)

        if __test__:
            self.__consistency_check("duplicate", l)


    def level_entropy(self, l, complete=True, dense=False, multigraph=True,
                      norm=True, dl_ent=False):
        r"""Compute the description length of hierarchy level l.

        Parameters
        ----------
        l : ``int``
            Hierarchy level.
        complete : ``bool`` (optional, default: ``False``)
            If ``True``, the complete entropy will be returned, including constant
            terms not relevant to the block partition.
        dense : ``bool`` (optional, default: ``False``)
            If ``True``, the "dense" variant of the entropy will be computed.
        multigraph : ``bool`` (optional, default: ``True``)
            If ``True``, the multigraph entropy will be used.
        norm : ``bool`` (optional, default: ``True``)
            If ``True``, the entropy will be "normalized" by dividing by the
            number of edges.
        dl_ent : ``bool`` (optional, default: ``False``)
            If ``True``, the description length of the degree sequence will be
            approximated by its entropy.
        """

        bstate = self.levels[l]

        S = bstate.entropy(dl=True, partition_dl=True,
                           dense=dense or (l > 0 and self.deg_corr != "full"),
                           multigraph=multigraph or l > 0,
                           complete=complete or (l > 0 and self.deg_corr == "full"),
                           norm=norm,
                           dl_ent=dl_ent)

        return S

    def entropy(self, complete=True, dense=False, multigraph=True, norm=True,
                dl_ent=False):
        r"""Compute the description length of the entire hierarchy.

        Parameters
        ----------
        complete : ``bool`` (optional, default: ``False``)
            If ``True``, the complete entropy will be returned, including constant
            terms not relevant to the block partition.
        dense : ``bool`` (optional, default: ``False``)
            If ``True``, the "dense" variant of the entropy will be computed.
        multigraph : ``bool`` (optional, default: ``True``)
            If ``True``, the multigraph entropy will be used.
        norm : ``bool`` (optional, default: ``True``)
            If ``True``, the entropy will be "normalized" by dividing by the
            number of edges.
        dl_ent : ``bool`` (optional, default: ``False``)
            If ``True``, the description length of the degree sequence will be
            approximated by its entropy.
        """

        S = 0
        for l in range(len(self.levels)):
            S += self.level_entropy(l, complete=complete, dense=dense,
                                    multigraph=multigraph, norm=norm,
                                    dl_ent=dl_ent)
        return S

    def get_bstack(self):
        r"""Return the nested levels as individual graphs.

        This returns a list of :class:`~graph_tool.Graph` instances
        representing the inferred hierarchy at each level. Each graph has two
        internal vertex and edge property maps named "count" which correspond to
        the vertex and edge counts at the lower level, respectively. Additionally,
        an internal vertex property map named "b" specifies the block partition.
        """

        bstack = []
        for l, bstate in enumerate(self.levels):
            cg = bstate.g
            if l == 0:
                cg = GraphView(cg, skip_properties=True)
            cg.vp["b"] = bstate.b.copy()
            cg.ep["count"] = bstate.eweight
            if bstate.overlap:
                cg.vp["node_index"] = bstate.node_index.copy()

            bstack.append(cg)
            if bstate.N == 1:
                break
        if bstack[-1].num_vertices() > 1:
            cg = Graph(directed=bstack[-1].is_directed())
            cg.add_vertex()
            cg.vp["b"] = cg.new_vertex_property("int")
            e = cg.add_edge(0, 0)
            ew = cg.new_edge_property("int")
            ew[e] = self.levels[-1].E
            cg.ep["count"] = ew
            bstack.append(cg)
        return bstack


    def project_level(self, l):
        r"""Project the partition at level ``l`` onto the lowest level, and return the
        corresponding :class:`~graph_tool.community.BlockState` (or
        :class:`~graph_tool.community.OverlapBlockState`).  """
        if self.overlap != "full":
            clabel = b = self.levels[l].b.copy()
            while l - 1 >= 0:
                clabel = b
                b = self.levels[l - 1].b.copy()
                pmap(b, clabel)
                l -= 1
        else:
            b = self.levels[l].b.copy()
        state = self.levels[0].copy(b=b.a)
        return state

    def print_summary(self):
        for l, state in enumerate(self.levels):
            print("l: %d, N: %d, B: %d" % (l, state.N, state.B))


def nested_mcmc_sweep(state, beta=1., c=1., dl=False, sequential=True, verbose=False):
    r"""Performs a Markov chain Monte Carlo sweep on all levels of the hierarchy.

    Parameters
    ----------
    state : :class:`~graph_tool.community.NestedBlockState`
        The nested block state.
    beta : `float` (optional, default: `1.0`)
        The inverse temperature parameter :math:`\beta`.
    c : ``float`` (optional, default: ``1.0``)
        This parameter specifies how often fully random moves are attempted,
        instead of more likely moves based on the inferred block partition.
        For ``c == 0``, no fully random moves are attempted, and for ``c == inf``
        they are always attempted.
    dl : ``bool`` (optional, default: ``False``)
        If ``True``, the change in the whole description length will be
        considered after each vertex move, not only the entropy.
    sequential : ``bool`` (optional, default: ``True``)
        If ``True``, the move attempts on the vertices are done in sequential
        random order. Otherwise a total of `N` moves attempts are made, where
        `N` is the number of vertices, where each vertex can be selected with
        equal probability.
    verbose : ``bool`` (optional, default: ``False``)
        If ``True``, verbose information is displayed.

    Returns
    -------

    dS_moves : list of (``float``, ``int``) tuples
       The entropy difference (per edge) and number of accepted block membership
       moves after a full sweep for each level.


    Notes
    -----

    This algorithm performs a Markov chain Monte Carlo sweep on each level of the
    network, via the function :func:`~graph_tool.community.mcmc_sweep`.

    This algorithm has a worse-case complexity of :math:`O(E \times L)`, where
    :math:`E` is the number of edges in the network, and :math:`L` is the depth
    of the hierarchy.

    Examples
    --------
    .. testsetup:: nested_mcmc

       gt.seed_rng(42)
       np.random.seed(42)

    .. doctest:: nested_mcmc

       >>> g = gt.collection.data["polbooks"]
       >>> state = gt.NestedBlockState(g, Bs=[10, 5, 3, 2, 1], deg_corr=True)
       >>> ret = gt.nested_mcmc_sweep(state)
       >>> print(ret)
       [(0.0, 0), (0.0, 0), (0.0, 0), (0.0, 0), (0.0, 0)]

    References
    ----------

    .. [peixoto-efficient-2014] Tiago P. Peixoto, "Efficient Monte Carlo and greedy
       heuristic for the inference of stochastic block models", Phys. Rev. E 89, 012804 (2014),
       :doi:`10.1103/PhysRevE.89.012804`, :arxiv:`1310.4378`.
    .. [peixoto-hierarchical-2014] Tiago P. Peixoto, "Hierarchical block structures and high-resolution
       model selection in large networks ", Phys. Rev. X 4, 011047 (2014), :doi:`10.1103/PhysRevX.4.011047`,
       :arxiv:`1310.4377`.
    .. [peixoto-model-2014] Tiago P. Peixoto, "Model selection and hypothesis
       testing for large-scale network models with overlapping groups",
       :arxiv:`1409.3059`.
    """

    rets = []
    for l, bstate in enumerate(state.levels):
        if verbose:
            print("Level:", l, "N:", bstate.N, "B:", bstate.B)

        # constraint partitions not to invalidate upper layers
        if l < len(state.levels) - 1:
            clabel = state._NestedBlockState__project_partition(l, l + 1)
        else:
            clabel = bstate.g.new_vertex_property("int")

        # propagate externally imposed clabel at the bottom
        cclabel = state._NestedBlockState__propagate_clabel(l)
        cclabel.a += clabel.a * (cclabel.a.max() + 1)
        continuous_map(cclabel)

        bstate.clabel = cclabel

        ret = mcmc_sweep(bstate, beta=beta, c=c, dl=dl,
                         dense = l > 0 and state.deg_corr != "full",
                         multigraph = l > 0,
                         sequential=sequential, verbose=verbose)

        bstate.clabel.a = 0

        rets.append(ret)
    return rets

def replace_level(l, state, min_B=None, max_B=None, max_b=None, nsweeps=10,
                  nmerge_sweeps=10, r=2, c=0, epsilon=0., sequential=True,
                  dl=False, dense=False, multigraph=True, sparse_thresh=100,
                  verbose=False, checkpoint=None, minimize_state=None,
                  dl_ent=False):
    r"""Replaces level l with another state with a possibly different number of
    groups. This may change not only the state at level l, but also the one at
    level l + 1, which needs to be 'rebuilt' because of the label changes at
    level l."""

    if __test__:
        state._NestedBlockState__consistency_check("(before) replace level", l)

    bstate = state.levels[l]
    g = bstate.g
    base_g = g if not bstate.overlap else bstate.base_g
    eweight = bstate.eweight if not bstate.overlap else None

    if l > 0 or min_B is None:
        if l + 1 < len(state.levels):
            min_B = state.levels[l + 1].B
        else:
            min_B = 1
    if l > 0 or max_B is None:
        max_B = bstate.N

    min_B = max(min_B, state.clabel.a.max() + 1)

    # constraint partitions not to invalidate upper layers
    if l < len(state.levels) - 1:
        clabel = state._NestedBlockState__project_partition(l, l + 1)
    else:
        clabel = bstate.g.new_vertex_property("int")

    assert min_B <= max_B, (min_B, max_B, bstate.B, bstate.N, g.num_vertices(), state.clabel.a.max() + 1, clabel.a.max() + 1, l)

    # propagate externally imposed clabel at the bottom
    cclabel = state._NestedBlockState__propagate_clabel(l)
    assert cclabel.a.max() <= state.clabel.a.max(),  (cclabel.a.max(), state.clabel.a.max())
    cclabel.a += clabel.a * (cclabel.a.max() + 1)
    continuous_map(cclabel)

    min_B = max(min_B, cclabel.a.max() + 1)

    assert min_B <= max_B, (min_B, max_B, bstate.B, bstate.N, g.num_vertices(), cclabel.a.max() + 1, state.clabel.a.max() + 1, clabel.a.max() + 1, l)

    if __test__:
        assert bstate._BlockState__check_clabel(), "invalid clabel before minimize!"

    nested_dl = l < len(state.levels) - 1

    Sb = get_b_dl(state.levels[l],
                  dense=(l > 0 and state.deg_corr != "full") or dense,
                  multigraph=l > 0 or multigraph,
                  nested_dl=nested_dl,
                  nested_overlap=state.overlap == "full",
                  dl_ent=dl_ent)

    if __test__:
        assert clabel.a.max() + 1 <= min_B

    if l == 0 and max_b is not None:
        istate = state.levels[0].copy(b=max_b.copy(),
                                      clabel=cclabel.a if state.overlap else cclabel)
        assert istate._BlockState__check_clabel(), "invalid clabel!"
        if istate.B > max_B:
            istate = multilevel_minimize(istate, B=max_B, nsweeps=nsweeps,
                                         nmerge_sweeps=nmerge_sweeps, c=c,
                                         r=r, dl=dl, sequential=sequential,
                                         adaptive_sweeps=True,
                                         greedy_cooling=True,
                                         epsilon=epsilon,
                                         verbose=verbose=="full")
        init_states = [istate]
    else:
        init_states = None

    res = minimize_blockmodel_dl(base_g, eweight=eweight,
                                 deg_corr=bstate.deg_corr,
                                 nsweeps=nsweeps,
                                 nmerge_sweeps=nmerge_sweeps,
                                 c=c, r=r, sequential=sequential,
                                 adaptive_sweeps=True, greedy_cooling=True,
                                 epsilon=epsilon,
                                 max_B=max_B,
                                 min_B=min_B,
                                 clabel=cclabel if not bstate.overlap else cclabel.a,
                                 max_BE=bstate.max_BE,
                                 dl=dl,
                                 dense=(l > 0 and state.deg_corr != "full") or dense,
                                 multigraph=l > 0 or multigraph,
                                 sparse_heuristic=base_g.num_vertices() > sparse_thresh,
                                 nested_dl=nested_dl,
                                 overlap=bstate.overlap,
                                 nested_overlap=state.overlap == "full",
                                 nonoverlap_compare=False,
                                 nonoverlap_init=False,
                                 init_states=init_states,
                                 verbose=verbose=="full",
                                 ##exaustive=g.num_vertices() <= 100,
                                 #minimize_state=minimize_state.minimize_state, >>>>>> HERE <<<<<
                                 checkpoint=checkpoint,
                                 dl_ent=dl_ent)

    if __test__:
        assert (res.clabel.a == cclabel.a).all(), (res.clabel.a, cclabel.a)
        assert res._BlockState__check_clabel(), "invalid clabel after minimize!"

    res.clabel.a = 0
    b = res.b

    Sf = get_b_dl(res,
                  dense=(l > 0 and state.deg_corr != "full") or dense,
                  multigraph=l > 0 or multigraph,
                  nested_dl=nested_dl,
                  nested_overlap=state.overlap == "full",
                  dl_ent=dl_ent)

    kept = False
    if Sf - Sb >= -1e-10:
        kept = True
    if Sf - Sb == 0 and bstate.B != state.levels[l].B:
        kept = False

    if kept:
        Sf_rej = Sf
        Sf = Sb
    else:
        # rebuild current level
        state._NestedBlockState__rebuild_level(l, b=b, clabel=clabel)

    if verbose:
        print("level", l, ": resizing", bstate.B, "->", state.levels[l].B, ", dS:", Sf - Sb, end="")
        if kept:
            print(" [kept, rejected (%d, %g) vs (%d, %g)]" % (res.B, Sf_rej, bstate.B, Sb))
        else:
            print()

    if __test__:
        state._NestedBlockState__consistency_check("replace level", l)

    dS = Sf - Sb
    return dS, kept


class NestedMinimizeState(object):
    r"""This object stores information regarding the current entropy minimization
    state, so that the algorithms can resume previously started runs.
    This object can be saved to disk via the :mod:`pickle` interface."""

    def __init__(self):
        self.minimize_state = MinimizeState()
        self.l = 0
        self.bs = []
        self.done = []
        self.init = True

    def clear(self):
        self.minimize_state.clear()
        self.l = 0
        del self.bs[:]
        del self.done[:]

    def sync(self, state):
        if len(self.bs) == 0:
            for s in state.levels:
                self.bs.append(array(s.b.fa))
        while len(self.done) < len(state.levels):
            self.done.append(False)

    def delete(self, l):
        del self.done[l]
        del self.bs[l]

    def insert(self, l, state):
        self.done.insert(l + 1, False)
        ba = array(state.levels[l].b.fa)
        self.bs.insert(l + 1, ba)

    def mark_level(self, l, done, state):
        while len(state.levels) > len(self.bs):
            ba = array(state.levels[len(self.bs)].b.fa)
            self.bs.append(ba)
            self.done.append(False)

        self.done[l] = done
        if done:
            self.bs[l] = array(state.levels[l].b.fa)

    def clear_mstate(self):
        self.minimize_state.clear()

def get_checkpoint_wrap(checkpoint, state, minimize_state, dl_ent):
    S_total = state.entropy(complete=True, dl_ent=dl_ent)

    if checkpoint is not None:
        def check_wrap(bstate, Sb, delta, nmoves, ms):
            l = minimize_state.l
            bstate = state.levels[l]
            S_l = bstate.entropy()
            S = S_total - S_l + Sb
            if bstate is None:
                checkpoint(None, S, delta, nmoves, minimize_state)
            else:
                checkpoint(state, S, delta, nmoves, minimize_state)
        chkp = check_wrap
    else:
        chkp = None

    return chkp


def nested_tree_sweep(state, min_B=None, max_B=None, max_b=None, nsweeps=10,
                      epsilon=0., r=2., nmerge_sweeps=10, c=0, dl=False,
                      dense=False, multigraph=True, sequential=True,
                      sparse_thresh=100, checkpoint=None, minimize_state=None,
                      frozen_levels=None, verbose=False, **kwargs):
    r"""Performs one greedy sweep in the entire hierarchy tree, attempting to
    decrease its description length.

    Parameters
    ----------
    state : :class:`~graph_tool.community.NestedBlockState`
        The nested block state.
    min_B : ``int`` (optional, default: ``None``)
        Minimum number of blocks at the lowest level.
    max_B : ``int`` (optional, default: ``None``)
        Maximum number of blocks at the lowest level.
    max_b : ``int`` (optional, default: ``None``)
        Block partition used for the maximum number of blocks at the lowest
        level.
    nsweeps : ``int`` (optional, default: ``10``)
        The number of sweeps done after each merge step to reach the local
        minimum.
    epsilon : ``float`` (optional, default: ``0``)
        Converge criterion for ``adaptive_sweeps``.
    r : ``float`` (optional, default: ``2.``)
        Agglomeration ratio for the merging steps. Each merge step will attempt
        to find the best partition into :math:`B_{i-1} / r` blocks, where
        :math:`B_{i-1}` is the number of blocks in the last step.
    nmerge_sweeps : `int` (optional, default: `10`)
        The number of merge sweeps done, where in each sweep a better merge
        candidate is searched for every block.
    c : ``float`` (optional, default: ``1.0``)
        This parameter specifies how often fully random moves are attempted,
        instead of more likely moves based on the inferred block partition.
        For ``c == 0``, no fully random moves are attempted, and for ``c == inf``
        they are always attempted.
    dense : ``bool`` (optional, default: ``False``)
        If ``True``, the "dense" variant of the entropy will be computed.
    sequential : ``bool`` (optional, default: ``True``)
        If ``True``, the move attempts on the vertices are done in sequential
        random order. Otherwise a total of `N` moves attempts are made, where
        `N` is the number of vertices, where each vertex can be selected with
        equal probability.
    sparse_thresh : ``int`` (optional, default: ``100``)
        If the number of nodes in the higher level multigraphs exceeds this
        number, the sparse entropy will be used to find the best partition,
        but the dense entropy will be used to compare different partitions.
    checkpoint : function (optional, default: ``None``)
        If provided, this function will be called after each call to
        :func:`mcmc_sweep`. This can be used to store the current state, so it
        can be continued later. The function must have the following signature:

        .. code-block:: python

            def checkpoint(state, S, delta, nmoves, minimize_state):
                ...

        where `state` is either a :class:`~graph_tool.community.NestedBlockState`
        instance or ``None``, `S` is the current entropy value, `delta` is
        the entropy difference in the last MCMC sweep, and `nmoves` is the
        number of accepted block membership moves. The ``minimize_state``
        argument is a :class:`NestedMinimizeState` instance which specifies
        the current state of the algorithm, which can be stored via :mod:`pickle`,
        and supplied via the ``minimize_state`` option below to continue from an
        interrupted run.

        This function will also be called when the MCMC has finished for the
        current value of :math:`B`, in which case ``state == None``, and the
        remaining parameters will be zero, except the last.
    minimize_state : :class:`NestedMinimizeState` (optional, default: ``None``)
        If provided, this will specify an exact point of execution from which
        the algorithm will continue. The expected object is a
        :class:`NestedMinimizeState` instance which will be passed to the
        callback of the ``checkpoint`` option above, and  can be stored by
        :mod:`pickle`.
    frozen_levels : :class:`list` (optional, default: ``None``)
        List of levels (``int``s) which will remain unmodified during the
        algorithm.
    verbose : ``bool`` (optional, default: ``False``)
        If ``True``, verbose information is displayed.

    Returns
    -------

    dS : ``float``
       The description length difference (per edge) after the move.


    Notes
    -----

    This algorithm performs a constrained agglomerative heuristic on each level
    of the network, via the function :func:`~graph_tool.community.multilevel_minimize`.

    This algorithm has worst-case complexity of :math:`O(N\ln^2 N \times L)`,
    where  :math:`N` is the number of nodes in the network, and :math:`L` is
    the depth of the hierarchy.

    References
    ----------

    .. [peixoto-efficient-2014] Tiago P. Peixoto, "Efficient Monte Carlo and greedy
       heuristic for the inference of stochastic block models", Phys. Rev. E 89, 012804 (2014),
       :doi:`10.1103/PhysRevE.89.012804`, :arxiv:`1310.4378`.
    .. [peixoto-hierarchical-2014] Tiago P. Peixoto, "Hierarchical block structures and high-resolution
       model selection in large networks ", Phys. Rev. X 4, 011047 (2014), :doi:`10.1103/PhysRevX.4.011047`,
       :arxiv:`1310.4377`.
    .. [peixoto-model-2014] Tiago P. Peixoto, "Model selection and hypothesis
       testing for large-scale network models with overlapping groups",
       :arxiv:`1409.3059`.
    """

    dl_ent = kwargs.get("dl_ent", False)

    if minimize_state is None:
        minimize_state = NestedMinimizeState()
    mstate = minimize_state
    mstate.sync(state)

    args = dict(state=state, nsweeps=nsweeps, nmerge_sweeps=nmerge_sweeps, r=r,
                c=c, epsilon=epsilon, sequential=sequential, dl=dl, dense=dense,
                multigraph=multigraph, sparse_thresh=sparse_thresh,
                min_B=min_B, max_B=max_B, max_b=max_b, checkpoint=checkpoint,
                minimize_state=minimize_state, dl_ent=dl_ent)

    #_Si = state.entropy(dense=dense, multigraph=dense)
    dS = 0

    if frozen_levels is None:
        frozen_levels = set()

    while mstate.l >= 0:
        l = mstate.l

        if mstate.done[l]:
            if verbose:
                print("level", l, ": skipping", state.levels[l].B)
            mstate.l -= 1
            continue

        Si = state.entropy(dl_ent=dl_ent)

        kept = True

        if l in frozen_levels:
            kept = False

        # replace level
        if kept:
            ddS, kept = replace_level(l, verbose=verbose, **args)
            dS += ddS
            mstate.clear_mstate()

        if __test__:
            if kept:
                assert abs(state.entropy(dl_ent=dl_ent) - Si) < 1e-6, "inconsistent replace at level %d (%g, %g)" % (l, state.entropy(), Si)
            state._NestedBlockState__consistency_check("replace level", l)

        # delete level
        if (kept and l > 0 and l < len(state.levels) - 1 and
            not (min_B is not None and l == 1 and state.levels[l].B < min_B)):
            Si = state.entropy(dl_ent=dl_ent)

            bstates = [state.levels[l-1], state.levels[l], state.levels[l + 1]]

            state._NestedBlockState__delete_level(l)
            #replace_level(l, **args)

            Sf = state.entropy(dl_ent=dl_ent)

            mstate.clear_mstate()
            if Sf > Si:
                state.levels[l - 1] = bstates[0]
                state.levels.insert(l, bstates[1])
                state.levels[l + 1] = bstates[2]
            else:
                kept = False
                dS += Sf - Si

                mstate.delete(l)

                if verbose:
                    print("level", l, ": deleted", (bstates[1].N, bstates[1].B), ", dS:", Sf - Si, len(state.levels))

            if __test__:
                if kept:
                    assert abs(state.entropy(dl_ent=dl_ent) - Si) < 1e-6, "inconsistent delete at level %d (%g, %g)" % (l, state.entropy(), Si)
                state._NestedBlockState__consistency_check("delete complete", l)

        # insert new level (duplicate and replace)
        if kept and l > 0:
            Si = state.entropy(dl_ent=dl_ent)

            bstates = [state.levels[l].copy()]
            if l < len(state.levels) - 1:
                bstates.append(state.levels[l + 1].copy())
            if l < len(state.levels) - 2:
                bstates.append(state.levels[l + 2].copy())

            state._NestedBlockState__duplicate_level(l)

            replace_level(l + 1, verbose=False, **args)

            Sf = state.entropy(dl_ent=dl_ent)

            mstate.clear_mstate()
            if Sf >= Si:
                del state.levels[l + 1]
                for j in range(len(bstates)):
                    state.levels[l + j] = bstates[j]
                if bstates[-1].B == 1:
                    del state.levels[l + len(bstates):]
            else:
                kept = False
                dS += Sf - Si

                mstate.insert(l, state)

                l += 1

                if verbose:
                    print("level", l, ": inserted", state.levels[l].B, ", dS:", Sf - Si)

            if __test__:
                state._NestedBlockState__consistency_check("delete", l)
                if kept:
                    assert abs(state.entropy(dl_ent=dl_ent) - Si) < 1e-8, "inconsistent delete at level %d (%g, %g)" % (l, state.entropy(), Si)

        mstate.mark_level(l, done=True, state=state)
        if not kept:
            if l + 1 < len(state.levels):
                mstate.mark_level(l + 1, done=False, state=state)
            if l > 0:
                mstate.mark_level(l - 1, done=False, state=state)
            l += 1
        else:
            if ((l + 1 < len(state.levels) and not mstate.done[l + 1]) or
                (l + 1 == len(state.levels) and state.levels[l].B > 1)):
                l += 1
            else:
                l -= 1

        if l >= len(state.levels):
            l = len(state.levels) - 1

        # create a new level at the top with B=1, if necessary
        if l == len(state.levels) - 1 and state.levels[l].B > 1:
            NB = state.levels[l].B if not state.overlap else 2 * state.levels[l].E
            state._NestedBlockState__rebuild_level(l, b=state.levels[l].g.new_vertex_property("int"))
            mstate.mark_level(l + 1, done=False, state=state)
            l += 1

        mstate.l = l
        if checkpoint is not None:
            checkpoint(None, 0, 0, 0, mstate)

        if __test__:
            state._NestedBlockState__consistency_check("tree sweep step", l)

    # _Sf = state.entropy(dense=dense, multigraph=dense, dl_ent=dl_ent)

    return dS



def init_nested_state(g, Bs, deg_corr=True, overlap=False, dl=False, dense=False,
                      multigraph=True, eweight=None, vweight=None,
                      clabel=None, nsweeps=10, epsilon=0., r=2, nmerge_sweeps=10,
                      c=0, sequential=True, sparse_thresh=100,
                      checkpoint=None, minimize_state=None, max_BE=1000,
                      verbose=False, **kwargs):
    r"""Initializes a nested block hierarchy with sizes given by ``Bs``.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        The graph being modelled.
    Bs : list of ``int`` (optional, default: ``None``)
        Number of blocks for each hierarchy level.
    deg_corr : ``bool`` (optional, default: ``True``)
        If ``True``, the degree-corrected version of the blockmodel ensemble will
        be used in the bottom level, otherwise the traditional variant will be used.
    dense : ``bool`` (optional, default: ``False``)
        If ``True``, the "dense" variant of the entropy will be computed.
    eweight : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Edge weights (i.e. multiplicity).
    vweight : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Vertex weights (i.e. multiplicity).
    nsweeps : ``int`` (optional, default: ``10``)
        The number of sweeps done after each merge step to reach the local
        minimum.
    epsilon : ``float`` (optional, default: ``0``)
        Converge criterion for ``adaptive_sweeps``.
    r : ``float`` (optional, default: ``2.``)
        Agglomeration ratio for the merging steps. Each merge step will attempt
        to find the best partition into :math:`B_{i-1} / r` blocks, where
        :math:`B_{i-1}` is the number of blocks in the last step.
    nmerge_sweeps : `int` (optional, default: `10`)
        The number of merge sweeps done, where in each sweep a better merge
        candidate is searched for every block.
    c : ``float`` (optional, default: ``0.``)
        This parameter specifies how often fully random moves are attempted,
        instead of more likely moves based on the inferred block partition.
        For ``c == 0``, no fully random moves are attempted, and for ``c == inf``
        they are always attempted.
    sequential : ``bool`` (optional, default: ``True``)
        If ``True``, the move attempts on the vertices are done in sequential
        random order. Otherwise a total of `N` moves attempts are made, where
        `N` is the number of vertices, where each vertex can be selected with
        equal probability.
    sparse_thresh : ``int`` (optional, default: ``100``)
        If the number of nodes in the higher level multigraphs exceeds this
        number, the sparse entropy will be used to find the best partition,
        but the dense entropy will be used to compare different partitions.
    checkpoint : function (optional, default: ``None``)
        If provided, this function will be called after each call to
        :func:`mcmc_sweep`. This can be used to store the current state, so it
        can be continued later. The function must have the following signature:

        .. code-block:: python

            def checkpoint(state, S, delta, nmoves, minimize_state):
                ...

        where `state` is either a :class:`~graph_tool.community.NestedBlockState`
        instance or ``None``, `S` is the current entropy value, `delta` is
        the entropy difference in the last MCMC sweep, and `nmoves` is the
        number of accepted block membership moves. The ``minimize_state``
        argument is a :class:`NestedMinimizeState` instance which specifies
        the current state of the algorithm, which can be stored via :mod:`pickle`,
        and supplied via the ``minimize_state`` option below to continue from an
        interrupted run.

        This function will also be called when the MCMC has finished for the
        current value of :math:`B`, in which case ``state == None``, and the
        remaining parameters will be zero, except the last.
    minimize_state : :class:`NestedMinimizeState` (optional, default: ``None``)
        If provided, this will specify an exact point of execution from which
        the algorithm will continue. The expected object is a
        :class:`NestedMinimizeState` instance which will be passed to the
        callback of the ``checkpoint`` option above, and  can be stored by
        :mod:`pickle`.
    verbose : ``bool`` (optional, default: ``False``)
        If ``True``, verbose information is displayed.

    Returns
    -------

    state : :class:`~graph_tool.community.NestedBlockState`
       The initialized nested state.


    Notes
    -----

    This algorithm performs an agglomerative heuristic on each level  of the
    network, via the function :func:`~graph_tool.community.multilevel_minimize`.

    This algorithm has worst-case complexity of :math:`O(N\ln^2 N \times L)`,
    where  :math:`N` is the number of nodes in the network, and :math:`L` is
    the depth of the hierarchy.

    References
    ----------

    .. [peixoto-efficient-2014] Tiago P. Peixoto, "Efficient Monte Carlo and greedy
       heuristic for the inference of stochastic block models", Phys. Rev. E 89, 012804 (2014),
       :doi:`10.1103/PhysRevE.89.012804`, :arxiv:`1310.4378`.
    .. [peixoto-hierarchical-2014] Tiago P. Peixoto, "Hierarchical block structures and high-resolution
       model selection in large networks ", Phys. Rev. X 4, 011047 (2014), :doi:`10.1103/PhysRevX.4.011047`,
       :arxiv:`1310.4377`.
    .. [peixoto-model-2014] Tiago P. Peixoto, "Model selection and hypothesis
       testing for large-scale network models with overlapping groups",
       :arxiv:`1409.3059`.
    """

    dl_ent = kwargs.get("dl_ent", False)

    if minimize_state is None:
        minimize_state = NestedMinimizeState()
    mstate = minimize_state

    state = NestedBlockState(g, eweight=eweight, vweight=vweight, Bs=[1],
                             deg_corr=deg_corr, overlap=overlap, clabel=clabel)

    chkp = get_checkpoint_wrap(checkpoint, state, minimize_state, dl_ent)

    bg = g
    ecount = eweight

    for l, B in enumerate(Bs):
        ba = None
        if l < len(mstate.bs):
            ba = mstate.bs[l]
        else:
            if l == 0:
                if state.overlap:
                    bstate = OverlapBlockState(bg, B=bg.num_vertices(), #b=bg.vertex_index.copy("int"),
                                               vweight=vweight,
                                               eweight=ecount,
                                               deg_corr=deg_corr != False,
                                               #clabel=clabel,
                                               max_BE=max_BE)
                else:
                    bstate = BlockState(bg, B=bg.num_vertices(), #b=bg.vertex_index.copy("int"),
                                        vweight=vweight,
                                        eweight=ecount,
                                        deg_corr=deg_corr != False,
                                        #clabel=clabel,
                                        max_BE=max_BE)
            else:
                bstate = state.levels[l-1].get_block_state(b=ba,
                                                           overlap=overlap == "full",
                                                           deg_corr=deg_corr == "full")[0]

            if B == 1:
                bstate.copy(b=bstate.g.new_vertex_property("int").a)
            else:
                bstate = multilevel_minimize(bstate, B, nsweeps=nsweeps, epsilon=epsilon,
                                             r=r, nmerge_sweeps=nmerge_sweeps,
                                             greedy=True, c=c, dl=dl,
                                             dense=(l > 0 and g.num_vertices() < sparse_thresh) or dense,
                                             multigraph=(l > 0 and g.num_vertices() < sparse_thresh) or multigraph,
                                             sequential=sequential, verbose=verbose,
                                             checkpoint=chkp,
                                             minimize_state=minimize_state.minimize_state)
            ba = array(bstate.b.fa)
            mstate.bs.append(ba)
            minimize_state.clear_mstate()

        state._NestedBlockState__rebuild_level(len(state.levels) - 1, b=ba)
        bg = state.levels[l].bg
        ecount = state.levels[l].mrs


    for l, B in enumerate(Bs):
        if l + 1 < len(state.levels):
            assert state.levels[l].B == state.levels[l + 1].N

    minimize_state.clear()
    mstate.sync(state)
    if checkpoint is not None:
        checkpoint(None, 0, 0, 0, mstate)

    return state

def minimize_nested_blockmodel_dl(g, Bs=None, bs=None, min_B=None, max_B=None,
                                  max_b=None, deg_corr=True, overlap=False,
                                  nonoverlap_init=True, dl=True,
                                  multigraph=True, dense=False, eweight=None,
                                  vweight=None, clabel=None, frozen_levels=None,
                                  nsweeps=10, epsilon=1e-3, c=0,
                                  nmerge_sweeps=10, r=2, sparse_thresh=100,
                                  sequential=True, verbose=False,
                                  checkpoint=None, minimize_state=None,
                                  **kwargs):
    r"""Find the block hierarchy of an unspecified size which minimizes the description
    length of the network, according to the nested stochastic blockmodel ensemble which
    best describes it.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph being used.
    Bs : list of ``int`` (optional, default: ``None``)
        Initial number of blocks for each hierarchy level.
    bs : list of :class:`~graph_tool.PropertyMap` or :class:`~numpy.ndarray` instances (optional, default: ``None``)
        Initial block labels on the vertices, for each hierarchy level.
    min_B : ``int`` (optional, default: ``None``)
        Minimum number of blocks at the lowest level.
    max_B : ``int`` (optional, default: ``None``)
        Maximum number of blocks at the lowest level.
    max_b : ``int`` (optional, default: ``None``)
        Block partition used for the maximum number of blocks at the lowest
        level.
    deg_corr : ``bool`` (optional, default: ``True``)
        If ``True``, the degree-corrected version of the blockmodel ensemble
        will be used in the bottom level, otherwise the traditional variant will
        be used.
    overlap : ``bool`` (optional, default: ``False``)
        If ``True``, the mixed-membership version of the blockmodel will be used.
    nonoverlap_init : ``bool`` (optional, default: ``True``)
        If ``True``, and `´overlap == True``, the minimization starts by first
        fitting the non-overlapping model, and using that as a starting state.
    dl : ``bool`` (optional, default: ``True``)
        If ``True``, the change in the whole description length will be
        considered after each vertex move, not only the entropy.
    multigraph : ``bool`` (optional, default: ``False``)
            If ``True``, the multigraph entropy will be used.
    dense : ``bool`` (optional, default: ``False``)
        If ``True``, the "dense" variant of the entropy will be computed.
    eweight : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Edge weights (i.e. multiplicity).
    vweight : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Vertex weights (i.e. multiplicity).
    clabel : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Constraint labels on the vertices. If supplied, vertices with different
        label values will not be clustered in the same group.
    frozen_levels : :class:`list` (optional, default: ``None``)
        List of levels (``int``s) which will remain unmodified during the
        algorithm.
    nsweeps : ``int`` (optional, default: ``10``)
        The number of sweeps done after each merge step to reach the local
        minimum.
    epsilon : ``float`` (optional, default: ``0``)
        The number of sweeps necessary for the local minimum will
        be estimated to be enough so that no more than ``epsilon * N`` nodes
        changes their states in the last ``nsweeps`` sweeps.
    c : ``float`` (optional, default: ``0.``)
        This parameter specifies how often fully random moves are attempted,
        instead of more likely moves based on the inferred block partition.
        For ``c == 0``, no fully random moves are attempted, and for ``c == inf``
        they are always attempted.
    nmerge_sweeps : `int` (optional, default: `10`)
        The number of merge sweeps done, where in each sweep a better merge
        candidate is searched for every block.
    r : ``float`` (optional, default: ``2.``)
        Agglomeration ratio for the merging steps. Each merge step will attempt
        to find the best partition into :math:`B_{i-1} / r` blocks, where
        :math:`B_{i-1}` is the number of blocks in the last step.
    sparse_thresh : ``int`` (optional, default: ``100``)
        If the number of blocks at some level is larger than this value, the
        sparse entropy will be used to find the best partition, but the dense
        entropy will be used to compare different partitions.
    sequential : ``bool`` (optional, default: ``True``)
        If ``True``, the move attempts on the vertices are done in sequential
        random order. Otherwise a total of `N` moves attempts are made, where
        `N` is the number of vertices, where each vertex can be selected with
        equal probability.
    checkpoint : function (optional, default: ``None``)
        If provided, this function will be called after each call to
        :func:`mcmc_sweep`. This can be used to store the current state, so it
        can be continued later. The function must have the following signature:

        .. code-block:: python

            def checkpoint(state, L, delta, nmoves, minimize_state):
                ...

        where `state` is either a :class:`~graph_tool.community.NestedBlockState`
        instance or ``None``, `L` is the current description length, `delta` is
        the entropy difference in the last MCMC sweep, and `nmoves` is the
        number of accepted block membership moves. The ``minimize_state``
        argument is a :class:`~graph_tool.community.NestedMinimizeState`
        instance which specifies the current state of the algorithm, which can
        be stored via :mod:`pickle`, and supplied via the ``minimize_state``
        option below to continue from an interrupted run.

        This function will also be called when the MCMC has finished for the
        current value of :math:`B`, in which case ``state == None``, and the
        remaining parameters will be zero, except the last.
    minimize_state : :class:`MinimizeState` (optional, default: ``None``)
        If provided, this will specify an exact point of execution from which
        the algorithm will continue. The expected object is a
        :class:`~graph_tool.community.NestedMinimizeState`
        instance which will be passed to the callback of the ``checkpoint``
        option above, and  can be stored by :mod:`pickle`.
    verbose : ``bool`` (optional, default: ``False``)
        If ``True``, verbose information is displayed.

    Returns
    -------
    state : :class:`~graph_tool.community.NestedBlockState`
        The nested block state.

    Notes
    -----

    This algorithm attempts to find a block partition hierarchy of an unspecified size
    which minimizes the description length of the network,

    .. math::

       \Sigma = \mathcal{L}_{t/c} + \mathcal{S}_n,

    where :math:`\mathcal{S}_{n}` is the nested blockmodel entropy given by

    .. math::

       \mathcal{S}_n  = \mathcal{S}_{t/c}(\{e^0_{rs}\}, \{n^0_r\}) + \sum_{l=1}^LS_m(\{e^l_{rs}\}, \{n^l_r\}).

    with :math:`\mathcal{S}_{t/c}` and :math:`\mathcal{S}_{m}` described in the docstring of
    :meth:`~graph_tool.community.BlockState.entropy`, and :math:`\{e^l_{rs}\}`
    and :math:`\{n^l_r\}` are the edge and node counts at hierarchical level :math:`l`.
    Additionally :math:`\mathcal{L}_{t/c}` is the information necessary to
    describe the block partitions, i.e. :math:`\mathcal{L}_t=\sum_{l=0}^L\mathcal{L}^l_t`, with

    .. math::

        \mathcal{L}^l_t = \ln\left(\!\!{B_l\choose B_{l-1}}\!\!\right) + \ln B_{l-1}! - \sum_r \ln n_r^l!.


    See [peixoto-hierarchical-2014]_ for details on the algorithm.

    This algorithm has a complexity of :math:`O(N \ln^2 N)`, where :math:`N`
    is the number of nodes in the network.

    Examples
    --------
    .. testsetup:: nested_mdl

       gt.seed_rng(42)
       np.random.seed(42)

    .. doctest:: nested_mdl

       >>> g = gt.collection.data["power"]
       >>> state = gt.minimize_nested_blockmodel_dl(g, deg_corr=True)
       >>>
       >>> # route edges according to the hierarchy
       >>> bstack = state.get_bstack()
       >>> t = gt.get_hierarchy_tree(bstack)[0]
       >>> tpos = pos = gt.radial_tree_layout(t, t.vertex(t.num_vertices() - 1), weighted=True)
       >>> cts = gt.get_hierarchy_control_points(g, t, tpos)
       >>> pos = g.own_property(tpos)
       >>>
       >>> b = bstack[0].vp["b"]
       >>> gt.graph_draw(g, pos=pos, vertex_fill_color=b, vertex_shape=b, edge_control_points=cts,
       ...               edge_color=[0, 0, 0, 0.3], vertex_anchor=0, output="power_nested_mdl.pdf")
       <...>

    .. testcleanup:: nested_mdl

       gt.graph_draw(g, pos=pos, vertex_fill_color=b, vertex_shape=b, edge_control_points=cts, edge_color=[0, 0, 0, 0.3], vertex_anchor=0, output="power_nested_mdl.png")

    .. figure:: power_nested_mdl.*
       :align: center

       Block partition of a power-grid network, which minimizes the description
       length of the network according to the nested (degree-corrected) stochastic blockmodel.


    .. doctest:: nested_mdl_overlap

       >>> g = gt.collection.data["power"]
       >>> state = gt.minimize_nested_blockmodel_dl(g, deg_corr=True, overlap=True, dl=False)
       >>> bv, *rest, bc = state.levels[0].get_overlap_blocks()
       >>>
       >>> # get the dominant group for each node for the layout
       >>> n_r = np.zeros(state.levels[0].B)
       >>> b = g.new_vertex_property("int")
       >>> d = g.new_vertex_property("int")
       >>> for v in g.vertices():
       ...     i = bc[v].a.argmax()
       ...     b[v] = bv[v][i]
       ...     n_r[b[v]] += 1
       ...     d[v] = len(bv[v])
       >>> g.vp["b"] = b
       >>> bstack = state.get_bstack()
       >>> bstack[0] = g
       >>>
       >>> # route edges according to the hierarchy
       >>> t = gt.get_hierarchy_tree(bstack)[0]
       >>> tpos = pos = gt.radial_tree_layout(t, t.vertex(t.num_vertices() - 1), weighted=True)
       >>> cts = gt.get_hierarchy_control_points(g, t, tpos)
       >>> pos = g.own_property(tpos)
       >>>
       >>> eg = gt.get_block_edge_gradient(g, state.levels[0].get_edge_blocks())
       >>> gt.graph_draw(g, pos=pos, vertex_pie_fractions=bc,
       ...               vertex_pie_colors=bv, vertex_shape="pie",
       ...               vertex_size=10, vertex_anchor=0, vorder=d,
       ...               edge_gradient=eg, edge_control_points=cts,
       ...               output="power_nested_mdl_overlap.pdf")
       <...>

    .. testcleanup:: nested_mdl_overlap

       gt.graph_draw(g, pos=pos, vertex_pie_fractions=bc, vertex_pie_colors=bv, vertex_shape="pie", edge_gradient=eg, edge_control_points=cts, vertex_size=10, vertex_anchor=0, vorder=d, output="power_nested_mdl_overlap.png")

    .. figure:: power_nested_mdl_overlap.*
       :align: center

       Overlapping block partition of a power-grid network, which minimizes the
       description length of the network according to the nested
       (degree-corrected, overlapping) stochastic blockmodel.



    References
    ----------
    .. [peixoto-hierarchical-2014] Tiago P. Peixoto, "Hierarchical block structures and high-resolution
       model selection in large networks ", Phys. Rev. X 4, 011047 (2014), :doi:`10.1103/PhysRevX.4.011047`,
       :arxiv:`1310.4377`.
    .. [peixoto-efficient-2014] Tiago P. Peixoto, "Efficient Monte Carlo and greedy
       heuristic for the inference of stochastic block models", Phys. Rev. E 89, 012804 (2014),
       :doi:`10.1103/PhysRevE.89.012804`, :arxiv:`1310.4378`.
    .. [peixoto-model-2014] Tiago P. Peixoto, "Model selection and hypothesis
       testing for large-scale network models with overlapping groups",
       :arxiv:`1409.3059`.
    """

    dl_ent = kwargs.get("dl_ent", False)

    if minimize_state is None:
        minimize_state = NestedMinimizeState()

    if overlap and nonoverlap_init and minimize_state.init and bs is None:
        if verbose:
            print("Non-overlapping initialization...")
        state = minimize_nested_blockmodel_dl(g, Bs=Bs, bs=bs,
                                              min_B=min_B,
                                              max_B=max_B,
                                              deg_corr=deg_corr, overlap=False,
                                              dl=dl, dense=dense,
                                              multigraph=multigraph,
                                              eweight=eweight,
                                              vweight=vweight,
                                              clabel=clabel if isinstance(clabel, PropertyMap) else None,
                                              nsweeps=nsweeps,
                                              epsilon=epsilon, c=c,
                                              nmerge_sweeps=nmerge_sweeps, r=r,
                                              sparse_thresh=sparse_thresh,
                                              sequential=sequential,
                                              verbose=verbose,
                                              checkpoint=checkpoint,
                                              minimize_state=minimize_state,
                                              dl_ent=dl_ent)
        if overlap != "full":
            if clabel is not None:
                bstate = state.levels[0].copy(overlap=True, clabel=clabel)
            else:
                bstate = state.levels[0].copy(overlap=True,
                                              clabel=g.new_vertex_property("int"))
            unilevel_minimize(bstate, nsweeps=nsweeps,
                              epsilon=epsilon, c=c,
                              nmerge_sweeps=nmerge_sweeps,
                              dl=dl,
                              sequential=sequential)
            bs = [s.b.a for s in state.levels]
            bs[0] = bstate.b.a

        else:
            bstates = [bstate.copy(overlap=True) for bstate in state.levels]
            bs = [s.b.a for s in bstates]
        if nonoverlap_init != "full":
            bs = [bs[0], zeros(bs[0].max() + 1, dtype=bs[0].dtype)]

        Bs = [b.max() + 1 for b in bs]
        max_B = Bs[0]
        max_b = bs[0].copy()

        minimize_state.clear()
        minimize_state.init = False

        if verbose:
            print("Overlapping minimization starting from:")
            state.print_summary()

    if Bs is None:
        if minimize_state is not None:
            Bs = [ba.max() + 1 for ba in minimize_state.bs]
            if len(Bs) == 0:
                Bs = [1]
        else:
            Bs = [1]

    if bs is not None:
        Bs = [ba.max() + 1 for ba in bs]

    if Bs[-1] > 1:
        Bs += [1]

    if bs is None:
        state = init_nested_state(g, Bs=Bs, deg_corr=deg_corr,
                                  overlap=overlap, eweight=eweight,
                                  vweight=vweight, clabel=clabel,
                                  verbose=verbose,
                                  nsweeps=nsweeps, nmerge_sweeps=nmerge_sweeps,
                                  dl=dl, dense=dense, multigraph=multigraph,
                                  epsilon=epsilon, sparse_thresh=sparse_thresh,
                                  checkpoint=checkpoint,
                                  minimize_state=minimize_state,
                                  dl_ent=dl_ent)
    else:
        state = NestedBlockState(g, bs=bs, deg_corr=deg_corr, overlap=overlap,
                                 eweight=eweight, vweight=vweight, clabel=clabel)

    minimize_state.sync(state)

    chkp = get_checkpoint_wrap(checkpoint, state, minimize_state, dl_ent)

    dS = nested_tree_sweep(state,
                           min_B=min_B,
                           max_B=max_B,
                           max_b=max_b,
                           verbose=verbose,
                           nsweeps=nsweeps,
                           nmerge_sweeps=nmerge_sweeps,
                           r=r, epsilon=epsilon,
                           dense=dense, dl=dl,
                           multigraph=multigraph,
                           sparse_thresh=sparse_thresh,
                           checkpoint=chkp,
                           minimize_state=minimize_state,
                           frozen_levels=frozen_levels,
                           dl_ent=dl_ent)

    return state

def get_hierarchy_tree(bstack):
    r"""Obtain the nested hierarchical levels as a tree.

    This transforms a list of :class:`~graph_tool.Graph` instances, given by the
    ``bstack`` parameter, representing the inferred hierarchy at each level, into
    a single :class:`~graph_tool.Graph` containing the hierarchy tree.

    It is expected for each graph in ``bstack`` to contain an internal vertex
    property map named "b" which contains the node partition at the
    corresponding level.

    Returns
    -------

    tree : :class:`~graph_tool.Graph`
       A directed graph, where vertices are blocks, and a directed edge points
       to an upper to a lower level in the hierarchy.
    label : :class:`~graph_tool.PropertyMap`
       A vertex property map containing the block label for each node.
    """

    g = bstack[0]
    b = g.vp["b"]
    bstack = bstack[1:]

    t = Graph()

    t.add_vertex(g.num_vertices())
    label = t.vertex_index.copy("int")

    last_pos = 0
    for l, s in enumerate(bstack):
        pos = t.num_vertices()
        t.add_vertex(s.num_vertices())
        label.a[-s.num_vertices():] = arange(s.num_vertices())

        for v in g.vertices():
            w = t.vertex(int(v) + last_pos)
            u = t.vertex(b[v] + pos)
            t.add_edge(u, w)

        last_pos = pos
        g = s
        if g.num_vertices() == 1:
            break
        b = g.vp["b"]
    return t, label
