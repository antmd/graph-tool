#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# graph_tool -- a general graph manipulation python module
#
# Copyright (C) 2006-2014 Tiago de Paula Peixoto <tiago@skewed.de>
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


class NestedBlockState(object):
    r"""This class encapsulates the nested block state of a given graph.

    This must be instantiated and used by functions such as :func:`nested_mcmc_sweep`.

    The instances of this class contain a data member called ``levels``, which
    is a list of :class:`~graph_tool.community.BlockState` instances, containing
    the entire nested hierarchy.

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
    max_BE : ``int`` (optional, default: ``1000``)
        If the number of blocks exceeds this number, a sparse representation of
        the block graph is used, which is slightly less efficient, but uses less
        memory,
    """

    def __init__(self, g, eweight=None, vweight=None, bs=None,
                 Bs=None, deg_corr=True, max_BE=1000):
        L = len(Bs) if Bs is not None else len(bs)
        cg = g
        vcount = vweight
        ecount = eweight

        self.levels = []

        for i in range(L):
            Bl = Bs[i] if Bs is not None else None
            bl = None
            if bs is not None:
                if isinstance(bs[i], PropertyMap):
                    bl = cg.own_property(bs[i])
                else:
                    bl = cg.new_vertex_property("int")
                    bl.fa = bs[i]
            state = BlockState(cg, B=Bl, b=bl, eweight=ecount,
                               vweight=vcount,
                               deg_corr=deg_corr and i == 0,
                               max_BE=max_BE)
            self.levels.append(state)

            cg = state.bg
            vcount = None # no vweight for upper levels!
            ecount = state.mrs

    def __rebuild_level(self, l, bl=None, bmap=None):
        if bl is None:
            if l > 0:
                cg = self.levels[l - 1].bg
                ecount = self.levels[l - 1].mrs
                if bmap is not None: # translate b values
                    bl = cg.new_vertex_property("int")
                    bl.fa = bmap
                else:
                    bl = self.levels[l].b
                bl = cg.own_property(bl)
                state = BlockState(cg, b=bl, eweight=ecount,
                                   deg_corr=False, max_BE=self.levels[0].max_BE)
                if l < len(self.levels):
                    self.levels[l] = state
                else:
                    self.levels.append(state)
        else:
            cg = self.levels[l].g
            ecount = self.levels[l].eweight
            vcount = self.levels[l].vweight
            if isinstance(bl, PropertyMap):
                bl = cg.own_property(bl)
            else:
                b = cg.new_vertex_property("int")
                b.fa = bl
                bl = b
            state = BlockState(cg, b=bl, eweight=ecount, vweight=vcount,
                               deg_corr=self.levels[l].deg_corr,
                               max_BE=self.levels[0].max_BE)
            self.levels[l] = state
            if l == len(self.levels) - 1 and state.B > 1:
                self._NestedBlockState__rebuild_level(l + 1, bmap=zeros(state.B, dtype="int"))

    def __delete_level(self, l):
        if l == 0:
            raise ValueError("cannot delete level l=0")
        clabel = self.levels[l].b.copy()
        b = self.levels[l - 1].b.copy()
        assert(self.levels[l - 1].B == self.levels[l].N)
        libcommunity.vector_map(b.a, clabel.a)
        del self.levels[l]
        self._NestedBlockState__rebuild_level(l - 1, bl=b)
        assert(self.levels[l - 1].B == self.levels[l].N)

    def __duplicate_level(self, l):
        cg = self.levels[l].bg
        ecount = self.levels[l].mrs

        bstate = BlockState(cg,
                            B=cg.num_vertices(),
                            #b=cg.vertex_index.copy("int"),
                            eweight=ecount,
                            deg_corr=False,
                            max_BE=self.levels[l].max_BE)

        self.levels.insert(l + 1, bstate)
        assert(self.levels[l].B == self.levels[l + 1].N)


    def level_entropy(self, l, complete=False, random=False, dense=False,
                      multigraph=False):
        r"""Compute the description length of hierarchy level l.

        Parameters
        ----------
        l : ``int``
            Hierarchy level.
        complete : ``bool`` (optional, default: ``False``)
            If ``True``, the complete entropy will be returned, including constant
            terms not relevant to the block partition.
        random : ``bool`` (optional, default: ``False``)
            If ``True``, the entropy entropy corresponding to an equivalent random
            graph (i.e. no block partition) will be returned.
        dense : ``bool`` (optional, default: ``False``)
            If ``True``, the "dense" variant of the entropy will be computed.
        multigraph : ``bool`` (optional, default: ``False``)
            If ``True``, the multigraph entropy will be used. Only has an effect
            if ``dense == True``.
        """

        bstate = self.levels[l]
        B = bstate.B
        N = bstate.N
        nr = bstate.wr.a

        S = bstate.entropy(dl=False, dense=dense or l > 0,
                           multigraph=multigraph or l > 0,
                           complete=complete) + \
            partition_entropy(B=B, N=N, nr=nr) / bstate.E

        if l == 0 and complete:
            if bstate.deg_corr:
                S_seq = libcommunity.deg_entropy(bstate.g._Graph__graph,
                                                 _prop("v", bstate.g,
                                                       bstate.b),
                                                 bstate.B)
                S += S_seq / bstate.E

        return S

    def entropy(self, complete=False, random=False, dense=False,
                multigraph=False):
        r"""Compute the description length of the entire hierarchy.

        Parameters
        ----------
        complete : ``bool`` (optional, default: ``False``)
            If ``True``, the complete entropy will be returned, including constant
            terms not relevant to the block partition.
        random : ``bool`` (optional, default: ``False``)
            If ``True``, the entropy entropy corresponding to an equivalent random
            graph (i.e. no block partition) will be returned.
        dense : ``bool`` (optional, default: ``False``)
            If ``True``, the "dense" variant of the entropy will be computed.
        multigraph : ``bool`` (optional, default: ``False``)
            If ``True``, the multigraph entropy will be used. Only has an effect
            if ``dense == True``.
        """

        S = 0
        for l in range(len(self.levels)):
            S += self.level_entropy(l, complete, random, dense, multigraph)
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
        r"""Project the partition at level ``l`` onto the lowest level, and
        return the corresponding property map.
        """
        clabel = b = self.levels[l].b.copy()
        while l - 1 >= 0:
            clabel = b
            b = self.levels[l - 1].b.copy()
            libcommunity.vector_map(b.a, clabel.a)
            l -= 1
        return b


def nested_mcmc_sweep(state, beta=1., random_move=False, c=1., sequential=True,
                      verbose=False):
    r"""Performs a Markov chain Monte Carlo sweep on all levels of the hierarchy.

    Parameters
    ----------
    state : :class:`~graph_tool.community.NestedBlockState`
        The nested block state.
    beta : `float` (optional, default: `1.0`)
        The inverse temperature parameter :math:`\beta`.
    random_move : ``bool`` (optional, default: ``False``)
        If ``True``, the proposed moves will attempt to place the vertices in
        fully randomly-chosen blocks. If ``False``, the proposed moves will be
        chosen with a probability depending on the membership of the neighbours
        and the currently-inferred block structure.
    c : ``float`` (optional, default: ``1.0``)
        This parameter specifies how often fully random moves are attempted,
        instead of more likely moves based on the inferred block partition.
        For ``c == 0``, no fully random moves are attempted, and for ``c == inf``
        they are always attempted.
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
       [(-0.11881394061738723, 58), (0.0, 0), (0.0, 0), (-0.00108916046910437, 1), (0.0, 0)]

    References
    ----------

    .. [peixoto-efficient-2013] Tiago P. Peixoto, "Efficient Monte Carlo and greedy
       heuristic for the inference of stochastic block models", :arxiv:`1310.4378`.
    .. [peixoto-hierarchical-2013] Tiago P. Peixoto, "Hierarchical block structures and high-resolution
       model selection in large networks ", :arxiv:`1310.4377`.
    """

    rets = []
    for l, bstate in enumerate(state.levels):
        if verbose:
            print("Level:", l, "N:", bstate.N, "B:", bstate.B)
        ret = mcmc_sweep(bstate, beta=beta, c=c, dense = l > 0,
                         multigraph = l > 0,
                         sequential=sequential,
                         random_move=random_move, verbose=verbose)
        if l + 1 < len(state.levels):
            state._NestedBlockState__rebuild_level(l + 1)
        rets.append(ret)
    return rets


def replace_level(l, state, nsweeps=10, nmerge_sweeps=10, r=2,
                  c=0, epsilon=0., sequential=True, random_move=False,
                  dense=False, sparse_thresh=100, verbose=False,
                  checkpoint=None, minimize_state=None):
    bstate = state.levels[l]
    g = bstate.g

    if l + 1 < len(state.levels):
        min_B = state.levels[l + 1].B
    else:
        min_B = 1
    max_B = bstate.N

    if l + 1 < len(state.levels):
        clabel = bstate.b.copy()
        b_upper = state.levels[l + 1].b
        libcommunity.vector_map(clabel.a, b_upper.a)
    else:
        clabel = None

    nested_dl = l < len(state.levels) - 1

    Sb = get_state_dl(state.levels[l], dense=l > 0 or dense,
                      nested_dl=nested_dl, clabel=clabel)

    res = minimize_blockmodel_dl(g, eweight=bstate.eweight,
                                 deg_corr=bstate.deg_corr, nsweeps=nsweeps,
                                 c=c, random_move=random_move,
                                 sequential=sequential,
                                 r=r,
                                 adaptive_sweeps=True, greedy_cooling=True,
                                 epsilon=epsilon, nmerge_sweeps=nmerge_sweeps,
                                 max_B=max_B, min_B=min_B,
                                 clabel=clabel,
                                 max_BE=bstate.max_BE,
                                 dense=l > 0 or dense,
                                 sparse_heuristic=g.num_vertices() > sparse_thresh,
                                 nested_dl=nested_dl,
                                 verbose=False,
                                 ##exaustive=g.num_vertices() <= 100,
                                 minimize_state=minimize_state.minimize_state,
                                 checkpoint=checkpoint)
    b = res[0]

    # rebuild current level
    state._NestedBlockState__rebuild_level(l, bl=b.copy())

    Sf = res[1]

    kept = False
    if Sf - Sb >= -1e-10:
        kept = True
    if Sf - Sb == 0 and bstate.B != state.levels[l].B:
        kept = False

    if kept:
        state.levels[l] = bstate
        Sf = Sb
    else:
        # rebuild upper level
        if l + 1 < len(state.levels):
            if clabel is None:
                clabel = state.levels[l].b.copy()
                b_upper = state.levels[l + 1].b
                libcommunity.vector_map(clabel.a, b_upper.a)
            bmap = arange(b.fa.max() + 1, dtype=b.fa.dtype)
            libcommunity.vector_rmap(b.fa, bmap)
            libcommunity.vector_map(bmap, clabel.fa)

            state._NestedBlockState__rebuild_level(l + 1, bmap=bmap)

    if verbose:
        print("level", l, ": resizing", bstate.B, "->", state.levels[l].B, ", dS:", Sf - Sb,
              "[kept]" if kept else "" )

    dS = Sf - Sb
    return dS, kept


class NestedMinimizeState(object):
    r"""This object stores information regarding the current entropy minimization
    state, so that the algorithms can resume previously started runs.
    This object can be saved to disk via the :mod:`pickle` interface."""

    def __init__(self):
        self.minimize_state = MinimizeState()
        self.l = -1
        self.bs = []
        self.done = []

    def clear(self):
        self.minimize_state.clear()
        self.l = -1
        self.bs.clear()
        self.done.clear()

    def sync(self, state):
        if len(self.bs) == 0:
            for s in state.levels:
                self.bs.append(array(s.b.a))
        while len(self.done) < len(state.levels):
            self.done.append(False)
        if self.l == -1:
            self.l = len(state.levels) - 1

    def delete(self, l):
        del self.done[l]
        del self.bs[l]

    def insert(self, l, state):
        self.done.insert(l + 1, False)
        ba = state.levels[l].b.fa
        self.bs.insert(l + 1, array(ba))

    def mark_level(self, l, done, state):
        while len(state.levels) > len(self.bs):
            ba = state.levels[len(self.bs)].b.fa
            self.bs.append(array(ba))
            self.done.append(False)

        self.done[l] = done
        if done:
            self.bs[l] = array(state.levels[l].b.fa)

    def clear_mstate(self):
        self.minimize_state = MinimizeState()

def nested_tree_sweep(state, nsweeps=10, epsilon=0., r=2.,
                      nmerge_sweeps=10, random_move=False,
                      c=0, dense=False, sequential=True,
                      sparse_thresh=100, checkpoint=None,
                      minimize_state=None, verbose=False):
    r"""Performs one greedy sweep in the entire hierarchy tree, attempting to
    decrease its description length.

    Parameters
    ----------
    state : :class:`~graph_tool.community.NestedBlockState`
        The nested block state.
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
    random_move : ``bool`` (optional, default: ``False``)
        If ``True``, the proposed moves will attempt to place the vertices in
        fully randomly-chosen blocks. If ``False``, the proposed moves will be
        chosen with a probability depending on the membership of the neighbours
        and the currently-inferred block structure.
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

    .. [peixoto-efficient-2013] Tiago P. Peixoto, "Efficient Monte Carlo and greedy
       heuristic for the inference of stochastic block models", :arxiv:`1310.4378`.
    .. [peixoto-hierarchical-2013] Tiago P. Peixoto, "Hierarchical block structures
       and high-resolution model selection in large networks ", :arxiv:`1310.4377`.
    """

    if minimize_state is None:
        minimize_state = NestedMinimizeState()
    mstate = minimize_state
    mstate.sync(state)

    args = dict(state=state, nsweeps=nsweeps, nmerge_sweeps=nmerge_sweeps, r=r,
                c=c, epsilon=epsilon, sequential=sequential, random_move=random_move,
                dense=dense, sparse_thresh=sparse_thresh,
                checkpoint=checkpoint, minimize_state=minimize_state)

    #_Si = state.entropy(dense=dense, multigraph=dense)
    dS = 0

    while mstate.l >= 0:
        l = mstate.l
        if mstate.done[l]:
            if verbose:
                print("level", l, ": skipping", state.levels[l].B)


            mstate.l -= 1
            continue

        ddS, kept = replace_level(l, verbose=verbose, **args)
        dS += ddS
        mstate.clear_mstate()

        if l > 0:
            assert(state.levels[l - 1].B == state.levels[l].N)

        if l + 1 < len(state.levels):
            assert(state.levels[l].B == state.levels[l + 1].N)

        # delete level
        if kept and l > 0 and l < len(state.levels) - 1:
            Si = state.entropy()

            bstates = [state.levels[l-1], state.levels[l]]

            state._NestedBlockState__delete_level(l)
            replace_level(l, **args)

            assert(state.levels[l-1].B == state.levels[l].N)
            if l + 1 < len(state.levels):
                assert(state.levels[l].B == state.levels[l + 1].N)

            Sf = state.entropy()

            mstate.clear_mstate()
            if Sf > Si:
                state.levels[l - 1] = bstates[0]
                state.levels.insert(l, bstates[1])
            else:
                kept = False
                dS += Sf - Si

                mstate.delete(l)

                if verbose:
                    print("level", l, ": deleted", bstates[1].B, ", dS:", Sf - Si, len(state.levels))

        # insert new level
        if kept:
            Si = state.entropy()

            bstate = state.levels[l]

            state._NestedBlockState__duplicate_level(l)

            replace_level(l, verbose=False, **args)

            assert(state.levels[l].B == state.levels[l+1].N)
            if l > 0:
                assert(state.levels[l - 1].B == state.levels[l].N)

            Sf = state.entropy()

            mstate.clear_mstate()
            if Sf >= Si:
                del state.levels[l + 1]
                state.levels[l] = bstate
            else:
                kept = False
                dS += Sf - Si

                mstate.insert(l, state)

                l += 1

                if verbose:
                    print("level", l, ": inserted", state.levels[l].B, ", dS:", Sf - Si)


        if not kept:
            mstate.mark_level(l, done=True, state=state)
            if l + 1 < len(state.levels):
                mstate.mark_level(l + 1, done=False, state=state)
            if l > 0:
                mstate.mark_level(l - 1, done=False, state=state)

            l += 1
        else:
            l -= 1

        if l >= len(state.levels):
            l = len(state.levels) - 1

        if l == len(state.levels) - 1 and state.levels[l].B > 1:
            state._NestedBlockState__rebuild_level(l + 1, bmap=zeros(state.levels[l].B))
            mstate.mark_level(l + 1, done=False, state=state)
            l += 1

        mstate.l = l
        if checkpoint is not None:
            checkpoint(None, 0, 0, 0, mstate)

    # _Sf = state.entropy(dense=dense, multigraph=dense)

    return dS



def init_nested_state(g, Bs, deg_corr=True, dense=False,
                      eweight=None, vweight=None,
                      nsweeps=10, epsilon=0., r=2,
                      nmerge_sweeps=10, random_move=False,
                      c=0, sequential=True,
                      sparse_thresh=100, checkpoint=None,
                      minimize_state=None, max_BE=1000,
                      verbose=False):
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
    random_move : ``bool`` (optional, default: ``False``)
        If ``True``, the proposed moves will attempt to place the vertices in
        fully randomly-chosen blocks. If ``False``, the proposed moves will be
        chosen with a probability depending on the membership of the neighbours
        and the currently-inferred block structure.
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

    .. [peixoto-efficient-2013] Tiago P. Peixoto, "Efficient Monte Carlo and greedy
       heuristic for the inference of stochastic block models", :arxiv:`1310.4378`.
    .. [peixoto-hierarchical-2013] Tiago P. Peixoto, "Hierarchical block structures
       and high-resolution model selection in large networks ", :arxiv:`1310.4377`.
    """

    if minimize_state is None:
        minimize_state = NestedMinimizeState()
    mstate = minimize_state

    state = NestedBlockState(g, eweight=eweight, vweight=vweight, Bs=[1],
                             deg_corr=deg_corr)

    if checkpoint is not None:
        def check_wrap(bstate, S, delta, nmoves, ms):
            if bstate is None:
                checkpoint(None, S, delta, nmoves, minimize_state)
            else:
                checkpoint(state, S, delta, nmoves, minimize_state)
        chkp = check_wrap
    else:
        chkp = None

    bg = g
    ecount = eweight

    for l, B in enumerate(Bs):
        if l < len(mstate.bs):
            ba = mstate.bs[l]
        else:
            bstate = BlockState(bg, B=bg.num_vertices(), #b=bg.vertex_index.copy("int"),
                                vweight=vweight if l == 0 else None,
                                eweight=ecount,
                                deg_corr=deg_corr and l == 0, max_BE=max_BE)

            bstate = multilevel_minimize(bstate, B, nsweeps=nsweeps, epsilon=epsilon,
                                         r=r, nmerge_sweeps=nmerge_sweeps,
                                         greedy=True, c=c,
                                         dense=(l > 0 and g.num_vertices() < sparse_thresh) or dense,
                                         multigraph=l > 0 or dense,
                                         random_move=random_move,
                                         sequential=sequential, verbose=verbose,
                                         checkpoint=chkp,
                                         minimize_state=minimize_state.minimize_state)
            ba = array(bstate.b.fa)
            mstate.bs.append(ba)
            minimize_state.minimize_state.clear()

        state._NestedBlockState__rebuild_level(len(state.levels) - 1, bl=ba)
        bg = state.levels[l].bg
        ecount = state.levels[l].mrs

    minimize_state.clear()
    mstate.sync(state)
    if checkpoint is not None:
        checkpoint(None, 0, 0, 0, mstate)

    return state

def minimize_nested_blockmodel_dl(g, Bs=None, bs=None, deg_corr=True,
                                  dense=False, eweight=None, vweight=None,
                                  nsweeps=10, epsilon=0., random_move=False, c=0,
                                  nmerge_sweeps=10, r=2, sparse_thresh=100,
                                  sequential=True, verbose=False,
                                  checkpoint=None, minimize_state=None):
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
        The number of sweeps necessary for the local minimum will
        be estimated to be enough so that no more than ``epsilon * N`` nodes
        changes their states in the last ``nsweeps`` sweeps.
    random_move : ``bool`` (optional, default: ``False``)
        If ``True``, the proposed moves will attempt to place the vertices in
        fully randomly-chosen blocks. If ``False``, the proposed moves will be
        chosen with a probability depending on the membership of the neighbours
        and the currently-inferred block structure.
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
    bstack : list of :class:`~graph_tool.Graph` instances
        This is a list of :class:`~graph_tool.Graph` instances representing the
        inferred hierarchy at each level. Each graph has two internal vertex and
        edge property maps named "count" which correspond to the vertex and edge
        counts at the lower level, respectively. Additionally, an internal
        vertex property map named "b" specifies the block partition.
    min_dl : ``float``
       Minimum value of the description length (in `nats <http://en.wikipedia.org/wiki/Nat_%28information%29>`_ per edge).

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


    See [peixoto-hierarchical-2013]_ for details on the algorithm.

    This algorithm has a complexity of :math:`O(N \ln^2 N)`, where :math:`N`
    is the number of nodes in the network.

    Examples
    --------
    .. testsetup:: nested_mdl

       gt.seed_rng(42)
       np.random.seed(42)

    .. doctest:: nested_mdl

       >>> g = gt.collection.data["power"]
       >>> bstack, mdl = gt.minimize_nested_blockmodel_dl(g, deg_corr=True)
       >>> t = gt.get_hierarchy_tree(bstack)[0]
       >>> tpos = pos = gt.radial_tree_layout(t, t.vertex(t.num_vertices() - 1), weighted=True)
       >>> cts = gt.get_hierarchy_control_points(g, t, tpos)
       >>> pos = g.own_property(tpos)
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

    References
    ----------

    .. [peixoto-hierarchical-2013] Tiago P. Peixoto, "Hierarchical block structures and high-resolution
       model selection in large networks ", :arxiv:`1310.4377`.
    .. [peixoto-efficient-2013] Tiago P. Peixoto, "Efficient Monte Carlo and greedy
       heuristic for the inference of stochastic block models", :arxiv:`1310.4378`.
    """


    if Bs is None:
        Bs = [1]

        if minimize_state is not None:
            Bs = [ba.max() + 1 for ba in minimize_state.bs]

    if bs is not None:
        Bs = [ba.max() + 1 for ba in bs]

    if Bs[-1] > 1:
        Bs += [1]

    if minimize_state is None:
        minimize_state = NestedMinimizeState()

    if bs is None:
        state = init_nested_state(g, Bs=Bs, deg_corr=deg_corr,
                                  eweight=eweight, vweight=vweight,
                                  verbose=verbose, nsweeps=nsweeps,
                                  nmerge_sweeps=nmerge_sweeps,
                                  dense=dense,
                                  epsilon=epsilon,
                                  sparse_thresh=sparse_thresh,
                                  checkpoint=checkpoint,
                                  minimize_state=minimize_state)
    else:
        state = NestedBlockState(g, bs=bs, deg_corr=deg_corr, eweight=eweight,
                                 vweight=vweight)

    if checkpoint is not None:
        def check_wrap(bstate, S, delta, nmoves, ms):
            if bstate is not None:
                checkpoint(state, S, delta, nmoves, minimize_state)
            else:
                checkpoint(None, S, delta, nmoves, minimize_state)
        chkp = check_wrap
    else:
        chkp = None


    dS = nested_tree_sweep(state, verbose=verbose,
                            nsweeps=nsweeps,
                            nmerge_sweeps=nmerge_sweeps,
                            r=r,
                            epsilon=epsilon,
                            sparse_thresh=sparse_thresh,
                            checkpoint=chkp,
                            minimize_state=minimize_state)

    bstack = state.get_bstack()

    return bstack, state.entropy()


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
