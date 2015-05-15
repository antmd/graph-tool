#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# graph_tool -- a general graph manipulation python module
#
# Copyright (C) 2006-2016 Tiago de Paula Peixoto <tiago@skewed.de>
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
from scipy.special import gammaln
import contextlib

from .. import group_vector_property, ungroup_vector_property, Vector_size_t, \
    openmp_get_num_threads, openmp_enabled

from .. decorators import _wraps

from .. dl_import import dl_import
dl_import("from . import libgraph_tool_community as libcommunity")

from .. generation import graph_union
from .. stats import vertex_hist

from . blockmodel import *
from . blockmodel import _bm_test

class CovariateBlockState(BlockState):
    r"""This class encapsulates the (possibly overlapping) block state of a given
    graph, where the edges are divided into discrete layers.

    This must be instantiated and used by functions such as :func:`mcmc_sweep`.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be modelled.
    ec : :class:`~graph_tool.PropertyMap`
        Edge :class:`~graph_tool.PropertyMap` containing edge covariates that
        will split the network in discrete layers.
    eweight : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Edge multiplicities (for multigraphs or block graphs).
    vweight : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Vertex multiplicities (for block graphs).
    b : :class:`~graph_tool.PropertyMap` or :class:`numpy.ndarray` (optional, default: ``None``)
        Initial block labels on the vertices or half-edges. If not supplied, it
        will be randomly sampled.
    B : ``int`` (optional, default: ``None``)
        Number of blocks. If not supplied it will be either obtained from the
        parameter ``b``, or set to the maximum possible value according to the
        minimum description length.
    clabel : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Constraint labels on the vertices. If supplied, vertices with different
        label values will not be clustered in the same group.
    layers : ``bool`` (optional, default: ``False``)
        If ``layers==True``, the "independent layers" version of the model is
        used, instead of the "edge covariates" version.
    deg_corr : ``bool`` (optional, default: ``True``)
        If ``True``, the degree-corrected version of the blockmodel ensemble will
        be assumed, otherwise the traditional variant will be used.
    overlap : ``bool`` (optional, default: ``False``)
        If ``True``, the overlapping version of the moddel will be used.
    """

    def __init__(self, g, ec, eweight=None, vweight=None, b=None, B=None,
                 clabel=None, layers=False, deg_corr=True, overlap=False,
                 **kwargs):

        BlockState._state_ref_count += 1

        self.g = g

        if kwargs.get("ec_done", False):
            self.ec = ec
        else:
            self.ec = g.new_edge_property("int")
            libcommunity.ec_hist(g._Graph__graph, _prop("e", g, ec),
                                 _prop("e", g, self.ec))
            ec = self.ec

        self.C = ec.fa.max() + 1
        self.layers = layers

        if "max_BE" in kwargs:
            del kwargs["max_BE"]
        max_BE = 0

        if not overlap:
            total_state = BlockState(GraphView(g, skip_properties=True),
                                     b=b, B=B, eweight=eweight,
                                     vweight=vweight, clabel=clabel,
                                     deg_corr=deg_corr,
                                     max_BE=max_BE, **kwargs)
        else:
            total_state = OverlapBlockState(g, b=b, B=B, eweight=eweight,
                                            vweight=vweight, clabel=clabel,
                                            deg_corr=deg_corr,
                                            max_BE=max_BE, **kwargs)
            self.base_g = total_state.base_g
            self.g = total_state.g

        self.total_state = total_state

        if overlap:
            self.base_ec = ec.copy()
            ec = total_state.eindex.copy()
            pmap(ec, self.ec)
            self.ec = ec.copy("int")

        self.eweight = total_state.eweight
        self.vweight = total_state.vweight

        self.E = total_state.E
        self.N = total_state.N

        self.is_weighted = total_state.is_weighted

        self.b = total_state.b
        self.B = total_state.B
        self.clabel = total_state.clabel

        self.deg_corr = deg_corr
        self.overlap = overlap

        self.vc = group_vector_property([self.g.new_vertex_property("int", -1)])
        self.vmap = group_vector_property([self.g.vertex_index.copy("int")])

        self.gs = []
        self.bmap = libcommunity.bmap_t()

        for l in range(0, self.C):
            u = Graph(directed=g.is_directed())
            u.vp["b"] = u.new_vertex_property("int")
            u.vp["weight"] = u.new_vertex_property("int")
            u.ep["weight"] = u.new_edge_property("int")
            u.vp["brmap"] = u.new_vertex_property("int")
            u.vp["vmap"] = u.new_vertex_property("int")
            self.gs.append(u)

        libcommunity.split_graph(self.g._Graph__graph,
                                 _prop("e", self.g, self.ec),
                                 _prop("v", self.g, self.b),
                                 _prop("e", self.g, self.eweight),
                                 _prop("v", self.g, self.vweight),
                                 _prop("v", self.g, self.vc),
                                 _prop("v", self.g, self.vmap),
                                 [u._Graph__graph for u in self.gs],
                                 [_prop("v", u, u.vp["b"]) for u in self.gs],
                                 [_prop("e", u, u.ep["weight"]) for u in self.gs],
                                 [_prop("v", u, u.vp["weight"]) for u in self.gs],
                                 self.bmap,
                                 [_prop("v", u, u.vp["brmap"]) for u in self.gs],
                                 [_prop("v", u, u.vp["vmap"]) for u in self.gs])

        if not self.layers:
            total_state.master = True
            total_state.slave = False
        else:
            total_state.master = False
            total_state.slave = False
        total_state.free_blocks = Vector_size_t()
        total_state.g.vp["brmap"] = total_state.g.vertex_index.copy("int")
        total_state.g.vp["vmap"] = total_state.g.vertex_index.copy("int")

        self.states = [total_state]

        self.max_BE = max_BE
        for u in self.gs:
            state = self.__gen_state(u)
            self.states.append(state)
            if _bm_test():
                assert state.mrs.fa.sum() == state.eweight.fa.sum(), ("inconsistent mrs!", l)


        self.__bg = None
        self.__mrs = None
        self.__bec = None
        self.__dummy_bg = Graph(directed=g.is_directed())
        self.__dummy_bg.add_vertex(self.B)
        #self.wr = self.__dummy_bg.own_property(total_state.wr)
        self.wr = total_state.wr

        self.sweep_vertices = total_state.sweep_vertices
        self.emat = None

        self.overlap_stats = self.total_state.overlap_stats

        self.__layer_entropy = None

        if _bm_test():
            assert self.mrs.fa.sum() == self.eweight.fa.sum(), "inconsistent mrs!"

        # computation cache
        libcommunity.init_safelog(int(5 * max(self.E, self.N)))
        libcommunity.init_xlogx(int(5 * max(self.E, self.N)))
        libcommunity.init_lgamma(int(3 * max(self.E, self.N)))


    def __get_base_u(self, u):
        node_index = u.vp["vmap"].copy("int64_t")
        pmap(node_index, self.total_state.node_index)
        base_u, nindex, vcount, ecount = \
            condensation_graph(u, node_index,
                               self_loops=True,
                               parallel_edges=True)[:4]
        rindex = zeros(nindex.a.max() + 1, dtype="int64")
        reverse_map(nindex, rindex)
        pmap(node_index, rindex)
        base_u.vp["vmap"] = nindex
        return base_u, node_index

    def __gen_state(self, u):
        nt = 1
        if openmp_enabled():
            nt = openmp_get_num_threads()
        B = u.num_vertices() + 2 * nt
        #B = max(B, u.vp["b"].a.max() + 1 + 2 * nt)
        if not self.overlap:
            state = BlockState(u, b=u.vp["b"],
                               B=B,
                               eweight=u.ep["weight"],
                               vweight=u.vp["weight"],
                               deg_corr=self.deg_corr,
                               force_weighted=self.is_weighted,
                               max_BE=self.max_BE)
        else:
            base_u, node_index = self.__get_base_u(u)
            state = OverlapBlockState(u, b=u.vp["b"].a,
                                      B=B,
                                      vweight=u.vp["weight"],
                                      node_index=node_index,
                                      base_g=base_u,
                                      deg_corr=self.deg_corr,
                                      max_BE=self.max_BE)
        state.master = False
        if not self.layers:
            state.slave = True
        else:
            state.slave = False
        state.free_blocks = Vector_size_t()
        return state

    def __merge_decorator(func):
        @contextlib.contextmanager
        def context_wrapper(self, l_src, l_tgt, revert=False):
            return func(self, l_src, l_tgt, revert=revert)

        @_wraps(func)
        def wrapper(self, l_src, l_tgt, revert=False):
            gen = context_wrapper(self, l_src, l_tgt, revert=revert)
            if revert:
                return gen
            else:
                with gen:
                    pass
        return wrapper

    @__merge_decorator
    def merge_layers(self, l_src, l_tgt, revert=False):
        r"""Merge layer `l_src` into `l_tgt` and delete `l_src`.

        If ``revert == True``, this will return a context that can be used with
        the ``with`` statement, so that when the context is closed, the layers
        are returned to their original state, i.e.

        .. code::

            with state.merge_layers(3, 10, revert=True):
                # here layers 3 and 10 are merged
                pass
            # here layers 3 and 10 exist separately as they did before

        """
        if l_src == l_tgt:
            return
        u_src = self.gs[l_src]
        u_tgt = self.gs[l_tgt]
        s_src = self.states[l_src + 1]
        s_tgt = self.states[l_tgt + 1]

        if self.overlap:
            u_src_base = self.__get_base_u(u_src)[0]
            u_tgt_base = self.__get_base_u(u_tgt)[0]
        else:
            u_src_base = u_src
            u_tgt_base = u_tgt

        intersection = u_src_base.new_vertex_property("int64_t", -1)

        u_tgt_vmap = u_tgt_base.vp["vmap"]
        vmap = {}
        for v in u_tgt_base.vertices():
            vmap[u_tgt_vmap[v]] = v

        u_src_vmap = u_src_base.vp["vmap"]
        for v in u_src_base.vertices():
            w = u_src_vmap[v]
            if w in vmap:
                intersection[v] = int(vmap[w])

        if self.overlap:
            u_tgt_base.ep["b"] = self.states[l_tgt + 1].get_edge_blocks()
            u_src_base.ep["b"] = self.states[l_src + 1].get_edge_blocks()
        else:
            u_tgt_base.vp["b"] = self.states[l_tgt + 1].b
            u_src_base.vp["b"] = self.states[l_src + 1].b

        tgt_bmap = {}
        src_rbmap = {}
        r_max = 0
        for r in range(self.B):
            if self.bmap.has(l_tgt + 1, r):
                tgt_bmap[r] = self.bmap.get(l_tgt + 1, r)
                r_max = max(r_max, tgt_bmap[r])
            if self.bmap.has(l_src + 1, r):
                src_rbmap[self.bmap.get(l_src + 1, r)] = r

        r_missing = list(set(range(r_max)) - set(tgt_bmap.values()))
        r_max += 1

        if self.overlap:
            b = u_src_base.ep["b"].copy()
            for e in u_src_base.edges():
                nb = []
                for r in b[e]:
                    nb.append(src_rbmap[r])
                for i, r in enumerate(nb):
                    if r in tgt_bmap:
                        nb[i] = tgt_bmap[r]
                    else:
                        if len(r_missing) > 0:
                            rr = r_missing[0]
                            del r_missing[0]
                        else:
                            rr = r_max
                            r_max += 1
                        self.bmap.set(l_tgt + 1, r, rr)
                        nb[i] = rr
                        tgt_bmap[r] = rr
                b[e] = nb
            b_src = b
            b_tgt = u_tgt_base.ep["b"]
            u_tgt_base.ep["weight"] = u_tgt_base.new_edge_property("int", 1)
            u_tgt_base.vp["weight"] = u_tgt_base.new_vertex_property("int", 1)
            u_src_base.ep["weight"] = u_src_base.new_edge_property("int", 1)
            u_src_base.vp["weight"] = u_src_base.new_vertex_property("int", 1)
        else:
            b = u_src_base.vp["b"].copy()
            for v in u_src_base.vertices():
                r = src_rbmap[b[v]]
                if r in tgt_bmap:
                    b[v] = tgt_bmap[r]
                else:
                    if len(r_missing) > 0:
                        rr = r_missing[0]
                        del r_missing[0]
                    else:
                        rr = r_max
                        r_max += 1
                    self.bmap.set(l_tgt + 1, r, rr)
                    b[v] = rr
                    tgt_bmap[r] = rr
            b_src = b
            b_tgt = u_tgt_base.vp["b"]

        props = [(b_tgt, b_src),
                 (u_tgt_base.vp["vmap"], u_src_base.vp["vmap"]),
                 (u_tgt_base.vp["weight"], u_src_base.vp["weight"]),
                 (u_tgt_base.ep["weight"], u_src_base.ep["weight"])]

        if not self.overlap:
            props.append((u_tgt_base.vp["brmap"],
                          u_src_base.vp["brmap"]))

        u, props = graph_union(u_tgt_base, u_src_base,
                               intersection=intersection,
                               props=props,
                               include=False)

        if self.overlap:
            u.ep["b"] = props[0]
        else:
            u.vp["b"] = props[0]
            u.vp["brmap"] = props[4]

        u.vp["vmap"] = props[1]
        u.vp["weight"] = props[2]
        u.ep["weight"] = props[3]

        if self.overlap:
            u, b, node_index, half_edges, eindex = half_edge_graph(u, u.ep["b"],
                                                                   self.B)
            u.vp["vmap"] = node_index
            u.vp["weight"] = u.new_vertex_property("int", 1)
            u.vp["b"] = b
            self.gs[l_tgt] = u
            self.states[l_tgt + 1] = self.__gen_state(self.gs[l_tgt])
        else:
            self.gs[l_tgt] = u
            self.states[l_tgt + 1] = self.__gen_state(self.gs[l_tgt])

        del self.states[l_src + 1]
        del self.gs[l_src]

        old_ec = self.ec.copy()
        self.ec.a[self.ec.a == l_src] = l_tgt
        self.ec.a[self.ec.a > l_src] -= 1
        if self.overlap:
            old_base_ec = self.base_ec.copy()
            self.base_ec.a[self.base_ec.a == l_src] = l_tgt
            self.base_ec.a[self.base_ec.a > l_src] -= 1
        self.C -= 1
        old_bmap = self.bmap.copy()
        self.bmap.del_c(l_src + 1)
        self.__bg = None
        old_layer_entropy = self.__layer_entropy
        self.__layer_entropy = None

        yield

        if revert:
            self.gs.insert(l_src, u_src)
            self.gs[l_tgt] = u_tgt
            self.states.insert(l_src + 1, s_src)
            self.states[l_tgt + 1] = s_tgt
            self.ec.a[:] = old_ec.a
            if self.overlap:
                self.base_ec.a[:] = old_base_ec.a
            self.C += 1
            self.bmap = old_bmap
            self.__layer_entropy = old_layer_entropy

    def __getstate__(self):
        state = dict(g=self.g,
                     ec=self.ec,
                     layers=self.layers,
                     eweight=self.eweight,
                     vweight=self.vweight,
                     b=self.b,
                     B=self.B,
                     clabel=self.clabel,
                     deg_corr=self.deg_corr)
        return state

    def __setstate__(self, state):
        self.__init__(**state)
        return state

    def copy(self, b=None, B=None, deg_corr=None, clabel=None, overlap=None,
             layers=None, ec=None):
        r"""Copies the block state. The parameters override the state properties, and
         have the same meaning as in the constructor."""
        state = CovariateBlockState(self.g,
                                    ec=self.ec if ec is None else ec,
                                    eweight=self.eweight,
                                    vweight=self.vweight,
                                    b=self.b if b is None else b,
                                    B=(self.B if b is None else None) if B is None else B,
                                    clabel=(self.clabel.a if self.overlap else self.clabel) if clabel is None else clabel,
                                    deg_corr=self.deg_corr if deg_corr is None else deg_corr,
                                    overlap=self.overlap if overlap is None else overlap,
                                    layers=self.layers if layers is None else layers,
                                    base_g=self.base_g if self.overlap else None,
                                    half_edges=self.total_state.half_edges if self.overlap else None,
                                    node_index=self.total_state.node_index if self.overlap else None,
                                    eindex=self.total_state.eindex if self.overlap else None,
                                    ec_done=ec is None)

        if not state._BlockState__check_clabel():
            b = state.b.a + state.clabel.a * state.B
            continuous_map(b)
            state = state.copy(b=b)

            if _bm_test():
                assert state._BlockState__check_clabel()

        return state

    def __repr__(self):
        return "<CovariateBlockState object with %d %sblocks, %d %s,%s for graph %s, at 0x%x>" % \
            (self.B, "overlapping " if self.overlap else "",
             self.C, "layers" if self.layers else "edge covariates",
             " degree corrected," if self.deg_corr else "",
             str(self.base_g if self.overlap else self.g), id(self))

    def get_bclabel(self):
        r"""Returns a :class:`~graph_tool.PropertyMap`` corresponding to constraint
        labels for the block graph."""
        return self.total_state.get_bclabel()

    def get_bg(self):
        r"""Returns the block graph."""
        if self.__bg is not None:
            return self.__bg, self.__mrs, self.__bec
        bg = Graph(directed=self.g.is_directed())
        mrs = bg.new_edge_property("int")
        ec = bg.new_edge_property("int")

        for l in range(self.C):
            u = GraphView(self.g, efilt=self.ec.a == l)
            ug = get_block_graph(u, self.B, self.b, self.vweight, self.eweight)
            uec = ug.new_edge_property("int")
            uec.a = l
            bg, props = graph_union(bg, ug,
                                    props=[(mrs, ug.ep["count"]),
                                           (ec, uec)],
                                    intersection=ug.vertex_index,
                                    include=True)
            mrs = props[0]
            ec = props[1]

        self.__bg = bg
        self.__mrs = mrs
        self.__bec = ec
        return bg, mrs, ec

    bg = property(lambda self: self.get_bg()[0])
    mrs = property(lambda self: self.get_bg()[1])
    bec = property(lambda self: self.get_bg()[2])

    def get_block_state(self, b=None, vweight=False, deg_corr=False,
                        overlap=False, layers=True):
        r"""Returns a :class:`~graph_tool.community.CovariateBlockState`` corresponding
        to the block graph. The parameters have the same meaning as the in the
        constructor."""

        bg, mrs, ec = self.get_bg()
        state = CovariateBlockState(bg, ec, eweight=mrs,
                                    vweight=bg.own_property(self.wr.copy()) if vweight else None,
                                    b=bg.vertex_index.copy("int") if b is None else b,
                                    clabel=self.get_bclabel(),
                                    deg_corr=deg_corr,
                                    overlap=overlap,
                                    max_BE=self.max_BE,
                                    layers=self.layers if layers is None else layers,
                                    ec_done=True)
        n_map = self.b.copy()
        return state, n_map


    def get_edge_blocks(self):
        r"""Returns an edge property map which contains the block labels pairs for each
        edge."""
        if not self.overlap:
            raise ValueError("edge blocks only available if overlap == True")
        return self.total_state.get_edge_blocks()

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
        if not self.overlap:
            raise ValueError("overlap blocks only available if overlap == True")
        return self.total_state.get_overlap_blocks()

    def get_nonoverlap_blocks(self):
        r"""Returns a scalar-valued vertex property map with the block mixture
        represented as a single number."""

        if not self.overlap:
            return self.b.copy()
        else:
            return self.total_state.get_nonoverlap_blocks()

    def get_majority_blocks(self):
        r"""Returns a scalar-valued vertex property map with the majority block
        membership of each node."""

        if not self.overlap:
            return self.b.copy()
        else:
            return self.total_state.get_majority_blocks()

    def __get_layer_entropy(self):
        if self.__layer_entropy is None:
            if self.layers:
                # we need to include the membership of the nodes in each layer
                g = self.base_g if self.overlap else self.g
                ec = self.base_ec if self.overlap else self.ec
                be = group_vector_property([ec, ec])
                lstate = OverlapBlockState(g, b=be, deg_corr=False)
                self.__layer_entropy = lstate.entropy(dl=True, edges_dl=False) - lstate.entropy(dl=False)
            else:
                self.__layer_entropy = 0
        return self.__layer_entropy

    def entropy(self, complete=True, dl=False, partition_dl=True, edges_dl=True,
                degree_dl=True, dense=False, multigraph=True, norm=False,
                dl_ent=False, **kwargs):
        r"""Calculate the entropy associated with the current block partition.

        Parameters
        ----------
        complete : ``bool`` (optional, default: ``False``)
            If ``True``, the complete entropy will be returned, including constant
            terms not relevant to the block partition.
        dl : ``bool`` (optional, default: ``False``)
            If ``True``, the full description length will be returned.
        partition_dl : ``bool`` (optional, default: ``True``)
            If ``True``, and ``dl == True`` the partition description length
            will be considered.
        edges_dl : ``bool`` (optional, default: ``True``)
            If ``True``, and ``dl == True`` the edge matrix description length
            will be considered.
        degree_dl : ``bool`` (optional, default: ``True``)
            If ``True``, and ``dl == True`` the degree sequence description
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
        """

        S = 0
        if self.layers:
            if dl and partition_dl:
                S += self.total_state.entropy(dl=True, partition_dl=True,
                                              degree_dl=False, edges_dl=False,
                                              norm=False) - \
                     self.total_state.entropy(dl=False, norm=False)
            for state in self.states[1:]:
                S += state.entropy(complete=complete, dl=dl, partition_dl=False,
                                   edges_dl=False, degree_dl=degree_dl,
                                   dense=dense, multigraph=multigraph,
                                   norm=False, dl_ent=dl_ent, **kwargs)
        else:
            for l, state in enumerate(self.states):
                if l == 0:
                    S += state.entropy(complete=complete, dl=dl,
                                       partition_dl=partition_dl,
                                       edges_dl=False, degree_dl=degree_dl,
                                       dense=dense, multigraph=multigraph,
                                       norm=False, dl_ent=dl_ent, **kwargs)
                    S -= libcommunity.covariate_entropy(state.bg._Graph__graph,
                                                        _prop("e", state.bg,
                                                              state.mrs))
                    if multigraph:
                        S -= libcommunity.entropy_parallel(state.g._Graph__graph,
                                                           _prop("e", state.g,
                                                                 state.eweight))
                else:
                    S += libcommunity.covariate_entropy(state.bg._Graph__graph,
                                                        _prop("e", state.bg,
                                                              state.mrs))
                    if multigraph:
                        S += libcommunity.entropy_parallel(state.g._Graph__graph,
                                                           _prop("e", state.g,
                                                                 state.eweight))

        if dl and edges_dl:
            bstate = self.get_block_state(b=zeros(self.B), layers=True)[0]
            S += bstate.entropy(dense=True, multigraph=True, dl=False,
                                norm=False)

        if dl:
            S += self.__get_layer_entropy()

        if norm:
            S /= self.E
        return S

def init_layer_confined(g, ec):
    tmp_state = CovariateBlockState(g, ec=ec, B=g.num_vertices())
    tmp_state = tmp_state.copy(overlap=True)
    be = tmp_state.get_edge_blocks()
    ba = ungroup_vector_property(be, [0])[0]
    ba.a = ba.a + ec.a * (ba.a.max() + 1)
    continuous_map(ba)
    be = group_vector_property([ba, ba])
    return be