#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# graph_tool -- a general graph manipulation python module
#
# Copyright (C) 2006-2013 Tiago de Paula Peixoto <tiago@skewed.de>
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

from .. import _degree, _prop, Graph, GraphView, libcore, _get_rng
import random
from numpy import *
from scipy.optimize import fsolve, fminbound
import scipy.special
from collections import defaultdict

from .. dl_import import dl_import
dl_import("from . import libgraph_tool_community as libcommunity")


class BlockState(object):
    r"""This class encapsulates the block state of a given graph.

    This must be instantiated and used by functions such as :func:`mcmc_sweep`.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    eweight : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Edge weights (i.e. multiplicity).
    vweight : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Vertex weights (i.e. multiplicity).
    b : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Initial block labels on the vertices. If not supplied, it will be
        randomly sampled.
    B : ``int`` (optional, default: ``None``)
        Number of blocks. If not supplied it will be either obtained from the
        parameter ``b``, or set to the maximum possible value according to the
        minimum description length.
    clabel : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        This parameter provides a constraint label, such that vertices with
        different labels will not be allowed to belong to the same block. If not given,
        all labels will be assumed to be the same.
    deg_corr : ``bool`` (optional, default: ``True``)
        If ``True``, the degree-corrected version of the blockmodel ensemble will
        be assumed, otherwise the traditional variant will be used.

    """

    def __init__(self, g, eweight=None, vweight=None, b=None, clabel=None, B=None,
                 deg_corr=True):
        self.g = g
        self.bg = Graph(directed=g.is_directed())
        self.mrs = self.bg.new_edge_property("int")
        self.mrp = self.bg.new_vertex_property("int")
        if g.is_directed():
            self.mrm = self.bg.new_vertex_property("int")
        else:
            self.mrm = self.mrp
        self.wr = self.bg.new_vertex_property("int")

        if eweight is None:
            eweight = g.new_edge_property("int")
            eweight.a = 1
        elif eweight.value_type() != "int32_t":
            eweight = eweight.copy(value_type="int32_t")
        if vweight is None:
            vweight = g.new_vertex_property("int")
            vweight.a = 1
        elif vweight.value_type() != "int32_t":
            vweight = vweight.copy(value_type="int32_t")
        self.eweight = eweight
        self.vweight = vweight

        self.E = self.eweight.a.sum()
        self.N = self.vweight.a.sum()

        self.deg_corr = deg_corr

        if clabel is None:
            self.clabel = g.new_vertex_property("int")
            self.L = 1
        else:
            self.clabel = clabel
            self.L = int(self.clabel.a.max() + 1)

        if b is None:
            if B is None:
                B = get_max_B(self.N, self.E, directed=g.is_directed())
            B = int(ceil(B/float(self.L)) * self.L)
            b = g.new_vertex_property("int")
            b.a = random.randint(0, B / self.L, len(b.a))
            b.a = self.clabel.a + b.a * self.L
            self.b = b
        else:
            if B is None:
                B = int(b.a.max()) + 1
            B = int(ceil(B/float(self.L)) * self.L)
            self.b = b = b.copy(value_type="int32_t")

        cg, br, vcount, ecount = condensation_graph(g, b,
                                                    vweight=vweight,
                                                    eweight=eweight,
                                                    self_loops=True)[:4]
        self.bg.add_vertex(B)
        if b.a.max() >= B:
            raise ValueError("Maximum value of b is largest or equal to B!")

        self.vertices = libcommunity.get_vector(B)
        self.vertices.a = arange(B)

        self.mrp.a = 0
        self.mrm.a = 0
        self.mrs.a = 0
        for e in cg.edges():
            r = self.bg.vertex(br[e.source()])
            s = self.bg.vertex(br[e.target()])
            be = self.bg.add_edge(r, s)
            if self.bg.is_directed():
                self.mrs[be] = ecount[e]
            else:
                self.mrs[be] = ecount[e] if r != s else 2 * ecount[e]
            self.mrp[r] += ecount[e]
            self.mrm[s] += ecount[e]

        self.wr.a = 0
        for v in cg.vertices():
            r = self.bg.vertex(br[v])
            self.wr[r] = vcount[v]

        self.__regen_emat()
        self.__build_egroups()

    def __regen_emat(self):
        self.emat = libcommunity.create_emat(self.bg._Graph__graph, len(self.vertices))

    def __build_egroups(self):
        self.esrcpos = self.g.new_edge_property("vector<int>")
        self.etgtpos = self.g.new_edge_property("vector<int>")
        self.egroups = libcommunity.build_egroups(self.g._Graph__graph,
                                                self.bg._Graph__graph,
                                                _prop("v", self.g, self.b),
                                                _prop("e", self.g, self.eweight),
                                                _prop("e", self.g, self.esrcpos),
                                                _prop("e", self.g, self.etgtpos))

    def get_blocks(self):
        r"""Returns the property map which contains the block labels for each vertex."""
        return self.b

    def get_bg(self):
        r"""Returns the block graph."""
        return self.bg

    def get_ers(self):
        r"""Returns the edge property map of the block graph which contains the :math:`e_{rs}` matrix entries."""
        return self.mrs

    def get_er(self):
        r"""Returns the vertex property map of the block graph which contains the number
        :math:`e_r` of half-edges incident on block :math:`r`. If the graph is
        directed, a pair of property maps is returned, with the number of
        out-edges :math:`e^+_r` and in-edges :math:`e^-_r`, respectively."""
        if self.bg.is_directed():
            return self.mrp. self.mrm
        else:
            return self.mrp

    def get_nr(self):
        r"""Returns the vertex property map of the block graph which contains the block sizes :math:`n_r`."""
        return self.wr

    def get_eweight(self):
        r"""Returns the block edge counts associated with the block matrix
        :math:`e_{rs}`. For directed graphs it is identical to :math:`e_{rs}`,
        but for undirected graphs it is identical except for the diagonal, which
        is :math:`e_{rr}/2`."""
        eweight = self.mrs.copy()
        if not self.g.is_directed():
            for r in self.bg.vertices():
                e = self.bg.edge(r, r)
                if e is not None:
                    eweight[e] /= 2
        return eweight

    def get_clabel(self):
        r"""Obtain the constraint label associated with each block."""
        blabel = self.bg.vertex_index.copy(value_type="int")
        blabel.a = blabel.a % self.L
        return blabel

    def entropy(self, complete=False, random=False, dl=False):
        r"""Calculate the entropy per edge associated with the current block partition.

        Parameters
        ----------
        complete : ``bool`` (optional, default: ``False``)
            If ``True``, the complete entropy will be returned, including constant
            terms not relevant to the block partition.
        random : ``bool`` (optional, default: ``False``)
            If ``True``, the entropy entropy corresponding to an equivalent random
            graph (i.e. no block partition) will be returned.
        dl : ``bool`` (optional, default: ``False``)
            If ``True``, the full description length will be returned.


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
        \sum_se_{sr}` are the number of out- and in-edges adjacent to block
        :math:`r`, respectively.

        If ``complete == False`` only the last term of the equations above will
        be returned. If ``random == True`` it will be assumed that :math:`B=1`
        despite the actual :math:`e_{rs}` matrix.  If ``dl == True``, the
        description length :math:`\mathcal{L}_t` of the model will be returned
        as well, as described in :func:`model_entropy`. Note that for the
        degree-corrected version the description length is

        .. math::

            \mathcal{L}_c = \mathcal{L}_t - N\sum_kp_k\ln p_k,

        where :math:`p_k` is the fraction of nodes with degree :math:`p_k`, and
        we have instead :math:`k \to (k^-, k^+)` for directed graphs.

        Note that the value returned corresponds to the entropy `per edge`,
        i.e. :math:`(\mathcal{S}_{t/c}\; [\,+\, \mathcal{L}_{t/c}])/ E`.

        """

        E = self.E
        N = self.N

        if self.deg_corr:
            if self.g.is_directed():
                S_rand = E * log(E)
            else:
                S_rand = E * log(2 * E)
        else:
            ak = E / float(N) if self.g.is_directed() else  2 * E / float(N)
            S_rand = E * log (N / ak)

        if random:
            S = S_rand
        else:
            S = libcommunity.entropy(self.bg._Graph__graph,
                                   _prop("e", self.bg, self.mrs),
                                   _prop("v", self.bg, self.mrp),
                                   _prop("v", self.bg, self.mrm),
                                   _prop("v", self.bg, self.wr),
                                   self.deg_corr)

        if complete:
            if self.deg_corr:
                S_seq = 0
                hist = defaultdict(int)
                for v in self.g.vertices():
                    hist[(v.in_degree(), v.out_degree())] += 1
                for k, v in hist.items():
                    p = v / float(self.g.num_vertices())
                    S_seq -= p * log(p)
                S_seq *= self.g.num_vertices()
                S += S_seq

                S -= E
                for v in self.g.vertices():
                    S -= scipy.special.gammaln(v.out_degree() + 1)
                    if self.g.is_directed():
                        S -= scipy.special.gammaln(v.in_degree() + 1)
            else:
                S += E
        else:
            S -= S_rand


        if dl:
            if random:
                S += model_entropy(1, N, E, directed=self.g.is_directed()) * E
            else:
                S += model_entropy(len(self.vertices), N, E, directed=self.g.is_directed()) * E

        return S / E

    def __min_dl(self):
        return self.entropy(complete=False, dl=True)

    def dist(self, r, s):
        r"""Compute the "distance" between blocks `r` and `s`, i.e. the entropy
        difference after they are merged together."""
        return libcommunity.dist(self.bg._Graph__graph, int(r), int(s),
                                 _prop("e", self.bg, self.mrs),
                                 _prop("v", self.bg, self.mrp),
                                 _prop("v", self.bg, self.mrm),
                                 _prop("v", self.bg, self.wr),
                                 self.deg_corr)


    def join(self, r, s):
        r"""Merge blocks `r` and `s` into a single block."""
        libcommunity.join(self.bg._Graph__graph, int(r), int(s),
                          _prop("e", self.bg, self.mrs),
                          _prop("v", self.bg, self.mrp),
                          _prop("v", self.bg, self.mrm),
                          _prop("v", self.bg, self.wr),
                          self.deg_corr,
                          self.vertices)
        del self.egroups

    def remove_vertex(self, v):
        r"""Remove vertex `v` from its current block."""
        libcommunity.remove_vertex(self.g._Graph__graph,
                                   self.bg._Graph__graph,
                                   int(v),
                                   _prop("e", self.bg, self.mrs),
                                   _prop("v", self.bg, self.mrp),
                                   _prop("v", self.bg, self.mrm),
                                   _prop("v", self.bg, self.wr),
                                   _prop("v", self.g, self.b))
        del self.egroups


    def add_vertex(self, v, r):
        r"""Add vertex `v` to block `r`."""
        libcommunity.add_vertex(v.get_graph()._Graph__graph,
                                self.bg._Graph__graph,
                                int(v), int(r),
                                _prop("e", self.bg, self.mrs),
                                _prop("v", self.bg, self.mrp),
                                _prop("v", self.bg, self.mrm),
                                _prop("v", self.bg, self.wr),
                                _prop("v", self.g, self.b))
        del self.egroups

    def move_vertex(self, v, nr):
        r"""Move vertex `v` to block `r`, and return the entropy difference."""
        dS = libcommunity.move_vertex(self.g._Graph__graph,
                                      self.bg._Graph__graph,
                                      self.emat,
                                      int(v), int(nr),
                                      _prop("e", self.bg, self.mrs),
                                      _prop("v", self.bg, self.mrp),
                                      _prop("v", self.bg, self.mrm),
                                      _prop("v", self.bg, self.wr),
                                      _prop("v", self.g, self.b),
                                      self.deg_corr,
                                      _prop("e", self.bg, self.eweight),
                                      _prop("v", self.bg, self.vweight))
        del self.egroups
        return dS / float(self.E)


    def get_dist_matrix(self):
        r"""Return the distance matrix between all blocks. The distance is
        defined as the entropy difference after two blocks are merged."""
        dist_matrix = zeros((len(self.vertices), len(self.vertices)))
        for ri, r in enumerate(self.vertices):
            for si, s in enumerate(self.vertices):
                if si > ri:
                    dist_matrix[ri, si] = self.dist(r, s)
        for ri, r in enumerate(self.vertices):
            for si, s in enumerate(self.vertices):
                if si < ri:
                    dist_matrix[ri, si] = dist_matrix[si, ri]
        return dist_matrix


    def get_matrix(self, reorder=False, clabel=None, niter=0, ret_order=False):
        r"""Returns the block matrix.

        Parameters
        ----------
        reorder : ``bool`` (optional, default: ``False``)
            If ``True``, the matrix is reordered so that blocks which are
            'similar' are close together.
        clabel : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
            Constraint labels to be imposed during reordering. Only has
            effect if ``reorder == True``.
        niter : ``int`` (optional, default: `0`)
            Number of iterations performed to obtain the best ordering. If
            ``niter == 0`` it will automatically determined. Only has effect
            if ``reorder == True``.
        ret_order : ``bool`` (optional, default: ``False``)
            If ``True``, the vertex ordering is returned. Only has effect if
            ``reorder == True``.

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
           >>> m = state.get_matrix(reorder=True)
           >>> figure()
           <...>
           >>> matshow(m)
           <...>
           >>> savefig("bloc_mat.pdf")

        .. testcleanup:: get_matrix

           savefig("bloc_mat.png")

        .. figure:: bloc_mat.*
           :align: center

           A  5x5 block matrix.

       """
        B = len(self.vertices)
        vmap = {}
        for r in range(len(self.vertices)):
            vmap[self.vertices[r]] = r

        if reorder:
            if niter == 0:
                niter = max(10 * len(self.vertices), 100)

            states = []

            label = None
            states = [self]
            Bi = B / 2

            while Bi > 1:
                cblabel = states[-1].get_clabel()
                if clabel is not None and len(states) == 1:
                    clabel.a += (cblabel.a.max() + 1) * clabel.a
                state = BlockState(states[-1].bg, B=Bi,
                                   clabel=clabel,
                                   vweight=states[-1].wr,
                                   eweight=states[-1].get_eweight(),
                                   deg_corr=states[-1].deg_corr)

                for i in range(niter):
                    mcmc_sweep(state, beta=1.)
                for i in range(niter):
                    mcmc_sweep(state, beta=float("inf"))

                states.append(state)

                Bi /= 2

                if Bi < cblabel.a.max() + 1:
                    break


            vorder = list(range(len(states[-1].vertices)))
            for state in reversed(states[1:]):
                norder = [[] for i in range(len(state.vertices))]
                #print(state.L)
                for v in state.g.vertices():
                    pos = vorder.index(state.b[v])
                    norder[pos].append(int(v))
                vorder = [item for sublist in norder for item in sublist]
        else:
            vorder = self.vertices

        order_map = zeros(B, dtype="int")
        for i, v in enumerate(vorder):
            order_map[vmap[v]] = i

        m = zeros((B, B))
        rmap = {}
        for e in self.bg.edges():
            r, s = vmap[int(e.source())], vmap[int(e.target())]
            r = order_map[r]
            s = order_map[s]
            rmap[r] = int(e.source())
            rmap[s] = int(e.target())
            m[r, s] = self.mrs[e]
            if not self.bg.is_directed():
                m[s, r] = m[r, s]

        if ret_order:
            return m, rmap
        else:
            return m

    def modularity(self):
        r"""Computes the modularity of the current block structure."""
        Q = 0
        for vi in self.vertices:
            v = self.bg.vertex(vi)
            err = 0
            for e in v.out_edges():
                if e.target() == v:
                    err = self.mrs[e]
                    break
            er = 0
            if self.bg.is_directed():
                err /= float(self.E)
                er = self.mrp[v] * self.mrm[v] / float(self.E ** 2)
            else:
                err /= float(2 * self.E)
                er = self.mrp[v] ** 2 / float(4 * self.E ** 2)
            Q += err - er
        return Q

def model_entropy(B, N, E, directed):
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

       \mathcal{L}_t = \ln\Omega_m + N\ln B,


    where :math:`N\ln B` is the information necessary to describe the
    block partition.

    References
    ----------

    .. [peixoto-parsimonious-2012] Tiago P. Peixoto "Parsimonious module inference in large networks",
       (2012), :arxiv:`1212.4794`.

    """
    return libcommunity.SB(float(B), int(N), int(E), directed)


def Sdl(B, S, N, E, directed=False):
    return S + model_entropy(B, N, E, directed)


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
    .. [peixoto-parsimonious-2012] Tiago P. Peixoto "Parsimonious module inference in large networks",
       (2012), :arxiv:`1212.4794`.


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
    2.4199998721289937

    References
    ----------
    .. [peixoto-parsimonious-2012] Tiago P. Peixoto "Parsimonious module inference in large networks",
       (2012), :arxiv:`1212.4794`.

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

def min_dist(state, n=0):
    r"""Return the minimum distance between all blocks, and the block pair which minimizes it.

    The parameter `state` must be an instance of the
    :class:`~graph_tool.community.BlockState` class, and `n` is the number of
    block pairs to sample. If `n == 0` all block pairs are sampled.


    Examples
    --------

    .. testsetup:: min_dist

       gt.seed_rng(42)
       np.random.seed(42)

    .. doctest:: min_dist

       >>> g = gt.collection.data["polbooks"]
       >>> state = gt.BlockState(g, B=4, deg_corr=True)
       >>> for i in range(1000):
       ...     ds, nmoves = gt.mcmc_sweep(state)
       >>> gt.min_dist(state)
       (795.7694502418635, 2, 3)

    """
    min_d, r, s = libcommunity.min_dist(state.bg._Graph__graph, int(n),
                                        _prop("e", state.bg, state.mrs),
                                        _prop("v", state.bg, state.mrp),
                                        _prop("v", state.bg, state.mrm),
                                        _prop("v", state.bg, state.wr),
                                        state.deg_corr,
                                        state.vertices,
                                        _get_rng())
    return min_d, r, s


def mcmc_sweep(state, beta=1., sequential=True, verbose=False, vertices=None):
    r"""Performs a Monte Carlo Markov chain sweep on the network, to sample the block partition according to a probability :math:`\propto e^{-\beta \mathcal{S}_{t/c}}`, where :math:`\mathcal{S}_{t/c}` is the blockmodel entropy.

    Parameters
    ----------
    state : :class:`~graph_tool.community.BlockState`
        The block state.
    beta : `float` (optional, default: `1.0`)
        The inverse temperature parameter :math:`\beta`.
    sequential : ``bool`` (optional, default: ``True``)
        If ``True``, the move attempts on the vertices are done in sequential
        random order. Otherwise a total of `N` moves attempts are made, where
        `N` is the number of vertices, where each vertex can be selected with
        equal probability.
    verbose : ``bool`` (optional, default: ``False``)
        If ``True``, verbose information is displayed.

    Returns
    -------

    dS : `float`
       The entropy difference (per edge) after a full sweep.
    nmoves : ``int``
       The number of accepted block membership moves.


    Notes
    -----

    This algorithm performs a Monte Carlo Markov chain sweep on the network,
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
    the markov chain by proposing membership moves :math:`r\to s` with
    probability :math:`p(r\to s|t) \propto e_{ts} + 1`, where :math:`t` is the
    block label of a random neighbour of the vertex being moved. See
    [peixoto-parsimonious-2012]_ for more details.

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
    .. [peixoto-parsimonious-2012] Tiago P. Peixoto "Parsimonious module inference in large networks",
       (2012), :arxiv:`1212.4794`.
    """

    if len(state.vertices) == 1:
        return 0., 0

    if vertices is None:
        u = state.g
        u.stash_filter()
        vertices = libcommunity.get_vector(u.num_vertices())
        u.pop_filter()
        vertices.a = arange(len(vertices.a))

    dS, nmoves = libcommunity.move_sweep(state.g._Graph__graph,
                                         state.bg._Graph__graph,
                                         state.emat,
                                         _prop("e", state.bg, state.mrs),
                                         _prop("v", state.bg, state.mrp),
                                         _prop("v", state.bg, state.mrm),
                                         _prop("v", state.bg, state.wr),
                                         _prop("v", state.g, state.b),
                                         _prop("v", state.g, state.clabel),
                                         state.L,
                                         vertices,
                                         state.deg_corr,
                                         _prop("e", state.g, state.eweight),
                                         _prop("v", state.g, state.vweight),
                                         state.egroups,
                                         _prop("e", state.g, state.esrcpos),
                                         _prop("e", state.g, state.etgtpos),
                                         float(beta), sequential,
                                         verbose, _get_rng())
    return dS / state.E, nmoves




def mc_get_dl(state, nsweep, greedy, rng, checkpoint=None, anneal=1,
              verbose=False):

    if len(state.vertices) == 1:
        return state._BlockState__min_dl()

    S = state._BlockState__min_dl()

    if nsweep >= 0:
        for i in range(nsweep):
            delta, nmoves = mcmc_sweep(state, beta=1, rng=rng)
            if checkpoint is not None:
                checkpoint(state, S, delta, nmoves)
            S += delta
        if greedy:
            for i in range(nsweep):
                delta, nmoves = mcmc_sweep(state, beta=float("inf"))
                if checkpoint is not None:
                    checkpoint(state, S, delta, nmoves)
                    S += delta
    else:
        # adaptive mode
        min_dl = S
        max_dl = S
        count = 0
        time = 0
        bump = False

        beta = 1.

        last_min = min_dl

        if verbose:
            print("beta = %g" % beta)
        while True:
            delta, nmoves = mcmc_sweep(state, beta=beta)
            if checkpoint is not None:
                checkpoint(state, S, delta, nmoves)
            S += delta

            if S < min_dl:
                min_dl = S
                count = 0
            elif S > max_dl:
                max_dl = S
                count = 0
            else:
                count += 1

            if count > abs(nsweep):
                if not bump:
                    min_dl = max_dl = S
                    bump = True
                    count = 0
                else:
                    if anneal <= 1 or min_dl == last_min:
                        break
                    else:
                        beta *= anneal
                        count = 0
                        last_min = min_dl
                        if verbose:
                            print("beta = %g" % beta)

        min_dl = S
        count = 0
        while True:
            delta, nmoves = mcmc_sweep(state, beta=float("inf"))
            if checkpoint is not None:
                checkpoint(state, S, delta, nmoves)
            S += delta
            if S < min_dl:
                min_dl = S
                count = 0
            else:
                count += 1
            if count > abs(nsweep):
                break

    return state._BlockState__min_dl()

def get_b_dl(g, bs, bs_start, B, nsweep, anneal, greedy, clabel, deg_corr, rng,
             checkpoint=None, verbose=False):
    prev_dl = float("inf")
    if B in bs:
        return bs[B][0]
    elif B in bs_start:
        if verbose:
            print("starting from previous result for B=%d" % B)
        prev_dl, b = bs_start[B]
        state = BlockState(g, b=b.copy(), clabel=clabel, deg_corr=deg_corr)
    else:
        n_iter = 0
        bs_keys = [k for k in bs.keys() if type(k) != str]
        B_sup = max(bs_keys) if len(bs_keys) > 0 else B
        for Bi in bs_keys:
            if Bi > B and Bi < B_sup:
                B_sup = Bi
        if B_sup == B:
            state = BlockState(g, B=B, clabel=clabel, deg_corr=deg_corr)
            b = state.b
        else:
            if verbose:
                print("shrinking from", B_sup, "to", B)
            b = bs[B_sup][1].copy()

            cg, br, vcount, ecount = condensation_graph(g, b, self_loops=True)[:4]

            blabel = cg.new_vertex_property("int")
            if clabel is not None:
                blabel.a = br.a % (clabel.a.max() + 1)

            bg_state = BlockState(cg, B=B, clabel=blabel, vweight=vcount, eweight=ecount,
                                  deg_corr=deg_corr)

            mc_get_dl(bg_state, nsweep=nsweep, greedy=greedy, rng=rng,
                      checkpoint=checkpoint, anneal=anneal, verbose=verbose)

            ### FIXME: the following could be improved by moving it to the C++
            ### side
            bmap = {}
            for v in bg_state.g.vertices():
               bmap[br[v]] = v
            for v in g.vertices():
                b[v] = bg_state.b[bmap[b[v]]]

            state = BlockState(g, b=b, B=B, clabel=clabel, deg_corr=deg_corr)

    dl = mc_get_dl(state, nsweep=nsweep,  greedy=greedy, rng=rng,
                   checkpoint=checkpoint, anneal=anneal, verbose=verbose)

    if dl < prev_dl:
        bs[B] = [dl, state.b.copy()]
    else:
        bs[B] = bs_start[B]
        dl = prev_dl
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

def minimize_blockmodel_dl(g, deg_corr=True, nsweeps=100, adaptive_convergence=True,
                           anneal=1., greedy_cooling=True, max_B=None, min_B=1,
                           mid_B=None, b_cache=None, b_start=None, clabel=None,
                           checkpoint=None, verbose=False):
    r"""Find the block partition of an unspecified size which minimizes the description
    length of the network, according to the stochastic blockmodel ensemble which
    best describes it.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph being used.
    nsweeps : ``int`` (optional, default: `50`)
        Number of sweeps per value of `B` tried. If `adaptive_convergence ==
        True`, this corresponds to the number of sweeps observed to determine
        convergence, not the total number of sweeps performed, which can be much
        larger.
    adaptive_convergence : ``bool`` (optional, default: ``True``)
        If ``True``, the parameter ``nsweeps`` represents not the total number
        of sweeps performed per value of ``B``, but instead the number of sweeps
        without an improvement on the value of :math:`S_{c/t}` so that
        convergence is assumed.
    anneal : ``float`` (optional, default: ``1.``)
       Annealing factor which multiplies the inverse temperature :math:`\beta`
       after each convergence. If ``anneal <= 1.``, no annealing is performed.
    greedy_cooling : ``bool`` (optional, default: ``True``)
        If ``True``, a final abrupt cooling step is performed after the Markov
        chain has equilibrated.
    max_B : ``int`` (optional, default: ``None``)
        Maximum number of blocks tried. If not supplied, it will be
        automatically determined.
    min_B : ``int`` (optional, default: `1`)
        Minimum number of blocks tried.
    mid_B : ``int`` (optional, default: ``None``)
        Middle of the range which brackets the minimum. If not supplied, will be
        automatically determined.
    b_cache : :class:`dict` with ``int`` keys and (``float``, :class:`~graph_tool.PropertyMap`) values (optional, default: ``None``)
        If provided, this corresponds to a dictionary where the keys are the
        number of blocks, and the values are tuples containing two values: the
        description length and its associated vertex partition. Values present
        in this dictionary will not be computed, and will be used unmodified as
        the solution for the corresponding number of blocks. This can be used to
        continue from a previously unfinished run.
    b_start : :class:`dict` with ``int`` keys and (``float``, :class:`~graph_tool.PropertyMap`) values (optional, default: ``None``)
        Like `b_cache`, but the partitions present in the dictionary will be
        used as the starting point of the minimization.
    clabel : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Constraint labels on the vertices, such that vertices with different
        labels cannot belong to the same block.
    checkpoint : function (optional, default: ``None``)
        If provided, this function will be called after each call to
        :func:`mcmc_sweep`. This can be used to store the current state, so it
        can be continued later. The function must have the following signature:

        .. code-block:: python

            def checkpoint(state, L, delta, nmoves):
                ...

        where `state` is either a :class:`~graph_tool.community.BlockState`
        instance or ``None``, `L` is the current description length, `delta` is
        the entropy difference in the last MCMC sweep, and `nmoves` is the
        number of accepted block membership moves.

        This function will also be called when the MCMC has finished for the
        current value of :math:`B`, in which case ``state == None``, and the
        remaining parameters will be zero.
    verbose : ``bool`` (optional, default: ``False``)
        If ``True``, verbose information is displayed.

    Returns
    -------
    b : :class:`~graph_tool.PropertyMap`
       Vertex property map with the best block partition.
    min_dl : `float`
       Minimum value of the description length (in `nats <http://en.wikipedia.org/wiki/Nat_%28information%29>`_).
    b_cache : :class:`dict` with ``int`` keys and (``float``, :class:`~graph_tool.PropertyMap`) values
        Dictionary where the keys are the number of blocks visited during the
        algorithm, and the values are tuples containing two values: the
        description length and its associated vertex partition.

    Notes
    -----

    This algorithm attempts to find a block partition of an unspecified size
    which minimizes the description length of the network,

    .. math::

       \Sigma_{t/c} = \mathcal{S}_{t/c} + \mathcal{L}_{t/c},

    where :math:`\mathcal{S}_{t/c}` is the blockmodel entropy (as described in
    the docstring of :func:`mcmc_sweep` and :meth:`BlockState.entropy`) and
    :math:`\mathcal{L}_{t/c}` is the information necessary to describe the model
    (as described in the docstring of :func:`model_entropy` and
    :meth:`BlockState.entropy`).

    The algorithm works by minimizing the entropy :math:`\mathcal{S}_{t/c}` for
    specific values of :math:`B` via :func:`mcmc_sweep` (with :math:`\beta = 1`
    and :math:`\beta\to\infty`), and minimizing :math:`\Sigma_{t/c}` via an
    one-dimensional Fibonacci search on :math:`B`. See
    [peixoto-parsimonious-2012]_ for more details.

    This algorithm has a complexity of :math:`O(\tau E\ln B_{\text{max}})`,
    where :math:`E` is the number of edges in the network, :math:`\tau` is the
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
       >>> b, mdl, b_cache = gt.minimize_blockmodel_dl(g)
       >>> gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=b, vertex_shape=b, output="polbooks_blocks_mdl.pdf")
       <...>

    .. testcleanup:: mdl

       gt.graph_draw(g, pos=g.vp["pos"], vertex_fill_color=b, vertex_shape=b, output="polbooks_blocks_mdl.png")

    .. figure:: polbooks_blocks_mdl.*
       :align: center

       Block partition of a political books network, which minimizes the description
       lenght of the network according to the degree-corrected stochastic blockmodel.

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
    .. [peixoto-parsimonious-2012] Tiago P. Peixoto "Parsimonious module inference in large networks",
       (2012), :arxiv:`1212.4794`.

    """

    rng = _get_rng()

    if adaptive_convergence:
        nsweeps = -abs(nsweeps)

    if max_B is None:
        max_B = get_max_B(g.num_vertices(), g.num_edges(), g.is_directed())
        if verbose:
            print("max_B:", max_B)
    if min_B is None:
        min_B = 1

    if mid_B is None:
        mid_B = get_mid(min_B, max_B)


    greedy = greedy_cooling
    if b_start is None:
        b_start = {}
    bs = b_cache
    if bs is None:
        bs = {}

    while True:
        f_max = get_b_dl(g, bs, b_start, max_B, nsweeps, anneal, greedy, clabel,
                         deg_corr, rng, checkpoint, verbose)
        f_mid = get_b_dl(g, bs, b_start, mid_B, nsweeps, anneal, greedy, clabel,
                         deg_corr, rng, checkpoint, verbose)
        f_min = get_b_dl(g, bs, b_start, min_B, nsweeps, anneal, greedy, clabel,
                         deg_corr, rng, checkpoint, verbose)

        if verbose:
            print("bracket:", min_B, mid_B, max_B, f_min, f_mid, f_max)

        if checkpoint is not None:
            checkpoint(None, 0, 0, 0)

        if f_max > f_mid > f_min:
            max_B = mid_B
            mid_B = get_mid(min_B, mid_B)
            if verbose:
                print("bracket:", min_B, mid_B, max_B, f_min, f_mid, f_max)
        elif f_max < f_mid < f_min:
            min_d = mid_B
            mid_B = get_mid(mid_B, max_B)
            if verbose:
                print("bracket:", min_B, mid_B, max_B, f_min, f_mid, f_max)
        else:
            break


    while True:
        if max_B - mid_B > mid_B - min_B:
            x = get_mid(mid_B, max_B)
        else:
            x = get_mid(min_B, mid_B)

        f_x = get_b_dl(g, bs, b_start, x, nsweeps, anneal, greedy, clabel,
                       deg_corr, rng, checkpoint, verbose)
        f_mid = get_b_dl(g, bs, b_start, mid_B, nsweeps, anneal, greedy, clabel,
                         deg_corr, rng, checkpoint, verbose)

        if verbose:
            print("bisect: (", min_B, mid_B, max_B, ") ->", x, f_x) #, is_fibo((mid_B - min_B)), is_fibo((max_B - mid_B)))

        if max_B - mid_B <= 1:
            min_dl = float(inf)
            best_B = None
            for Bi in bs.keys():
                if bs[Bi][0] < min_dl:
                    min_dl = bs[Bi][0]
                    best_B = Bi
            if verbose:
                print("best:", best_B, min_dl)
            return bs[best_B][1], bs[best_B][0], bs

        if checkpoint is not None:
            checkpoint(None, 0, 0, 0)

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
       17.17007316834295
    """

    if p is None:
        p = state.g.new_edge_property("vector<int>")

    libcommunity.edge_marginals(state.g._Graph__graph,
                                state.bg._Graph__graph,
                                len(state.vertices),
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
       20.001666525168361
       >>> gt.graph_draw(g, pos=g.vp["pos"], vertex_shape="pie", vertex_pie_fractions=pv, output="polbooks_blocks_soft_B4.pdf")
       <...>

    .. testcleanup:: cvm

       gt.graph_draw(g, pos=g.vp["pos"], vertex_shape="pie", vertex_pie_fractions=pv, output="polbooks_blocks_soft_B4.png")

    .. figure:: polbooks_blocks_soft_B4.*
       :align: center

       "Soft" block partition of a political books network with :math:`B=4`.

    """
    B = len(state.vertices)

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
    B = len(state.vertices)
    H = 0
    pv =  state.g.new_vertex_property("vector<double>")

    H, sH, Hmf, sHmf  = libcommunity.bethe_entropy(state.g._Graph__graph,
                                                 len(state.vertices),
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
        Graph to be used.
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
       :align: left

       Block partition of a political books network with :math:`B=5`.

    .. figure:: polbooks_blocks_B5_cond.*
       :align: right

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
    u = GraphView(g, directed=True)
    libcommunity.community_network_vavg(u._Graph__graph,
                                        gp._Graph__graph,
                                        _prop("v", g, prop),
                                        _prop("v", gp, cprop),
                                        _prop("v", gp, vcount),
                                        _prop("v", g, vweight),
                                        avp)
    libcommunity.community_network_eavg(u._Graph__graph,
                                        gp._Graph__graph,
                                        _prop("v", g, prop),
                                        _prop("v", gp, cprop),
                                        _prop("e", gp, ecount),
                                        _prop("e", g, eweight),
                                        aep,
                                        self_loops)
    return gp, cprop, vcount, ecount, r_avp, r_aep
