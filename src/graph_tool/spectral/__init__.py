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
``graph_tool.spectral`` - Spectral properties
---------------------------------------------
"""

from .. core import _degree, _prop, Graph, _limit_args
from numpy import *
import scipy.sparse


__all__ = ["adjacency", "laplacian", "incidence"]

def adjacency(g, sparse=True, weight=None):
    if g.get_vertex_filter()[0] != None:
        index = g.new_vertex_property("int64_t")
        for i,v in enumerate(g.vertices()):
            index[v] = i
    else:
        index = g.vertex_index
    N = g.num_vertices()
    if sparse:
        m = scipy.sparse.lil_matrix((N,N))
    else:
        m = matrix(zeros((N,N)))
    for v in g.vertices():
        for e in v.out_edges():
            m[index[v],index[e.target()]] = 1 if weight == None else weight[e]
    if sparse:
        m = m.tocsr()
    return m

def _get_deg(v, deg, weight):
    if deg == "total":
        if weight == None:
            d = v.in_degree() + v.out_degree()
        else:
            d = sum(weight[e] for e in v.all_edges())
    elif deg == "in":
        if weight == None:
            d = v.in_degree()
        else:
            d = sum(weight[e] for e in v.in_edges())
    else:
        if weight == None:
            d = v.out_degree()
        else:
            d = sum(weight[e] for e in v.out_edges())
    return d

@_limit_args({"deg":["total", "in", "out"]})
def laplacian(g, deg="total", normalized=True, sparse=True, weight=None):
    if g.get_vertex_filter()[0] != None:
        index = g.new_vertex_property("int64_t")
        for i,v in enumerate(g.vertices()):
            index[v] = i
    else:
        index = g.vertex_index
    N = g.num_vertices()
    if sparse:
        m = scipy.sparse.lil_matrix((N,N))
    else:
        m = matrix(zeros((N,N)))
    for v in g.vertices():
        d = _get_deg(v, deg, weight)
        if not normalized:
            m[index[v], index[v]] = d
        elif d > 0:
            m[index[v], index[v]] = 1
        for e in v.out_edges():
            if not normalized:
                m[index[v],index[e.target()]] = (-1 if weight == None
                                                 else -weight[e])
            else:
                val = (d*_get_deg(e.target(),deg,weight))**(-0.5)
                m[index[v],index[e.target()]] = val
    if sparse:
        m = m.tocsr()
    return m

def incidence(g, sparse=True):
    if g.get_vertex_filter()[0] != None:
        index = g.new_vertex_property("int64_t")
        for i,v in enumerate(g.vertices()):
            index[v] = i
    else:
        index = g.vertex_index

    eindex = g.new_edge_property("int64_t")
    for i, e in enumerate(g.edges()):
        eindex[e] = i

    N = g.num_vertices()
    E = g.num_edges()
    if sparse:
        m = scipy.sparse.lil_matrix((N,E))
    else:
        m = matrix(zeros((N,E)))
    for v in g.vertices():
        if g.is_directed():
            for e in v.out_edges():
                m[index[v],eindex[e]] += -1
            for e in v.in_edges():
                m[index[v],eindex[e]] += 1
        else:
            for e in v.out_edges():
                m[index[v],eindex[e]] += 1
    if sparse:
        m = m.tocsr()
    return m
