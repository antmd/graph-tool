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

from .. dl_import import dl_import
dl_import("import libgraph_tool_stats")

from .. core import _degree, _prop
from numpy import *

__all__ = ["vertex_hist", "edge_hist", "label_components",
           "label_parallel_edges",  "label_self_loops"]

def vertex_hist(g, deg, bins=[1], float_count=True):
    ret = libgraph_tool_stats.\
          get_vertex_histogram(g._Graph__graph, _degree(g, deg), bins)
    return [array(ret[0], dtype="float64") if float_count else ret[0], ret[1]]

def edge_hist(g, eprop, bins=[1], float_count=True):
    ret = libgraph_tool_stats.\
          get_edge_histogram(g._Graph__graph, _prop("e", g, eprop), bins)
    return [array(ret[0], dtype="float64") if float_count else ret[0], ret[1]]

def label_components(g, vprop=None):
    if vprop == None:
        vprop = g.new_vertex_property("int32_t")
    libgraph_tool_stats.\
          label_components(g._Graph__graph, _prop("v", g, vprop))
    return vprop

def label_parallel_edges(g, eprop):
    if eprop == None:
        eprop = g.new_edge_property("int32_t")
    libgraph_tool_stats.\
          label_parallel_edges(g._Graph__graph, _prop("e", g, eprop))
    return eprop

def label_self_loops(g, eprop):
    if eprop == None:
        eprop = g.new_edge_property("int32_t")
    libgraph_tool_stats.\
          label_self_loops(g._Graph__graph, _prop("e", g, eprop))
    return eprop
