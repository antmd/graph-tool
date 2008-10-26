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
dl_import("import libgraph_tool_clustering")

from .. core import _degree, _prop
from numpy import *

__all__ = ["local_clustering", "global_clustering", "extended_clustering"]

def local_clustering(g, prop=None):
    if prop == None:
        prop = g.new_vertex_property("double")
    libgraph_tool_clustering.extended_clustering(g._Graph__graph,
                                                 [_prop("v", g, prop)])
    return prop

def global_clustering(g):
    c = libgraph_tool_clustering.global_clustering(g._Graph__graph)
    return c

def extended_clustering(g, props=None, max_depth=3):
    if props == None:
        props = []
        for i in xrange(0, max_depth):
            props.append(g.new_vertex_property("double"))
    libgraph_tool_clustering.extended_clustering(g._Graph__graph,
                                                 [_prop("v", g, p) for p in props])
    return props
