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
``graph_tool.util`` - Graph Utilities
-------------------------------------
"""


from .. dl_import import dl_import
dl_import("import libgraph_tool_util")

from .. core import _degree, _prop

__all__ = ["find_vertex", "find_vertex_range", "find_edge", "find_edge_range"]

def find_vertex(g, deg, match):
    ret = libgraph_tool_util.\
          find_vertex_range(g._Graph__graph, _degree(g, deg), (match, match))
    return ret

def find_vertex_range(g, deg, range):
    ret = libgraph_tool_util.\
          find_vertex_range(g._Graph__graph, _degree(g, deg), range)
    return ret

def find_edge(g, prop, match):
    ret = libgraph_tool_util.\
          find_edge_range(g._Graph__graph, _prop("e", g, prop), (match, match))
    return ret

def find_edge_range(g, prop, range):
    ret = libgraph_tool_util.\
          find_edge_range(g._Graph__graph, _prop("e", g, prop), range)
    return ret
