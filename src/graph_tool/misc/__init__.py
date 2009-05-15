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
``graph_tool.misc`` - Miscellaneous functions
---------------------------------------------
"""

from .. dl_import import dl_import
dl_import("import libgraph_tool_misc")

from .. core import _prop
import random, sys
__all__ = ["random_rewire", "isomorphism"]

def random_rewire(g, strat="uncorrelated", self_loops = False,
                  parallel_edges = False, seed = 0):
    if seed != 0:
        seed = random.randint(0, sys.maxint)
    if g.is_reversed():
        was_reversed = True
    else:
        was_reversed = False
    g.set_reversed(False)
    libgraph_tool_misc.random_rewire(g._Graph__graph, strat, self_loops,
                                     parallel_edges, seed)
    if was_reversed:
        g.set_reversed(True)

def isomorphism(g1, g2, isomap=None):
    if isomap == None:
        isomap = g1.new_vertex_property("int32_t")
    return libgraph_tool_misc.\
           check_isomorphism(g1._Graph__graph,g2._Graph__graph,
                             _prop("v", g1, isomap))
