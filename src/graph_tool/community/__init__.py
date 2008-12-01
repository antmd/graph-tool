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
dl_import("import libgraph_tool_community")

from .. core import _degree, _prop
import random, sys

__all__ = ["community_structure", "modularity"]


def community_structure(g, gamma=1.0, corr="erdos", spins=None,
                        weights=None, max_iter=0, t_range=(100.0,0.01),
                        n_spins=0, verbose=False, history_file="", seed=0):
    if spins == None:
        spins = g.new_vertex_property("int32_t")
        new_spins = True
    else:
        new_spins = False
    if seed != 0:
        seed = random.randint(0, sys.maxint)
    libgraph_tool_community.community_structure(g._Graph__graph, gamma, corr,
                                                max_iter, t_range[1],
                                                t_range[0], n_spins, new_spins,
                                                seed, verbose, history_file,
                                                _prop("e", g, weights),
                                                _prop("v", g, spins))
    return spins

def modularity(g):
    c = libgraph_tool_clustering.global_clustering(g._Graph__graph)
    return c
