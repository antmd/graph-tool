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
``graph_tool.flow`` - Maximum flow algorithms
---------------------------------------------
"""

from .. dl_import import dl_import
dl_import("import libgraph_tool_flow")

from .. core import _prop, _check_prop_scalar, _check_prop_writable
__all__ = ["edmonds_karp_max_flow", "push_relabel_max_flow",
           "kolmogorov_max_flow", "max_cardinality_matching"]

def edmonds_karp_max_flow(g, source, target, capacity, residual=None):
    _check_prop_scalar(capacity, "capacity")
    if residual == None:
        residual = g.new_edge_property(capacity.value_type())
    _check_prop_scalar(residual, "residual")
    _check_prop_writable(residual, "residual")

    libgraph_tool_flow.\
           edmonds_karp_max_flow(g._Graph__graph, int(source), int(target),
                                 _prop("e", g, capacity),
                                 _prop("e", g, residual))
    return residual

def push_relabel_max_flow(g, source, target, capacity, residual=None):
    _check_prop_scalar(capacity, "capacity")
    if residual == None:
        residual = g.new_edge_property(capacity.value_type())
    _check_prop_scalar(residual, "residual")
    _check_prop_writable(residual, "residual")

    libgraph_tool_flow.\
           push_relabel_max_flow(g._Graph__graph, int(source), int(target),
                                 _prop("e", g, capacity),
                                 _prop("e", g, residual))
    return residual

def kolmogorov_max_flow(g, source, target, capacity, residual=None):
    _check_prop_scalar(capacity, "capacity")
    if residual == None:
        residual = g.new_edge_property(capacity.value_type())
    _check_prop_scalar(residual, "residual")
    _check_prop_writable(residual, "residual")

    libgraph_tool_flow.\
           kolmogorov_max_flow(g._Graph__graph, int(source), int(target),
                                 _prop("e", g, capacity),
                                 _prop("e", g, residual))
    return residual

def max_cardinality_matching(g, match=None):
    if match == None:
        match = g.new_edge_property("bool")
    _check_prop_scalar(match, "match")
    _check_prop_writable(match, "match")

    g.stash_filter(directed=True)
    g.set_directed(False)
    check = libgraph_tool_flow.\
            max_cardinality_matching(g._Graph__graph, _prop("e", g, match))
    g.pop_filter(directed=True)
    return match, check
