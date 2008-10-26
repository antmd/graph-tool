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
dl_import("import libgraph_tool_correlations")

from .. core import _degree, _prop
from numpy import *

__all__ = ["assortativity", "scalar_assortativity", "corr_hist",
           "combined_corr_hist", "avg_neighbour_corr", "avg_combined_corr"]

def assortativity(g, deg):
    return libgraph_tool_correlations.\
           assortativity_coefficient(g._Graph__graph, _degree(g, deg))

def scalar_assortativity(g, deg):
    return libgraph_tool_correlations.\
           scalar_assortativity_coefficient(g._Graph__graph,
                                            _degree(g, deg))

def corr_hist(g, deg1, deg2, bins=[[1],[1]], weight=None, float_count=True):
    ret = libgraph_tool_correlations.\
          vertex_correlation_histogram(g._Graph__graph, _degree(g, deg1),
                                       _degree(g, deg2), _prop("e", g, weight),
                                       bins[0], bins[1])
    return [array(ret[0], dtype="float64") if float_count else ret[0],
            [ret[1][0], ret[1][1]]]

def combined_corr_hist(g, deg1, deg2, bins=[[1],[1]], float_count=True):
    ret = libgraph_tool_correlations.\
          vertex_combined_correlation_histogram(g._Graph__graph,
                                                _degree(g, deg1),
                                                _degree(g, deg2),
                                                bins[0], bins[1])
    return [array(ret[0], dtype="float64") if float_count else ret[0],
            [ret[1][0], ret[1][1]]]

def avg_neighbour_corr(g, deg1, deg2, bins=[1], weight=None):
    ret = libgraph_tool_correlations.\
          vertex_avg_correlation(g._Graph__graph, _degree(g, deg1),
                                 _degree(g, deg2), _prop("e", g, weight),
                                 bins)
    return [ret[0], ret[1], ret[2][0]]

def avg_combined_corr(g, deg1, deg2, bins=[1]):
    ret = libgraph_tool_correlations.\
          vertex_avg_combined_correlation(g._Graph__graph, _degree(g, deg1),
                                          _degree(g, deg2), bins)
    return [ret[0], ret[1], ret[2][0]]
