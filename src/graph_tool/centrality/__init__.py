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
dl_import("import libgraph_tool_centrality")

from .. core import _prop
import numpy

__all__ = ["pagerank", "betweenness", "central_point_dominance", "eigentrust",
           "absolute_trust"]

def pagerank(g, damping=0.8, prop=None, epslon=1e-6, max_iter=0,
             ret_iter=False):
    if prop == None:
        prop = g.new_vertex_property("double")
    ic = libgraph_tool_centrality.\
            get_pagerank(g._Graph__graph, _prop("v", g, prop), damping, epslon,
                         max_iter)
    if ret_iter:
        return prop, ic
    else:
        return prop

def betweenness(g, eprop=None, vprop=None, weight=None, norm=True):
    if vprop == None:
        vprop = g.new_vertex_property("double")
    if eprop == None:
        eprop = g.new_edge_property("double")
    if weight != None and weight.value_type() != eprop.value_type():
        nw = g.new_edge_property(eprop.value_type())
        g.copy_property(weight, nw)
        weight = nw
    libgraph_tool_centrality.\
            get_betweenness(g._Graph__graph, _prop("e", g, weight),
                            _prop("e", g, eprop), _prop("v", g, vprop), norm)
    return vprop, eprop

def central_point_dominance(g, betweenness):
    return libgraph_tool_centrality.\
               central_point_dominance(g._Graph__graph,
                                       _prop("v", g, betweenness))

def eigentrust(g, trust_map, vprop=None, epslon=1e-6, max_iter=0,
               ret_iter=False):
    if vprop == None:
        vprop = g.new_vertex_property("double")
    libgraph_tool_centrality.\
            get_eigentrust(g._Graph__graph, _prop("e", g, trust_map),
                           _prop("v", g, vprop), epslon, max_iter)
    return vprop

def absolute_trust(g, trust_map, vprop=None, epslon=0.1, max_iter=0,
                   seed=0):
    if seed != 0:
        seed = numpy.random.randint(0, sys.maxint)
    if vprop == None:
        vprop = g.new_vertex_property("vector<double>")
    ic = libgraph_tool_centrality.\
            get_absolute_trust(g._Graph__graph, _prop("e", g, trust_map),
                               _prop("v", g, vprop), epslon, max_iter, seed)
    return vprop

