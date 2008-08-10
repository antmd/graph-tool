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

import sys,os
import matplotlib.cm
import matplotlib.colors
from .. core import _degree, _prop, PropertyMap
from .. decorators import _limit_args
import gv

def graph_draw(g, pos=None, size=(15,15), pin=False, layout="neato",
               maxiter=None, ratio="fill", overlap=False, splines=False,
               mode="ipsep", penwidth=1.0, eweight=None, ewidth=None, gprops={},
               vprops={}, eprops={}, vcolor=None, ecolor=None,
               vcmap=matplotlib.cm.jet, vnorm=True, ecmap=matplotlib.cm.jet,
               enorm=True, output=None, output_format=None):
    """Draw a graph using graphviz."""

    if g.is_directed():
        gvg = gv.digraph("G")
    else:
        gvg = gv.graph("G")

    # main graph properties
    gv.setv(gvg,"outputorder", "edgesfirst")
    gv.setv(gvg,"mode", mode)
    if overlap == "false" and layout == "neato" and mode == "ipsep":
        overlap = "ipsep"
    if overlap:
        gv.setv(gvg,"overlap", "true")
    elif isinstance(overlap,str):
        gv.setv(gvg,"overlap", overlap)
    if splines:
        gv.setv(gvg,"splines", "true")
    gv.setv(gvg,"ratio", str(ratio))
    gv.setv(gvg,"size", "%f,%f" % (size[0]/2.54,size[0]/2.54)) # centimeters
    if maxiter != None:
        gv.setv(gvg,"maxiter", str(maxiter))

    # apply all user supplied properties
    for k,val in vprops:
        if isinstance(val, PropertyMap):
            gv.setv(gvg, k, str(val[g]))
        else:
            gv.setv(gvg, k, str(val))

    # normalize color properties
    if vcolor != None and not isinstance(vcolor, str):
        minmax = [float("inf"), -float("inf")]
        for v in g.vertices():
            c = vcolor[v]
            minmax[0] = min(c,minmax[0])
            minmax[1] = max(c,minmax[1])
        print minmax
        if vnorm:
            vnorm = matplotlib.colors.normalize(vmin=minmax[0], vmax=minmax[1])

    if ecolor != None and not isinstance(ecolor, str):
        minmax = [float("inf"), -float("inf")]
        for e in g.edges():
            c = ecolor[v]
            minmax[0] = min(c,minmax[0])
            minmax[1] = max(c,minmax[1])
        if enorm:
            enorm = matplotlib.colors.normalize(vmin=minmax[0], vmax=minmax[1])

    nodes = []
    edges = []

    # add nodes
    for v in g.vertices():
        n = gv.node(gvg,str(g.vertex_index[v]))
        gv.setv(n, "width", "0.1")
        gv.setv(n, "height", "0.1")
        gv.setv(n, "height", "0.1")
        gv.setv(n, "style", "filled")
        gv.setv(n, "color", "black")
        # apply color
        if vcolor != None:
            if isinstance(vcolor,str):
                gv.setv(n, "fillcolor", vcolor)
            else:
                color = tuple([int(c*255.0) for c in vcmap(vnorm(vcolor[v]))])
                gv.setv(n, "fillcolor", "#%.2x%.2x%.2x%.2x" % color)
        else:
            gv.setv(n, "fillcolor", "red")
        gv.setv(n, "label", "")

        # user supplied position
        if pos != None:
            gv.setv(n, "pos", "%f,%f" % (pos[0][v],pos[1][v]))
            gv.setv(n, "pin", str(pin))

        # apply all user supplied properties
        for k,val in vprops:
            if isinstance(val, PropertyMap):
                gv.setv(n, k, str(val[v]))
            else:
                gv.setv(n, k, str(val))
        nodes.append(n)
    for e in g.edges():
        ge = gv.edge(nodes[g.vertex_index[e.source()]],
                     nodes[g.vertex_index[e.target()]])
        gv.setv(ge, "arrowsize", "0.3")
        gv.setv(ge, "arrowhead", "vee")
        # apply color
        if ecolor != None:
            if isinstance(ecolor,str):
                gv.setv(ge, "color", ecolor)
            else:
                color = tuple([int(c*255.0) for c in ecmap(enorm(ecolor[e]))])
                gv.setv(ge, "color", "#%.2x%.2x%.2x%.2x" % color)

        # apply weight
        if eweight != None:
            if isinstance(eweight, PropertyMap):
                gv.setv(ge, "weight", str(eweight[e]))
            else:
                gv.setv(ge, "weight", str(eweight))

        # apply width
        if ewidth != None:
            if isinstance(ewidth, PropertyMap):
                gv.setv(ge, "penwidth", str(ewidth[e]))
            else:
                gv.setv(ge, "penwidth", str(ewidth))

        # apply all user supplied properties
        for k,v in eprops:
            if isinstance(v, PropertyMap):
                gv.setv(ge, k, str(v[e]))
            else:
                gv.setv(ge, k, str(v))
        edges.append(ge)

    gv.layout(gvg, layout)
    gv.render(gvg, "dot", "/dev/null") # retrieve postitions

    if pos == None:
        pos = (g.new_vertex_property("double"), g.new_vertex_property("double"))
    for n in xrange(0, len(nodes)):
        p = gv.getv(nodes[n], "pos")
        p = p.split(",")
        pos[0][g.vertex(n)] = float(p[0])
        pos[1][g.vertex(n)] = float(p[1])
    if output == None:
        output = ""
        output_format = "xlib"
    elif output_format == None:
        outout_format = output.split(".")[-1]

    # if using xlib we need to fork the process, otherwise good ol' graphviz
    # will call exit() when the window is closed
    if output_format == "xlib":
        pid = os.fork()
        if pid == 0:
            gv.render(gvg, output_format, output)
            sys.exit(1) # since we forked, it's good to be sure
        os.wait()
    else:
        gv.render(gvg, output_format, output)
    gv.rm(gvg)
    return pos
