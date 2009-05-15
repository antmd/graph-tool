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
``graph_tool.draw`` - Graph Drawing
-----------------------------------
"""

import sys, os, os.path, time, warnings
from .. core import _degree, _prop, PropertyMap
from .. decorators import _limit_args

try:
    import gv
except ImportError:
    warnings.warn("error importing gv module... graph_draw() will not work.",
                  ImportWarning)
try:
    import matplotlib.cm
    import matplotlib.colors
except ImportError:
    warnings.warn("error importing matplotlib module... " + \
                  "graph_draw() will not work.", ImportWarning)

def graph_draw(g, pos=None, size=(15,15), pin=False, layout="neato",
               maxiter=None, ratio="fill", overlap=False, splines=False,
               mode="major", vsize=0.1, penwidth=1.0, eweight=None, ewidth=None,
               gprops={}, vprops={}, eprops={}, vcolor=None, ecolor=None,
               vcmap=matplotlib.cm.jet, vnorm=True, ecmap=matplotlib.cm.jet,
               enorm=True, output="", output_format="auto", returngv=False,
               fork=False, seed=0):
    """Draw a graph using graphviz."""

    if output != "":
        output = os.path.expanduser(output)
        # check opening file for writing, since graphview will bork if it is not
        # possible to open file
        if os.path.dirname(output) != "" and \
               not os.access(os.path.dirname(output), os.W_OK):
            raise IOError("cannot write to " + os.path.dirname(output))

    if g.is_directed():
        gvg = gv.digraph("G")
    else:
        gvg = gv.graph("G")

    # main graph properties
    gv.setv(gvg,"outputorder", "edgesfirst")
    gv.setv(gvg,"mode", mode)
    if overlap == False:
        if layout == "neato" and mode == "ipsep":
            overlap = "ipsep"
        else:
            overlap = "false"
    else:
        overlap = "true"
    if isinstance(overlap,str):
        gv.setv(gvg,"overlap", overlap)
    if splines:
        gv.setv(gvg,"splines", "true")
    gv.setv(gvg,"ratio", str(ratio))
    gv.setv(gvg,"size", "%f,%f" % (size[0]/2.54,size[1]/2.54)) # centimeters
    if maxiter != None:
        gv.setv(gvg,"maxiter", str(maxiter))
    if seed != 0:
        if type(seed) == int:
            gv.setv(gvg, "start", "%d" % seed)
        else:
            gv.setv(gvg, "start", seed)

    # apply all user supplied properties
    for k,val in gprops.iteritems():
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
        if minmax[0] == minmax[1]:
            minmax[1] += 1
        if vnorm:
            vnorm = matplotlib.colors.normalize(vmin=minmax[0], vmax=minmax[1])

    if ecolor != None and not isinstance(ecolor, str):
        minmax = [float("inf"), -float("inf")]
        for e in g.edges():
            c = ecolor[e]
            minmax[0] = min(c,minmax[0])
            minmax[1] = max(c,minmax[1])
        if minmax[0] == minmax[1]:
            minmax[1] += 1
        if enorm:
            enorm = matplotlib.colors.normalize(vmin=minmax[0], vmax=minmax[1])

    nodes = []
    edges = []

    # add nodes
    for v in g.vertices():
        n = gv.node(gvg,str(g.vertex_index[v]))
        if type(vsize) != tuple:
            vw = vh = vsize
        else:
            vw, vh = vsize
        if type(vw) == PropertyMap:
            vw = vw[v]
        if type(vh) == PropertyMap:
            vh = vh[v]

        if type(vw) == str and vw == "in":
            vw = v.in_degree()
        if type(vw) == str and vw == "out":
            vw = v.out_degree()
        if type(vw) == str and vw == "total":
            vw = v.in_degree() + v.out_degree()

        if type(vh) == str and vh == "in":
            vh = v.in_degree()
        if type(vh) == str and vh == "out":
            vh = v.out_degree()
        if type(vh) == str and vh == "total":
            vh = v.in_degree() + v.out_degree()

        gv.setv(n, "width", "%g" % vw)
        gv.setv(n, "height", "%g" % vh)
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
        for k,val in vprops.iteritems():
            if isinstance(val, PropertyMap):
                gv.setv(n, k, str(val[v]))
            else:
                gv.setv(n, k, str(val))
        nodes.append(n)
    for e in g.edges():
        ge = gv.edge(nodes[g.vertex_index[e.source()]],
                     nodes[g.vertex_index[e.target()]])
        gv.setv(ge, "arrowsize", "0.3")
        if g.is_directed():
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
        for k,v in eprops.iteritems():
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

    if output_format == "auto":
        if output == "":
            output_format = "xlib"
        else:
            output_format = output.split(".")[-1]

    # if using xlib we need to fork the process, otherwise good ol' graphviz
    # will call exit() when the window is closed
    if output_format == "xlib" or fork:
        pid = os.fork()
        if pid == 0:
            gv.render(gvg, output_format, output)
            os._exit(0) # since we forked, it's good to be sure
        if output_format != "xlib":
            os.wait()
    else:
        gv.render(gvg, output_format, output)

    if returngv:
        return pos, gv
    else:
        gv.rm(gvg)
        del gvg
        return pos
