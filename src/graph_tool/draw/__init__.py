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
``graph_tool.draw`` - Graph drawing
-----------------------------------
"""

import sys, os, os.path, time, warnings
from .. core import _degree, _prop, PropertyMap, _check_prop_vector,\
     _check_prop_scalar, _check_prop_writable, group_vector_property,\
     ungroup_vector_property
from .. decorators import _limit_args
import numpy.random
from numpy import *

from .. dl_import import dl_import
dl_import("import libgraph_tool_layout")

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

__all__ = ["graph_draw", "arf_layout", "random_layout"]

def graph_draw(g, pos=None, size=(15, 15), pin=False, layout= "neato",
               maxiter=None, ratio= "fill", overlap="prism", sep=None,
               splines=False, vsize=0.1, penwidth=1.0, elen=None, gprops={},
               vprops={}, eprops={}, vcolor=None, ecolor=None,
               vcmap=matplotlib.cm.jet, vnorm=True, ecmap=matplotlib.cm.jet,
               enorm=True, output= "", output_format= "auto", returngv=False,
               fork=False, seed=0):
    r"""Draw a graph using graphviz.

    Parameters
    ----------
    g : Graph
        Graph to be used.
    pos : tuple of PropertyMaps (optional, default: None)
        Vertex property maps containing the x and y coordinates of the vertices.
    size : tuple of scalars (optional, default: (15,15))
        Size (in centimeters) of the canvas.
    pin : bool (default: False)
        If True, the vertices are not moved from their initial position.
    layout : string (default: "neato")
        Layout engine to be used. Possible values are "neato", "fdp", "dot",
        "circo", "twopi" and "arf".
    maxiter : int (default: None)
        If specified, limits the maximum number of iterations.
    ratio : string or float (default: "fill")
        Sets the aspect ratio (drawing height/drawing width) for the
        drawing. Note that this is adjusted before the 'size' attribute
        constraints are enforced.

        If ratio is numeric, it is taken as the desired aspect ratio. Then, if
        the actual aspect ratio is less than the desired ratio, the drawing
        height is scaled up to achieve the desired ratio; if the actual ratio is
        greater than that desired ratio, the drawing width is scaled up.

        If ratio = "fill" and the size attribute is set, node positions are
        scaled, separately in both x and y, so that the final drawing exactly
        fills the specified size.

        If ratio = "compress" and the size attribute is set, dot attempts to
        compress the initial layout to fit in the given size. This achieves a
        tighter packing of nodes but reduces the balance and symmetry.
        This feature only works in dot.

        If ratio = "expand", the size attribute is set, and both the width and
        the height of the graph are less than the value in size, node positions
        are scaled uniformly until at least one dimension fits size exactly.
        Note that this is distinct from using size as the desired size, as here
        the drawing is expanded before edges are generated and all node and text
        sizes remain unchanged.

        If ratio = "auto", the page attribute is set and the graph cannot be
        drawn on a single page, then size is set to an ``ideal'' value. In
        particular, the size in a given dimension will be the smallest integral
        multiple of the page size in that dimension which is at least half the
        current size. The two dimensions are then scaled independently to the
        new size. This feature only works in dot.
    overlap : bool or string (default: "prism")
        Determines if and how node overlaps should be removed. Nodes are first
        enlarged using the sep attribute. If True, overlaps are retained. If
        the value is "scale", overlaps are removed by uniformly scaling in x and
        y. If the value is False, node overlaps are removed by a Voronoi-based
        technique. If the value is "scalexy", x and y are separately scaled to
        remove overlaps.

        If sfdp is available, one can set overlap to "prism" to use a proximity
        graph-based algorithm for overlap removal. This is the preferred
        technique, though "scale" and False can work well with small graphs.
        This technique starts with a small scaling up, controlled by the
        overlap_scaling attribute, which can remove a significant portion of the
        overlap. The prism option also accepts an optional non-negative integer
        suffix. This can be used to control the number of attempts made at
        overlap removal. By default, overlap="prism" is equivalent to
        overlap="prism1000". Setting overlap="prism0" causes only the scaling
        phase to be run.

        If the value is "compress", the layout will be scaled down as much as
        possible without introducing any overlaps, obviously assuming there are
        none to begin with.
    sep : float (default: None)
        Specifies margin to leave around nodes when removing node overlap. This
        guarantees a minimal non-zero distance between nodes.
    splines : bool (default: False)
        If True, the edges are drawn as splines and routed around the vertices.
    vsize : float or PropertyMap (default: 0.1)
        Default vertex size (width and height).
    penwidth : float  or PropertyMap (default: 1.0)
        Specifies the width of the pen, in points, used to draw lines and
        curves, including the boundaries of edges and clusters. It has no effect
        on text.
    elen : float or PropertyMap (default: None)
        Preferred edge length, in inches.
    gprops : dict (default: {})
        Additional graph properties, as a dictionary. The keys are the property
        names, and the values must be convertible to string.
    vprops : dict (default: {})
        Additional vertex properties, as a dictionary. The keys are the property
        names, and the values must be convertible to string, or vertex property
        maps, with values convertible to strings.
    eprops : dict (default: {})
        Additional edge properties, as a dictionary. The keys are the property
        names, and the values must be convertible to string, or edge property
        maps, with values convertible to strings.
    vcolor : string or PropertyMap (default: None)
        Drawing color for vertices. If the valued supplied is a property map,
        the values must be scalar types, whose color values are obtained from
        the 'vcmap' argument.
    ecolor : string or PropertyMap (default: None)
        Drawing color for edges. If the valued supplied is a property map,
        the values must be scalar types, whose color values are obtained from
        the 'ecmap' argument.
    vcmap : matplotlib.colors.Colormap (default: matplotlib.cm.jet)
        Vertex color map.
    vnorm : bool (default: True)
        Normalize vertex color values to the [0,1] range.
    ecmap : matplotlib.colors.Colormap (default: matplotlib.cm.jet)
        Edge color map.
    enorm : bool (default: True)
        Normalize edge color values to the [0,1] range.
    output : string (default: "")
        Output file name.
    output_format : string (default: "auto")
        Output file format. Possible values are "auto", "xlib", "ps", "svg",
        "svgz", "fig", "mif", "hpgl", "pcl", "png", "gif", "dia", "imap",
        "cmapx". If the value is "auto", the format is guessed from the 'output'
        parameter, or 'xlib' if it is empty. If the value is None, no output is
        produced.
    returngv : bool (default: False)
        Return the graph object used internally with the gv module.
    fork : bool (default: False)
        If true, the program is forked before drawing. This is used as a
        work-around for a bug in graphviz, where the exit() function is called,
        which would cause the calling program to end. This is always assumed
        'True', if output_format = 'xlib'.
    seed : int (default: 0)
        Seed for the random number generator. If the value 0, a different random
        value is used each time.

    Returns
    -------
    pos_x : PropertyMap
        Vertex property map with the x-coordinates of the vertices.
    pos_y : PropertyMap
        Vertex property map with the y-coordinates of the vertices.
    gv : gv.digraph or gv.graph (optional, only if returngv == True)
        Internally used graphviz graph.


    Notes
    -----
    This function is a wrapper for the graphviz_ python
    routines. Extensive additional documentation for the graph, vertex and edge
    properties is available at: http://www.graphviz.org/doc/info/attrs.html.


    Examples
    --------
    >>> from numpy import *
    >>> from numpy.random import seed, zipf
    >>> seed(42)
    >>> g = gt.random_graph(1000, lambda: min(zipf(2.4), 40),
    ...                     lambda i,j: exp(abs(i-j)), directed=False)
    >>> # extract largest component
    >>> comp = gt.label_components(g)
    >>> h = gt.vertex_hist(g, comp)
    >>> max_comp = h[1][list(h[0]).index(max(h[0]))]
    >>> g.remove_vertex_if(lambda v: comp[v] != max_comp)
    >>>
    >>> deg = g.degree_property_map("out")
    >>> deg.get_array()[:] = 2*(sqrt(deg.get_array()[:])*0.5 + 0.4)
    >>> ebet = gt.betweenness(g)[1]
    >>> ebet.get_array()[:] *= 4000
    >>> ebet.get_array()[:] += 10
    >>> gt.graph_draw(g, vsize=deg, vcolor=deg, elen=10, ecolor=ebet,
    ...               penwidth=ebet, overlap="prism", output="graph-draw.png")
    (...)

    .. figure:: graph-draw.png
        :align: center

        Kamada-Kawai force-directed layout of a graph with a power-law degree
        distribution, and dissortative degree correlation. The vertex size and
        color indicate the degree, and the edge color and width the edge
        betweeness centrality.

    References
    ----------
    .. _graphviz: http://www.graphviz.org


    """

    if output != "" and output != None:
        output = os.path.expanduser(output)
        # check opening file for writing, since graphviz will bork if it is not
        # possible to open file
        if os.path.dirname(output) != "" and \
               not os.access(os.path.dirname(output), os.W_OK):
            raise IOError("cannot write to " + os.path.dirname(output))

    if g.is_directed():
        gvg = gv.digraph("G")
    else:
        gvg = gv.graph("G")

    if layout == "arf":
        layout = "neato"
        pos = arf_layout(g, pos=pos)
        pin = True

    if pos != None:
        # copy user-supplied property
        if isinstance(pos, PropertyMap):
            pos = ungroup_vector_property(g, pos, [0,1])
        else:
            pos = (g.copy_property(pos[0]), g.copy_property(pos[1]))

    # main graph properties
    gv.setv(gvg,"outputorder", "edgesfirst")
    gv.setv(gvg,"mode", "major")
    if overlap == False:
        overlap = "false"
    else:
        overlap = "true"
    if isinstance(overlap,str):
        gv.setv(gvg,"overlap", overlap)
    if sep != None:
        gv.setv(gvg,"sep", str(sep))
    if splines:
        gv.setv(gvg,"splines", "true")
    gv.setv(gvg,"ratio", str(ratio))
    gv.setv(gvg,"size", "%f,%f" % (size[0]/2.54,size[1]/2.54)) # centimeters
    if maxiter != None:
        gv.setv(gvg,"maxiter", str(maxiter))

    if seed == 0:
        seed = numpy.random.randint(sys.maxint)
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
        else:
            vnorm = lambda x: x

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
        else:
            enorm = lambda x: x

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

        # apply edge length
        if elen != None:
            if isinstance(elen, PropertyMap):
                gv.setv(ge, "len", str(elen[e]))
            else:
                gv.setv(ge, "len", str(elen))

        # apply width
        if penwidth != None:
            if isinstance(penwidth, PropertyMap):
                gv.setv(ge, "penwidth", str(penwidth[e]))
            else:
                gv.setv(ge, "penwidth", str(penwidth))

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
        elif output != None:
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
    elif output != None:
        gv.render(gvg, output_format, output)

    # I don't get this, but it seems necessary
    pos[0].get_array()[:] /= 100
    pos[1].get_array()[:] /= 100

    pos = group_vector_property(g, pos)

    if returngv:
        return pos, gv
    else:
        gv.rm(gvg)
        del gvg
        return pos

def random_layout(g, shape=None, pos=None, dim=2):
    if pos == None:
        pos = [g.new_vertex_property("double") for i in xrange(dim)]

    if isinstance(pos, PropertyMap) and "vector" in pos.value_type():
        pos = ungroup_vector_property(pos)

    if shape == None:
        shape = (sqrt(g.num_vertices()), sqrt(g.num_vertices()))

    for i in xrange(dim):
        _check_prop_scalar(pos[i], name="pos[%d]" % i)
        _check_prop_writable(pos[i], name="pos[%d]" % i)
        a = pos[i].get_array()
        a[:] = numpy.random.random(len(a))*shape[i]

    pos = group_vector_property(g, pos)
    return pos

def arf_layout(g, weight=None, d=0.1, a=10, dt=0.001, epsilon=1e-6,
               max_iter=1000, pos=None, dim=2):
    if pos == None:
        pos = random_layout(g, dim=dim)
    _check_prop_vector(pos, name="pos", floating=True)

    g.stash_filter(directed=True)
    g.set_directed(False)
    libgraph_tool_layout.arf_layout(g._Graph__graph, _prop("v", g, pos),
                                    _prop("e", g, weight), d, a, dt, max_iter,
                                    epsilon, dim)
    g.pop_filter()
    return pos
