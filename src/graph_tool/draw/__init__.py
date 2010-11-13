#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# graph_tool -- a general graph manipulation python module
#
# Copyright (C) 2007-2010 Tiago de Paula Peixoto <tiago@skewed.de>
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

Summary
+++++++

.. autosummary::
   :nosignatures:

   graph_draw
   arf_layout
   random_layout

Contents
++++++++
"""

import sys, os, os.path, time, warnings, tempfile
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
    from pylab import imread
except ImportError:
    warnings.warn("error importing matplotlib module... " + \
                  "graph_draw() will not work.", ImportWarning)

__all__ = ["graph_draw", "arf_layout", "random_layout"]


def graph_draw(g, pos=None, size=(15, 15), pin=False, layout=None,
               maxiter=None, ratio="fill", overlap="prism", sep=None,
               splines=False, vsize=0.1, penwidth=1.0, elen=None, gprops={},
               vprops={}, eprops={}, vcolor=None, ecolor=None,
               vcmap=matplotlib.cm.jet, vnorm=True, ecmap=matplotlib.cm.jet,
               enorm=True, vorder=None, eorder=None, output="",
               output_format="auto", returngv=False, fork=False,
               return_bitmap=False, seed=0):
    r"""Draw a graph using graphviz.

    Parameters
    ----------
    g : Graph
        Graph to be used.
    pos : PropertyMap or tuple of PropertyMaps (optional, default: None)
        Vertex property maps containing the x and y coordinates of the vertices.
    size : tuple of scalars (optional, default: (15,15))
        Size (in centimeters) of the canvas.
    pin : bool (default: False)
        If True, the vertices are not moved from their initial position.
    layout : string (default: "neato" if g.num_vertices() <= 1000 else "sfdp")
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
        drawn on a single page, then size is set to an "ideal" value. In
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
    vsize : float, PropertyMap, or tuple (default: 0.1)
        Default vertex size (width and height). If a tuple is specified, the
        first value should be a property map, and the second is a scale factor.
    penwidth : float, PropertyMap or tuple (default: 1.0)
        Specifies the width of the pen, in points, used to draw lines and
        curves, including the boundaries of edges and clusters. It has no effect
        on text. If a tuple is specified, the first value should be a property
        map, and the second is a scale factor.
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
    vorder : PropertyMap (default: None)
        Scalar vertex property map which specifies the order with which vertices
        are drawn.
    eorder : PropertyMap (default: None)
        Scalar edge property map which specifies the order with which edges
        are drawn.
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
        If True, the program is forked before drawing. This is used as a
        work-around for a bug in graphviz, where the exit() function is called,
        which would cause the calling program to end. This is always assumed
        'True', if output_format = 'xlib'.
    return_bitmap : bool (default: False)
        If True, a bitmap (:class:`~numpy.ndarray`) of the rendered graph is
        returned.

    Returns
    -------
    pos : PropertyMap
        Vector vertex property map with the x and y coordinates of the vertices.
    gv : gv.digraph or gv.graph (optional, only if returngv == True)
        Internally used graphviz graph.


    Notes
    -----
    This function is a wrapper for the [graphviz] python
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
    <...>

    .. figure:: graph-draw.png
        :align: center

        Kamada-Kawai force-directed layout of a graph with a power-law degree
        distribution, and dissortative degree correlation. The vertex size and
        color indicate the degree, and the edge color and width the edge
        betweeness centrality.

    References
    ----------
    .. [graphviz] http://www.graphviz.org

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

    if layout is None:
        layout = "neato" if g.num_vertices() <= 1000 else "sfdp"

    if layout == "arf":
        layout = "neato"
        pos = arf_layout(g, pos=pos)
        pin = True

    if pos != None:
        # copy user-supplied property
        if isinstance(pos, PropertyMap):
            pos = ungroup_vector_property(pos, [0,1])
        else:
            pos = (g.copy_property(pos[0]), g.copy_property(pos[1]))

    if type(vsize) == tuple:
        s = g.new_vertex_property("double")
        g.copy_property(vsize[0], s)
        s.a *= vsize[1]
        vsize = s

    if type(penwidth) == tuple:
        s = g.new_edge_property("double")
        g.copy_property(penwidth[0], s)
        s.a *= penwidth[1]
        penwidth = s

    # main graph properties
    gv.setv(gvg, "outputorder", "edgesfirst")
    gv.setv(gvg, "mode", "major")
    if overlap == False:
        overlap = "false"
    else:
        overlap = "true"
    if isinstance(overlap, str):
        gv.setv(gvg, "overlap", overlap)
    if sep != None:
        gv.setv(gvg, "sep", str(sep))
    if splines:
        gv.setv(gvg, "splines", "true")
    gv.setv(gvg, "ratio", str(ratio))
    # size is in centimeters... convert to inches
    gv.setv(gvg, "size", "%f,%f" % (size[0] / 2.54, size[1] / 2.54))
    if maxiter != None:
        gv.setv(gvg, "maxiter", str(maxiter))

    seed = numpy.random.randint(sys.maxint)
    gv.setv(gvg, "start", "%d" % seed)

    # apply all user supplied properties
    for k, val in gprops.iteritems():
        if isinstance(val, PropertyMap):
            gv.setv(gvg, k, str(val[g]))
        else:
            gv.setv(gvg, k, str(val))

    # normalize color properties
    if vcolor != None and not isinstance(vcolor, str):
        minmax = [float("inf"), -float("inf")]
        for v in g.vertices():
            c = vcolor[v]
            minmax[0] = min(c, minmax[0])
            minmax[1] = max(c, minmax[1])
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
            minmax[0] = min(c, minmax[0])
            minmax[1] = max(c, minmax[1])
        if minmax[0] == minmax[1]:
            minmax[1] += 1
        if enorm:
            enorm = matplotlib.colors.normalize(vmin=minmax[0], vmax=minmax[1])
        else:
            enorm = lambda x: x

    nodes = {}

    # add nodes
    if vorder != None:
        vertices = sorted(g.vertices(), lambda a, b: cmp(vorder[a], vorder[b]))
    else:
        vertices = g.vertices()
    for v in vertices:
        n = gv.node(gvg, str(g.vertex_index[v]))

        if type(vsize) == PropertyMap:
            vw = vh = vsize[v]
        else:
            vw = vh = vsize

        gv.setv(n, "width", "%g" % vw)
        gv.setv(n, "height", "%g" % vh)
        gv.setv(n, "style", "filled")
        gv.setv(n, "color", "black")
        # apply color
        if vcolor != None:
            if isinstance(vcolor, str):
                gv.setv(n, "fillcolor", vcolor)
            else:
                color = tuple([int(c * 255.0) for c in vcmap(vnorm(vcolor[v]))])
                gv.setv(n, "fillcolor", "#%.2x%.2x%.2x%.2x" % color)
        else:
            gv.setv(n, "fillcolor", "red")
        gv.setv(n, "label", "")

        # user supplied position
        if pos != None:
            gv.setv(n, "pos", "%f,%f" % (pos[0][v], pos[1][v]))
            gv.setv(n, "pin", str(pin))

        # apply all user supplied properties
        for k, val in vprops.iteritems():
            if isinstance(val, PropertyMap):
                gv.setv(n, k, str(val[v]))
            else:
                gv.setv(n, k, str(val))
        nodes[v] = n

    # add edges
    if eorder != None:
        edges = sorted(g.edges(), lambda a, b: cmp(eorder[a], eorder[b]))
    else:
        edges = g.edges()
    for e in edges:
        ge = gv.edge(nodes[e.source()],
                     nodes[e.target()])
        gv.setv(ge, "arrowsize", "0.3")
        if g.is_directed():
            gv.setv(ge, "arrowhead", "vee")

        # apply color
        if ecolor != None:
            if isinstance(ecolor, str):
                gv.setv(ge, "color", ecolor)
            else:
                color = tuple([int(c * 255.0) for c in ecmap(enorm(ecolor[e]))])
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
        for k, v in eprops.iteritems():
            if isinstance(v, PropertyMap):
                gv.setv(ge, k, str(v[e]))
            else:
                gv.setv(ge, k, str(v))

    gv.layout(gvg, layout)
    gv.render(gvg, "dot", "/dev/null")  # retrieve positions

    if pos == None:
        pos = (g.new_vertex_property("double"), g.new_vertex_property("double"))
    for n, n_gv in nodes.iteritems():
        p = gv.getv(n_gv, "pos")
        p = p.split(",")
        pos[0][n] = float(p[0])
        pos[1][n] = float(p[1])

    # I don't get this, but it seems necessary
    pos[0].a /= 100
    pos[1].a /= 100

    pos = group_vector_property(pos)

    if return_bitmap:
        # This is a not-so-nice hack which obtains an image buffer from a png
        # file. It is a pity that graphviz does not give access to its internal
        # buffers.
        tmp = tempfile.mkstemp(suffix=".png")[1]
        gv.render(gvg, "png", tmp)
        img = imread(tmp)
        os.remove(tmp)
    else:
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
                os._exit(0)  # since we forked, it's good to be sure
            if output_format != "xlib":
                os.wait()
        elif output != None:
            gv.render(gvg, output_format, output)

    ret = [pos]
    if return_bitmap:
        ret.append(img)

    if returngv:
        ret.append(gv)
    else:
        gv.rm(gvg)
        del gvg

    if len(ret) > 1:
        return tuple(ret)
    else:
        return ret[0]


def random_layout(g, shape=None, pos=None, dim=2):
    r"""Performs a random layout of the graph.

    Parameters
    ----------
    g : Graph
        Graph to be used.
    shape : tuple (optional, default: None)
        Rectangular shape of the bounding area. If None, a square of linear size
        :math:`\sqrt{N}` is used.
    pos : PropertyMap (optional, default: None)
        Vector vertex property maps where the coordinates should be stored.
    dim : int (optional, default: 2)
        Number of coordinates per vertex.

    Returns
    -------
    pos : A vector vertex property map
        Vertex property map with the coordinates of the vertices.

    Notes
    -----
    This algorithm has complexity :math:`O(V)`.
    """

    if pos == None:
        pos = [g.new_vertex_property("double") for i in xrange(dim)]

    if isinstance(pos, PropertyMap) and "vector" in pos.value_type():
        pos = ungroup_vector_property(pos)

    if shape == None:
        shape = [sqrt(g.num_vertices())] * dim

    for i in xrange(dim):
        _check_prop_scalar(pos[i], name="pos[%d]" % i)
        _check_prop_writable(pos[i], name="pos[%d]" % i)
        a = pos[i].get_array()
        a[:] = numpy.random.random(len(a)) * shape[i]

    pos = group_vector_property(g, pos)
    return pos


def arf_layout(g, weight=None, d=0.1, a=10, dt=0.001, epsilon=1e-6,
               max_iter=1000, pos=None, dim=2):
    r"""Calculate the ARF spring-block layout of the graph.

    Parameters
    ----------
    g : Graph
        Graph to be used.
    weight : PropertyMap (optional, default: None)
        An edge property map with the respective weights.
    d : float (optional, default: 0.1)
        Opposing force between vertices.
    a : float (optional, default: 10)
        Attracting force between adjacent vertices.
    dt : float (optional, default: 0.001)
        Iteration step size.
    epsilon : float (optional, default: 1e-6)
        Convergence criterion.
    max_iter : int (optional, default: 1000)
        Maximum number of iterations. If this value is 0, it runs until
        convergence.
    pos : PropertyMap (optional, default: None)
        Vector vertex property maps where the coordinates should be stored.
    dim : int (optional, default: 2)
        Number of coordinates per vertex.

    Returns
    -------
    pos : A vector vertex property map
        Vertex property map with the coordinates of the vertices.

    Notes
    -----
    This algorithm is defined in [geipel-self-organization-2007]_, and has
    complexity :math:`O(V^2)`.

    Examples
    --------
    >>> from numpy.random import seed, zipf
    >>> seed(42)
    >>> g = gt.random_graph(100, lambda: 3, directed=False)
    >>> t = gt.min_spanning_tree(g)
    >>> g.set_edge_filter(t)
    >>> pos = gt.graph_draw(g, output=None) # initial configuration
    >>> pos = gt.arf_layout(g, pos=pos, max_iter=0)
    >>> gt.graph_draw(g, pos=pos, pin=True, output="graph-draw-arf.png")
    <...>

    .. figure:: graph-draw-arf.png
        :align: center

        ARF layout of a minimum spanning tree of a random graph.

    References
    ----------
    .. [geipel-self-organization-2007] Markus M. Geipel, "Self-Organization
       applied to Dynamic Network Layout" , International Journal of Modern
       Physics C vol. 18, no. 10 (2007), pp. 1537-1549, arXiv:0704.1748v5
    .. _arf: http://www.sg.ethz.ch/research/graphlayout
    """

    if pos == None:
        if dim != 2:
            pos = random_layout(g, dim=dim)
        else:
            pos = graph_draw(g, output=None)
    _check_prop_vector(pos, name="pos", floating=True)

    g.stash_filter(directed=True)
    try:
        g.set_directed(False)
        libgraph_tool_layout.arf_layout(g._Graph__graph, _prop("v", g, pos),
                                        _prop("e", g, weight), d, a, dt,
                                        max_iter, epsilon, dim)
    finally:
        g.pop_filter(directed=True)
    return pos
