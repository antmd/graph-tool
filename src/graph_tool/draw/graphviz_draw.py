#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# graph_tool -- a general graph manipulation python module
#
# Copyright (C) 2006-2015 Tiago de Paula Peixoto <tiago@skewed.de>
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

from __future__ import division, absolute_import, print_function

import sys
import os
import os.path
import time
import warnings
import ctypes
import ctypes.util
import tempfile
from .. import PropertyMap, group_vector_property, ungroup_vector_property
import numpy.random
import copy

from .. draw import arf_layout

try:
    import matplotlib.cm
    import matplotlib.colors
except ImportError:
    msg = "Error importing matplotlib module... graphviz_draw() will not work"
    warnings.filterwarnings("always", msg, ImportWarning)
    warnings.warn(msg, ImportWarning)
    raise

try:
    libname = ctypes.util.find_library("c")
    libc = ctypes.CDLL(libname)
    if hasattr(libc, "open_memstream"):
        libc.open_memstream.restype = ctypes.POINTER(ctypes.c_char)
except OSError:
    msg = "Error importing C standard library... graphviz_draw() will not work."
    warnings.filterwarnings("always", msg, ImportWarning)
    warnings.warn(msg, ImportWarning)
    pass


try:
    libname = ctypes.util.find_library("gvc")
    if libname is None:
        raise OSError()
    libgv = ctypes.CDLL(libname, ctypes.RTLD_GLOBAL)
    # properly set the return types of certain functions
    ptype = ctypes.POINTER(ctypes.c_char)
    libgv.gvContext.restype = ptype
    libgv.agopen.restype = ptype
    libgv.agnode.restype = ptype
    libgv.agedge.restype = ptype
    libgv.agget.restype = ptype
    libgv.agstrdup_html.restype = ptype
    # create a context to use the whole time (if we keep freeing and recreating
    # it, we will hit a memory leak in graphviz)
    gvc = libgv.gvContext()

    try:
        gv_new_api = True
        libgv_directed = libgv.Agdirected
        libgv_undirected = libgv.Agundirected
    except AttributeError:
        gv_new_api = False
        libgv_directed = 1
        libgv_undirected = 0

except OSError:
    msg = "Error importing graphviz C library (libgvc)... graphviz_draw() will not work."
    warnings.filterwarnings("always", msg, ImportWarning)
    warnings.warn(msg, ImportWarning)


def htmlize(val):
    if len(val) >= 2 and val[0] == "<" and val[-1] == ">":
        return ctypes.string_at(libgv.agstrdup_html(val[1:-1]))
    return val


def aset(elem, attr, value):
    v = htmlize(str(value)).encode("utf8")
    libgv.agsafeset(elem, str(attr).encode("utf8"), v, v)


def aget(elem, attr):
    s = ctypes.string_at(libgv.agget(elem,
                                     str(attr).encode("utf8")))
    return s.decode("utf8")


def graphviz_draw(g, pos=None, size=(15, 15), pin=False, layout=None,
                  maxiter=None, ratio="fill", overlap=True, sep=None,
                  splines=False, vsize=0.105, penwidth=1.0, elen=None,
                  gprops={}, vprops={}, eprops={}, vcolor="#a40000",
                  ecolor="#2e3436", vcmap=None, vnorm=True, ecmap=None,
                  enorm=True, vorder=None, eorder=None, output="",
                  output_format="auto", fork=False, return_string=False):
    r"""Draw a graph using graphviz.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be drawn.
    pos : :class:`~graph_tool.PropertyMap` or tuple of :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Vertex property maps containing the x and y coordinates of the vertices.
    size : tuple of scalars (optional, default: ``(15,15)``)
        Size (in centimeters) of the canvas.
    pin : bool or :class:`~graph_tool.PropertyMap` (default: ``False``)
        If ``True``, the vertices are not moved from their initial position. If
        a :class:`~graph_tool.PropertyMap` is passed, it is used to pin nodes
        individually.
    layout : string (default: ``"neato" if g.num_vertices() <= 1000 else "sfdp"``)
        Layout engine to be used. Possible values are ``"neato"``, ``"fdp"``,
        ``"dot"``, ``"circo"``, ``"twopi"`` and ``"arf"``.
    maxiter : int (default: ``None``)
        If specified, limits the maximum number of iterations.
    ratio : string or float (default: ``"fill"``)
        Sets the aspect ratio (drawing height/drawing width) for the
        drawing. Note that this is adjusted before the ``size`` attribute
        constraints are enforced.

        If ``ratio`` is numeric, it is taken as the desired aspect ratio. Then,
        if the actual aspect ratio is less than the desired ratio, the drawing
        height is scaled up to achieve the desired ratio; if the actual ratio is
        greater than that desired ratio, the drawing width is scaled up.

        If ``ratio == "fill"`` and the size attribute is set, node positions are
        scaled, separately in both x and y, so that the final drawing exactly
        fills the specified size.

        If ``ratio == "compress"`` and the size attribute is set, dot attempts
        to compress the initial layout to fit in the given size. This achieves a
        tighter packing of nodes but reduces the balance and symmetry.  This
        feature only works in dot.

        If ``ratio == "expand"``, the size attribute is set, and both the width
        and the height of the graph are less than the value in size, node
        positions are scaled uniformly until at least one dimension fits size
        exactly.  Note that this is distinct from using size as the desired
        size, as here the drawing is expanded before edges are generated and all
        node and text sizes remain unchanged.

        If ``ratio == "auto"``, the page attribute is set and the graph cannot
        be drawn on a single page, then size is set to an "ideal" value. In
        particular, the size in a given dimension will be the smallest integral
        multiple of the page size in that dimension which is at least half the
        current size. The two dimensions are then scaled independently to the
        new size. This feature only works in dot.
    overlap : bool or string (default: ``"prism"``)
        Determines if and how node overlaps should be removed. Nodes are first
        enlarged using the sep attribute. If ``True``, overlaps are retained. If
        the value is ``"scale"``, overlaps are removed by uniformly scaling in x
        and y. If the value is ``False``, node overlaps are removed by a
        Voronoi-based technique. If the value is ``"scalexy"``, x and y are
        separately scaled to remove overlaps.

        If sfdp is available, one can set overlap to ``"prism"`` to use a
        proximity graph-based algorithm for overlap removal. This is the
        preferred technique, though ``"scale"`` and ``False`` can work well with
        small graphs. This technique starts with a small scaling up, controlled
        by the overlap_scaling attribute, which can remove a significant portion
        of the overlap. The prism option also accepts an optional non-negative
        integer suffix. This can be used to control the number of attempts made
        at overlap removal. By default, ``overlap == "prism"`` is equivalent to
        ``overlap == "prism1000"``. Setting ``overlap == "prism0"`` causes only
        the scaling phase to be run.

        If the value is ``"compress"``, the layout will be scaled down as much
        as possible without introducing any overlaps, obviously assuming there
        are none to begin with.
    sep : float (default: ``None``)
        Specifies margin to leave around nodes when removing node overlap. This
        guarantees a minimal non-zero distance between nodes.
    splines : bool (default: ``False``)
        If ``True``, the edges are drawn as splines and routed around the
        vertices.
    vsize : float, :class:`~graph_tool.PropertyMap`, or tuple (default: ``0.105``)
        Default vertex size (width and height). If a tuple is specified, the
        first value should be a property map, and the second is a scale factor.
    penwidth : float, :class:`~graph_tool.PropertyMap` or tuple (default: ``1.0``)
        Specifies the width of the pen, in points, used to draw lines and
        curves, including the boundaries of edges and clusters. It has no effect
        on text. If a tuple is specified, the first value should be a property
        map, and the second is a scale factor.
    elen : float or :class:`~graph_tool.PropertyMap` (default: ``None``)
        Preferred edge length, in inches.
    gprops : dict (default: ``{}``)
        Additional graph properties, as a dictionary. The keys are the property
        names, and the values must be convertible to string.
    vprops : dict (default: ``{}``)
        Additional vertex properties, as a dictionary. The keys are the property
        names, and the values must be convertible to string, or vertex property
        maps, with values convertible to strings.
    eprops : dict (default: ``{}``)
        Additional edge properties, as a dictionary. The keys are the property
        names, and the values must be convertible to string, or edge property
        maps, with values convertible to strings.
    vcolor : string or :class:`~graph_tool.PropertyMap` (default: ``"#a40000"``)
        Drawing color for vertices. If the valued supplied is a property map,
        the values must be scalar types, whose color values are obtained from
        the ``vcmap`` argument.
    ecolor : string or :class:`~graph_tool.PropertyMap` (default: ``"#2e3436"``)
        Drawing color for edges. If the valued supplied is a property map,
        the values must be scalar types, whose color values are obtained from
        the ``ecmap`` argument.
    vcmap : :class:`matplotlib.colors.Colormap` (default: :class:`matplotlib.cm.jet`)
        Vertex color map.
    vnorm : bool (default: ``True``)
        Normalize vertex color values to the [0,1] range.
    ecmap : :class:`matplotlib.colors.Colormap` (default: :class:`matplotlib.cm.jet`)
        Edge color map.
    enorm : bool (default: ``True``)
        Normalize edge color values to the [0,1] range.
    vorder : :class:`~graph_tool.PropertyMap` (default: ``None``)
        Scalar vertex property map which specifies the order with which vertices
        are drawn.
    eorder : :class:`~graph_tool.PropertyMap` (default: ``None``)
        Scalar edge property map which specifies the order with which edges
        are drawn.
    output : string (default: ``""``)
        Output file name.
    output_format : string (default: ``"auto"``)
        Output file format. Possible values are ``"auto"``, ``"xlib"``,
        ``"ps"``, ``"svg"``, ``"svgz"``, ``"fig"``, ``"mif"``, ``"hpgl"``,
        ``"pcl"``, ``"png"``, ``"gif"``, ``"dia"``, ``"imap"``, ``"cmapx"``. If
        the value is ``"auto"``, the format is guessed from the ``output``
        parameter, or ``xlib`` if it is empty. If the value is ``None``, no
        output is produced.
    fork : bool (default: ``False``)
        If ``True``, the program is forked before drawing. This is used as a
        work-around for a bug in graphviz, where the ``exit()`` function is
        called, which would cause the calling program to end. This is always
        assumed ``True``, if ``output_format == 'xlib'``.
    return_string : bool (default: ``False``)
        If ``True``, a string containing the rendered graph as binary data is
        returned (defaults to png format).

    Returns
    -------
    pos : :class:`~graph_tool.PropertyMap`
        Vector vertex property map with the x and y coordinates of the vertices.
    gv : gv.digraph or gv.graph (optional, only if ``returngv == True``)
        Internally used graphviz graph.


    Notes
    -----
    This function is a wrapper for the [graphviz]_ routines. Extensive additional
    documentation for the graph, vertex and edge properties is available at:
    http://www.graphviz.org/doc/info/attrs.html.


    Examples
    --------
    .. testcode::
       :hide:

       np.random.seed(42)
       gt.seed_rng(42)
       from numpy import sqrt

    >>> g = gt.price_network(1500)
    >>> deg = g.degree_property_map("in")
    >>> deg.a = 2 * (sqrt(deg.a) * 0.5 + 0.4)
    >>> ebet = gt.betweenness(g)[1]
    >>> gt.graphviz_draw(g, vcolor=deg, vorder=deg, elen=10,
    ...                  ecolor=ebet, eorder=ebet, output="graphviz-draw.pdf")
    <...>

    .. testcode::
       :hide:

       gt.graphviz_draw(g, vcolor=deg, vorder=deg, elen=10,
                        ecolor=ebet, eorder=ebet, output="graphviz-draw.png")

    .. figure:: graphviz-draw.*
        :align: center

        Kamada-Kawai force-directed layout of a Price network with 1500
        nodes. The vertex size and color indicate the degree, and the edge color
        corresponds to the edge betweeness centrality

    References
    ----------
    .. [graphviz] http://www.graphviz.org

    """

    if output != "" and output is not None:
        output = os.path.expanduser(output)
        # check opening file for writing, since graphviz will bork if it is not
        # possible to open file
        if os.path.dirname(output) != "" and \
               not os.access(os.path.dirname(output), os.W_OK):
            raise IOError("cannot write to " + os.path.dirname(output))

    has_layout = False
    try:
        if gv_new_api:
            gvg = libgv.agopen("G".encode("utf8"),
                               libgv_directed if g.is_directed() else libgv_undirected,
                               None)
        else:
            gvg = libgv.agopen("G".encode("utf8"),
                               libgv_directed if g.is_directed() else libgv_undirected)

        if layout is None:
            if pin == False:
                layout = "neato" if g.num_vertices() <= 1000 else "sfdp"
            else:
                layout = "neato"

        if layout == "arf":
            layout = "neato"
            pos = arf_layout(g, pos=pos)
            pin = True

        if pos is not None:
            # copy user-supplied property
            if isinstance(pos, PropertyMap):
                pos = ungroup_vector_property(pos, [0, 1])
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
        aset(gvg, "outputorder", "edgesfirst")
        aset(gvg, "mode", "major")
        if type(overlap) is bool:
            overlap = "true" if overlap else "false"
        else:
            overlap = str(overlap)
        aset(gvg, "overlap", overlap)
        if sep is not None:
            aset(gvg, "sep", sep)
        if splines:
            aset(gvg, "splines", "true")
        aset(gvg, "ratio", ratio)
        # size is in centimeters... convert to inches
        aset(gvg, "size", "%f,%f" % (size[0] / 2.54, size[1] / 2.54))
        if maxiter is not None:
            aset(gvg, "maxiter", maxiter)

        seed = numpy.random.randint(sys.maxsize)
        aset(gvg, "start", "%d" % seed)

        # apply all user supplied graph properties
        for k, val in gprops.items():
            if isinstance(val, PropertyMap):
                aset(gvg, k, val[g])
            else:
                aset(gvg, k, val)

        # normalize color properties
        if (isinstance(vcolor, PropertyMap) and
            vcolor.value_type() != "string"):
            minmax = [float("inf"), -float("inf")]
            for v in g.vertices():
                c = vcolor[v]
                minmax[0] = min(c, minmax[0])
                minmax[1] = max(c, minmax[1])
            if minmax[0] == minmax[1]:
                minmax[1] += 1
            if vnorm:
                vnorm = matplotlib.colors.Normalize(vmin=minmax[0], vmax=minmax[1])
            else:
                vnorm = lambda x: x

        if (isinstance(ecolor, PropertyMap) and
            ecolor.value_type() != "string"):
            minmax = [float("inf"), -float("inf")]
            for e in g.edges():
                c = ecolor[e]
                minmax[0] = min(c, minmax[0])
                minmax[1] = max(c, minmax[1])
            if minmax[0] == minmax[1]:
                minmax[1] += 1
            if enorm:
                enorm = matplotlib.colors.Normalize(vmin=minmax[0],
                                                    vmax=minmax[1])
            else:
                enorm = lambda x: x

        if vcmap is None:
            vcmap = matplotlib.cm.jet

        if ecmap is None:
            ecmap = matplotlib.cm.jet

        # add nodes
        if vorder is not None:
            vertices = sorted(g.vertices(), key = lambda a: vorder[a])
        else:
            vertices = g.vertices()
        for v in vertices:
            if gv_new_api:
                n = libgv.agnode(gvg, str(int(v)).encode("utf8"))
            else:
                n = libgv.agnode(gvg, str(int(v)).encode("utf8"), True)

            if type(vsize) == PropertyMap:
                vw = vh = vsize[v]
            else:
                vw = vh = vsize

            aset(n, "shape", "circle")
            aset(n, "width", "%g" % vw)
            aset(n, "height", "%g" % vh)
            aset(n, "style", "filled")
            aset(n, "color", "#2e3436")
            # apply color
            if isinstance(vcolor, str):
                aset(n, "fillcolor", vcolor)
            else:
                color = vcolor[v]
                if isinstance(color, str):
                    aset(n, "fillcolor", color)
                else:
                    color = tuple([int(c * 255.0) for c in vcmap(vnorm(color))])
                    aset(n, "fillcolor", "#%.2x%.2x%.2x%.2x" % color)
            aset(n, "label", "")

            # user supplied position
            if pos is not None:
                if isinstance(pin, bool):
                    pin_val = pin
                else:
                    pin_val = pin[v]
                aset(n, "pos", "%f,%f%s" % (pos[0][v], pos[1][v],
                                            "!" if pin_val else ""))
                aset(n, "pin", pin_val)

            # apply all user supplied properties
            for k, val in vprops.items():
                if isinstance(val, PropertyMap):
                    aset(n, k, val[v])
                else:
                    aset(n, k, val)

        # add edges
        if eorder is not None:
            edges = sorted(g.edges(), key = lambda a: eorder[a])
        else:
            edges = g.edges()
        for e in edges:
            if gv_new_api:
                ge = libgv.agedge(gvg,
                                  libgv.agnode(gvg, str(int(e.source())).encode("utf8"), False),
                                  libgv.agnode(gvg, str(int(e.target())).encode("utf8"), False),
                                  str(g.edge_index[e]).encode("utf8"), True)
            else:
                ge = libgv.agedge(gvg,
                                  libgv.agnode(gvg, str(int(e.source())).encode("utf8")),
                                  libgv.agnode(gvg, str(int(e.target())).encode("utf8")))
            aset(ge, "arrowsize", "0.3")
            if g.is_directed():
                aset(ge, "arrowhead", "vee")

            # apply color
            if isinstance(ecolor, str):
                aset(ge, "color", ecolor)
            else:
                color = ecolor[e]
                if isinstance(color, str):
                    aset(ge, "color", color)
                else:
                    color = tuple([int(c * 255.0) for c in ecmap(enorm(color))])
                    aset(ge, "color", "#%.2x%.2x%.2x%.2x" % color)

            # apply edge length
            if elen is not None:
                if isinstance(elen, PropertyMap):
                    aset(ge, "len", elen[e])
                else:
                    aset(ge, "len", elen)

            # apply width
            if penwidth is not None:
                if isinstance(penwidth, PropertyMap):
                    aset(ge, "penwidth", penwidth[e])
                else:
                    aset(ge, "penwidth", penwidth)

            # apply all user supplied properties
            for k, v in eprops.items():
                if isinstance(v, PropertyMap):
                    aset(ge, k, v[e])
                else:
                    aset(ge, k, v)

        libgv.gvLayout(gvc, gvg, layout.encode("utf8"))
        has_layout = True
        retv = libgv.gvRender(gvc, gvg, "dot".encode("utf8"), None)  # retrieve positions only

        if pos == None:
            pos = (g.new_vertex_property("double"),
                   g.new_vertex_property("double"))
        for v in g.vertices():
            n = libgv.agnode(gvg, str(int(v)).encode("utf8"))
            p = aget(n, "pos")
            p = p.split(",")
            pos[0][v] = float(p[0])
            pos[1][v] = float(p[1])

        # I don't get this, but it seems necessary
        pos[0].a /= 100
        pos[1].a /= 100

        pos = group_vector_property(pos)

        if return_string:
            if output_format == "auto":
                output_format = "png"
            if hasattr(libc, "open_memstream"):
                buf = ctypes.c_char_p()
                buf_len = ctypes.c_size_t()
                fstream = libc.open_memstream(ctypes.byref(buf),
                                              ctypes.byref(buf_len))
                libgv.gvRender(gvc, gvg, output_format.encode("utf8"), fstream)
                libc.fclose(fstream)
                data = copy.copy(ctypes.string_at(buf, buf_len.value))
                libc.free(buf)
            else:
                # write to temporary file, if open_memstream is not available
                output = tempfile.mkstemp()[1]
                libgv.gvRenderFilename(gvc, gvg, output_format.encode("utf8"),
                                       output.encode("utf8"))
                data = open(output).read()
                os.remove(output)
        else:
            if output_format == "auto":
                if output == "":
                    output_format = "xlib"
                elif output is not None:
                    output_format = output.split(".")[-1]

            # if using xlib we need to fork the process, otherwise good ol'
            # graphviz will call exit() when the window is closed
            if output_format == "xlib" or fork:
                pid = os.fork()
                if pid == 0:
                    libgv.gvRenderFilename(gvc, gvg, output_format.encode("utf8"),
                                           output.encode("utf8"))
                    os._exit(0)  # since we forked, it's good to be sure
                if output_format != "xlib":
                    os.wait()
            elif output is not None:
                libgv.gvRenderFilename(gvc, gvg, output_format.encode("utf8"),
                                       output.encode("utf8"))

        ret = [pos]
        if return_string:
            ret.append(data)

    finally:
        if has_layout:
            libgv.gvFreeLayout(gvc, gvg)
        libgv.agclose(gvg)

    if len(ret) > 1:
        return tuple(ret)
    else:
        return ret[0]
