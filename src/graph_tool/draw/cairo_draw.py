#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# graph_tool -- a general graph manipulation python module
#
# Copyright (C) 2007-2012 Tiago de Paula Peixoto <tiago@skewed.de>
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

import os
import warnings

try:
    import cairo
except ImportError:
    msg = "Error importing cairo. Graph drawing will not work."
    warnings.filterwarnings("always", msg, ImportWarning)
    warnings.warn(msg, ImportWarning)

try:
    import matplotlib.cm
    import matplotlib.colors
except ImportError:
    msg = "Error importing matplotlib module. Graph drawing will not work."
    warnings.filterwarnings("always", msg, ImportWarning)
    warnings.warn(msg, ImportWarning)

import numpy as np
import gzip
import bz2
import zipfile
import copy
from collections import defaultdict

from .. import GraphView, PropertyMap, ungroup_vector_property,\
     group_vector_property, _prop
from .. stats import label_parallel_edges, label_self_loops

from .. dl_import import dl_import
dl_import("import libgraph_tool_draw")
try:
    from libgraph_tool_draw import vertex_attrs, edge_attrs, vertex_shape,\
        edge_marker
except ImportError:
    msg = "Error importing cairo-based drawing library. " + \
        "Was graph-tool compiled with cairomm support?"
    warnings.filterwarnings("always", msg, ImportWarning)
    warnings.warn(msg, ImportWarning)

from .. draw import sfdp_layout, random_layout, _avg_edge_distance, \
    coarse_graphs

_vdefaults = {
    "shape": "circle",
    "color": [0, 0, 0, 1],
    "fill_color": [0.640625, 0, 0, 0.9],
    "size": 5,
    "pen_width": 0.8,
    "halo": 0,
    "halo_color": [0., 0., 1., 0.5],
    "text": "",
    "text_color": [0., 0., 0., 1.],
    "text_position": -1.,
    "font_family": "serif",
    "font_slant": cairo.FONT_SLANT_NORMAL,
    "font_weight": cairo.FONT_WEIGHT_NORMAL,
    "font_size": 12.
    }

_edefaults = {
    "color": [0.1796875, 0.203125, 0.2109375, 0.8],
    "pen_width": 1,
    "start_marker": "none",
    "mid_marker": "none",
    "end_marker": "none",
    "marker_size": 4.,
    "control_points": [],
    }


def shape_from_prop(shape, enum):
    if isinstance(shape, PropertyMap):
        if shape.key_type() == "v":
            prop = shape.get_graph().new_vertex_property("int")
            descs = shape.get_graph().vertices()
        else:
            descs = shape.get_graph().edges()
            prop = shape.get_graph().new_edge_property("int")
        offset = min(enum.values.keys())
        vals = dict([(k - offset, v) for k, v in enum.values.items()])
        for v in descs:
            if shape.value_type() == "string":
                prop[v] = int(enum.__dict__[shape[v]])
            elif shape[v] in vals:
                prop[v] = int(vals[shape[v]])
            else:
                prop[v] = int(vals[hash(shape[v]) % len(vals)])
        return prop

    if isinstance(shape, str):
        return int(enum.__dict__[shape])
    else:
        return shape

    raise ValueError("Invalid value for attribute %s: %s" %
                     (repr(enum), repr(shape)))


def _convert(attr, val, cmap):
    if attr == vertex_attrs.shape:
        return shape_from_prop(val, vertex_shape)
    if attr in [edge_attrs.start_marker, edge_attrs.mid_marker,
                edge_attrs.end_marker]:
        return shape_from_prop(val, edge_marker)

    if attr in [vertex_attrs.color, vertex_attrs.fill_color,
                vertex_attrs.text_color, edge_attrs.color]:
        if isinstance(val, list):
            return val
        if isinstance(val, (tuple, np.ndarray)):
            return list(val)
        if isinstance(val, str):
            return list(matplotlib.colors.ColorConverter().to_rgba(val))
        if isinstance(val, PropertyMap):
            if val.value_type() in ["vector<double>", "vector<long double>"]:
                return val
            if val.value_type() in ["int32_t", "int64_t", "double",
                                    "long double", "unsigned long", "bool"]:
                if val.fa is None:
                    vrange = val[val.get_graph().vertex(0)]
                    vrange = [vrange, vrange]
                    for v in val.get_graph().vertices():
                        vrange[0] = min(vrange[0], val[v])
                        vrange[1] = max(vrange[1], val[v])
                else:
                    vrange = [val.fa.min(), val.fa.max()]
                cnorm = matplotlib.colors.normalize(vmin=vrange[0],
                                                    vmax=vrange[1])
                if val.key_type() == "v":
                    prop = val.get_graph().new_vertex_property("vector<double>")
                    descs = val.get_graph().vertices()
                else:
                    prop = val.get_graph().new_edge_property("vector<double>")
                    descs = val.get_graph().edges()
                for v in descs:
                    prop[v] = cmap(cnorm(val[v]))
                return prop
            if val.value_type() == "string":
                if val.key_type() == "v":
                    prop = val.get_graph().new_vertex_property("vector<double>")
                    for v in val.get_graph().vertices():
                        prop[v] = matplotlib.colors.ColorConverter().to_rgba(val[v])
                elif val.key_type() == "e":
                    prop = val.get_graph().new_edge_property("vector<double>")
                    for e in val.get_graph().edges():
                        prop[e] = matplotlib.colors.ColorConverter().to_rgba(val[e])
                return prop
        raise ValueError("Invalid value for attribute %s: %s" %
                         (repr(attr), repr(val)))
    return val


def _attrs(attrs, d, g, cmap):
    nattrs = {}
    defaults = {}
    for k, v in attrs.items():
        try:
            if d == "v":
                attr = vertex_attrs.__dict__[k]
            else:
                attr = edge_attrs.__dict__[k]
        except KeyError:
            warnings.warn("Unknown attribute: " + k, UserWarning)
            continue
        if isinstance(v, PropertyMap):
            nattrs[int(attr)] = _prop(d, g, _convert(attr, v, cmap))
        else:
            defaults[int(attr)] = _convert(attr, v, cmap)
    return nattrs, defaults


def get_attr(attr, d, attrs, defaults):
    if attr in attrs:
        p = attrs[attr]
    else:
        p = defaults[attr]
    if isinstance(p, PropertyMap):
        return p[d]
    else:
        return p


def position_parallel_edges(g, pos, loop_angle=float("nan"),
                            parallel_distance=1):
    lp = label_parallel_edges(GraphView(g, directed=False))
    ll = label_self_loops(g)
    g = GraphView(g, directed=True)
    if ((len(lp.fa) == 0 or lp.fa.max() == 0) and
        (len(ll.fa) == 0 or ll.fa.max() == 0)):
        return []
    else:
        spline = g.new_edge_property("vector<double>")
        libgraph_tool_draw.put_parallel_splines(g._Graph__graph,
                                                _prop("v", g, pos),
                                                _prop("e", g, lp),
                                                _prop("e", g, spline),
                                                loop_angle,
                                                parallel_distance)
        return spline


def parse_props(prefix, args):
    props = {}
    others = {}
    for k, v in args.items():
        if k.startswith(prefix + "_"):
            props[k.replace(prefix + "_", "")] = v
        else:
            others[k] = v
    return props, others


def cairo_draw(g, pos, cr, vprops=None, eprops=None, vorder=None, eorder=None,
               nodesfirst=False, vcmap=matplotlib.cm.jet,
               ecmap=matplotlib.cm.jet, loop_angle=float("nan"),
               parallel_distance=None, **kwargs):
    r"""
    Draw a graph to a :mod:`cairo` context.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be drawn.
    pos : :class:`~graph_tool.PropertyMap`
        Vector-valued vertex property map containing the x and y coordinates of
        the vertices.
    cr : :class:`~cairo.Context`
        A :class:`~cairo.Context` instance.
    vprops : dict (optional, default: ``None``)
        Dictionary with the vertex properties. Individual properties may also be
        given via the ``vertex_<prop-name>`` parameters, where ``<prop-name>`` is
        the name of the property.
    eprops : dict (optional, default: ``None``)
        Dictionary with the vertex properties. Individual properties may also be
        given via the ``edge_<prop-name>`` parameters, where ``<prop-name>`` is
        the name of the property.
    vorder : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        If provided, defines the relative order in which the vertices are drawn.
    eorder : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        If provided, defines the relative order in which the edges are drawn.
    nodesfirst : bool (optional, default: ``False``)
        If ``True``, the vertices are drawn first, otherwise the edges are.
    vcmap : :class:`matplotlib.colors.Colormap` (default: :class:`matplotlib.cm.jet`)
        Vertex color map.
    ecmap : :class:`matplotlib.colors.Colormap` (default: :class:`matplotlib.cm.jet`)
        Edge color map.
    loop_angle : float (optional, default: ``nan``)
        Angle used to draw self-loops. If ``nan`` is given, they will be placed
        radially from the center of the layout.
    parallel_distance : float (optional, default: ``None``)
        Distance used between parallel edges. If not provided, it will be
        determined automatically.
    vertex_* : :class:`~graph_tool.PropertyMap` or arbitrary types (optional, default: ``None``)
        Parameters following the pattern ``vertex_<prop-name>`` specify the
        vertex property with name ``<prop-name>``, as an alternative to the
        ``vprops`` parameter.
    edge_* : :class:`~graph_tool.PropertyMap` or arbitrary types (optional, default: ``None``)
        Parameters following the pattern ``edge_<prop-name>`` specify the edge
        property with name ``<prop-name>``, as an alternative to the ``eprops``
        parameter.

    """

    vprops = {} if vprops is None else copy.copy(vprops)
    eprops = {} if eprops is None else copy.copy(eprops)

    props, kwargs = parse_props("vertex", kwargs)
    vprops.update(props)
    props, kwargs = parse_props("edge", kwargs)
    eprops.update(props)
    for k in kwargs:
        warnings.warn("Unknown parameter: " + k, UserWarning)

    if "control_points" not in eprops:
        if parallel_distance is None:
            parallel_distance = vprops.get("size", _vdefaults["size"])
            if isinstance(parallel_distance, PropertyMap):
                parallel_distance = parallel_distance.fa.mean()
            parallel_distance /= 1.5
            M = cr.get_matrix()
            scale = transform_scale(M, 1,)
            parallel_distance /= scale
        eprops["control_points"] = position_parallel_edges(g, pos, loop_angle,
                                                           parallel_distance)
    if g.is_directed() and "end_marker" not in eprops:
        eprops["end_marker"] = "arrow"
    vattrs, vdefaults = _attrs(vprops, "v", g, vcmap)
    eattrs, edefaults = _attrs(eprops, "e", g, ecmap)
    vdefs = _attrs(_vdefaults, "v", g, vcmap)[1]
    vdefs.update(vdefaults)
    edefs = _attrs(_edefaults, "e", g, ecmap)[1]
    edefs.update(edefaults)

    if "control_points" not in eprops:
        if parallel_distance is None:
            parallel_distance = _defaults
        eprops["control_points"] = position_parallel_edges(g, pos, loop_angle,
                                                           parallel_distance)

    g = GraphView(g, directed=True)
    libgraph_tool_draw.cairo_draw(g._Graph__graph, _prop("v", g, pos),
                                  _prop("v", g, vorder), _prop("e", g, eorder),
                                  nodesfirst, vattrs, eattrs, vdefs, edefs, cr)


def graph_draw(g, pos=None, vprops=None, eprops=None, vorder=None, eorder=None,
               nodesfirst=False, output_size=(600, 600), fit_view=True,
               output=None, fmt="auto", **kwargs):
    r"""Draw a graph to screen or to a file using :mod:`cairo`.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be drawn.
    pos : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Vector-valued vertex property map containing the x and y coordinates of
        the vertices. If not given, it will be computed using :func:`sfdp_layout`.
    vprops : dict (optional, default: ``None``)
        Dictionary with the vertex properties. Individual properties may also be
        given via the ``vertex_<prop-name>`` parameters, where ``<prop-name>`` is
        the name of the property.
    eprops : dict (optional, default: ``None``)
        Dictionary with the vertex properties. Individual properties may also be
        given via the ``edge_<prop-name>`` parameters, where ``<prop-name>`` is
        the name of the property.
    vorder : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        If provided, defines the relative order in which the vertices are drawn.
    eorder : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        If provided, defines the relative order in which the edges are drawn.
    nodesfirst : bool (optional, default: ``False``)
        If ``True``, the vertices are drawn first, otherwise the edges are.
    output_size : tuple of scalars (optional, default: ``(600,600)``)
        Size of the drawing canvas. The units will depend on the output format
        (pixels for the screen, points for PDF, etc).
    fit_view : bool (optional, default: ``True``)
        If ``True``, the layout will be scaled to fit the entire display area.
    output : string (optional, default: ``None``)
        Output file name. If not given, the graph will be displayed via
        :func:`interactive_window`.
    fmt : string (default: ``"auto"``)
        Output file format. Possible values are ``"auto"``, ``"ps"``, ``"pdf"``,
        ``"svg"``, and ``"png"``. If the value is ``"auto"``, the format is
        guessed from the ``output`` parameter.
    vertex_* : :class:`~graph_tool.PropertyMap` or arbitrary types (optional, default: ``None``)
        Parameters following the pattern ``vertex_<prop-name>`` specify the
        vertex property with name ``<prop-name>``, as an alternative to the
        ``vprops`` parameter.
    edge_* : :class:`~graph_tool.PropertyMap` or arbitrary types (optional, default: ``None``)
        Parameters following the pattern ``edge_<prop-name>`` specify the edge
        property with name ``<prop-name>``, as an alternative to the ``eprops``
        parameter.
    **kwargs
        Any extra parameters are passed to :func:`~graph_tool.draw.interactive_window`,
        :class:`~graph_tool.draw.GraphWindow`, :class:`~graph_tool.draw.GraphWidget`
        and :func:`~graph_tool.draw.cairo_draw`.

    Returns
    -------
    pos : :class:`~graph_tool.PropertyMap`
        Vector vertex property map with the x and y coordinates of the vertices.
    selected : :class:`~graph_tool.PropertyMap` (optional, only if ``output is None``)
        Boolean-valued vertex property map marking the vertices which were
        selected interactively.

    Notes
    -----


    .. table:: **List of vertex properties**

        +---------------+---------------------------------------------------+------------------------+----------------------------------+
        | Name          | Description                                       | Accepted types         | Default Value                    |
        +===============+===================================================+========================+==================================+
        | shape         | The vertex shape. Can be one of the following     | ``str`` or ``int``     | ``"circle"``                     |
        |               | strings: "circle", "triangle", "square",          |                        |                                  |
        |               | "pentagon", "hexagon", "heptagon", "octagon"      |                        |                                  |
        |               | "double_circle", "double_triangle",               |                        |                                  |
        |               | "double_square", "double_pentagon",               |                        |                                  |
        |               | "double_hexagon", "double_heptagon",              |                        |                                  |
        |               | "double_octagon".                                 |                        |                                  |
        |               | Optionally, this might take a numeric value       |                        |                                  |
        |               | corresponding to position in the list above.      |                        |                                  |
        +---------------+---------------------------------------------------+------------------------+----------------------------------+
        | color         | Color used to stroke the lines of the vertex.     | ``str`` or list of     | ``[0., 0., 0., 1]``              |
        |               |                                                   | ``floats``             |                                  |
        +---------------+---------------------------------------------------+------------------------+----------------------------------+
        | fill_color    | Color used to fill the interior of the vertex.    | ``str`` or list of     | ``[0.640625, 0, 0, 0.9]``        |
        |               |                                                   | ``floats``             |                                  |
        +---------------+---------------------------------------------------+------------------------+----------------------------------+
        | size          | The size of the vertex, in the default units of   | ``float`` or ``int``   | ``5``                            |
        |               | the output format (normally either pixels or      |                        |                                  |
        |               | points).                                          |                        |                                  |
        +---------------+---------------------------------------------------+------------------------+----------------------------------+
        | pen_width     | Width of the lines used to draw the vertex, in    | ``float`` or ``int``   | ``0.8``                          |
        |               | the default units of the output format (normally  |                        |                                  |
        |               | either pixels or points).                         |                        |                                  |
        +---------------+---------------------------------------------------+------------------------+----------------------------------+
        | halo          | Whether to draw a circular halo around the        | ``bool``               | ``False``                        |
        |               | vertex.                                           |                        |                                  |
        +---------------+---------------------------------------------------+------------------------+----------------------------------+
        | halo_color    | Color used to draw the halo.                      | ``str`` or list of     | ``[0., 0., 1., 0.5]``            |
        |               |                                                   | ``floats``             |                                  |
        +---------------+---------------------------------------------------+------------------------+----------------------------------+
        | text          | Text to draw together with the vertex.            | ``str``                | ``""``                           |
        +---------------+---------------------------------------------------+------------------------+----------------------------------+
        | text_color    | Color used to draw the text.                      | ``str`` or list of     | ``[0., 0., 0., 1.]``             |
        |               |                                                   | ``floats``             |                                  |
        +---------------+---------------------------------------------------+------------------------+----------------------------------+
        | text_position | Position of the text relative to the vertex.      | ``float`` or ``int``   | ``-1``                           |
        |               | If the passed value is positive, it will          |                        |                                  |
        |               | correspond to an angle in radians, which will     |                        |                                  |
        |               | determine where the text will be placed outside   |                        |                                  |
        |               | the vertex. If the value is negative, the text    |                        |                                  |
        |               | will be placed inside the vertex.                 |                        |                                  |
        +---------------+---------------------------------------------------+------------------------+----------------------------------+
        | font_family   | Font family used to draw the text.                | ``str``                | ``"serif"``                      |
        +---------------+---------------------------------------------------+------------------------+----------------------------------+
        | font_slant    | Font slant used to draw the text.                 | ``cairo.FONT_SLANT_*`` | :data:`cairo.FONT_SLANT_NORMAL`  |
        +---------------+---------------------------------------------------+------------------------+----------------------------------+
        | font_weight   | Font weight used to draw the text.                | ``cairo.FONT_WEIGHT_*``| :data:`cairo.FONT_WEIGHT_NORMAL` |
        +---------------+---------------------------------------------------+------------------------+----------------------------------+
        | font_size     | Font size used to draw the text.                  | ``float`` or ``int``   | ``12``                           |
        +---------------+---------------------------------------------------+------------------------+----------------------------------+


    .. table:: **List of edge properties**

        +----------------+---------------------------------------------------+------------------------+----------------------------------+
        | Name           | Description                                       | Accepted types         | Default Value                    |
        +================+===================================================+========================+==================================+
        | color          | Color used to stroke the edge lines.              | ``str`` or list of     | ``[0.179, 0.203,0.210, 0.8]``    |
        |                |                                                   | floats                 |                                  |
        +----------------+---------------------------------------------------+------------------------+----------------------------------+
        | pen_width      | Width of the line used to draw the edge, in       | ``float`` or ``int``   | ``1.0``                          |
        |                | the default units of the output format (normally  |                        |                                  |
        |                | either pixels or points).                         |                        |                                  |
        +----------------+---------------------------------------------------+------------------------+----------------------------------+
        | start_marker,  | Edge markers. Can be one of "none", "arrow",      | ``str`` or ``int``     | ``-1``                           |
        | mid_marker,    | "circle", "square", "diamond", or "bar".          |                        |                                  |
        | end_marker     | Optionally, this might take a numeric value       |                        |                                  |
        |                | corresponding to position in the list above.      |                        |                                  |
        +----------------+---------------------------------------------------+------------------------+----------------------------------+
        | marker_size    | Size of edge markers, in units appropriate to the | ``float`` or ``int``   | ``4``                            |
        |                | output format (normally either pixels or points). |                        |                                  |
        +----------------+---------------------------------------------------+------------------------+----------------------------------+
        | control_points | Control points of a Bézier spline used to draw    | sequence of ``floats`` | ``[]``                           |
        |                | the edge.                                         |                        |                                  |
        +----------------+---------------------------------------------------+------------------------+----------------------------------+

    Examples
    --------
    >>> from numpy import *
    >>> from numpy.random import seed, zipf
    >>> seed(43)
    >>> g = gt.price_network(1500)
    >>> deg = g.degree_property_map("in")
    >>> deg.a = 4 * (sqrt(deg.a) * 0.5 + 0.4)
    >>> ebet = gt.betweenness(g)[1]
    >>> ebet.a /= ebet.a.max() / 10.
    >>> eorder = ebet.copy()
    >>> eorder.a *= -1
    >>> pos = gt.sfdp_layout(g)
    >>> control = g.new_edge_property("vector<double>")
    >>> for e in g.edges():
    ...     d = sqrt(sum((pos[e.source()].a - pos[e.target()].a) ** 2)) / 5
    ...     control[e] = [0.3, d, 0.7, d]
    >>> gt.graph_draw(g, pos=pos, vertex_size=deg, vertex_fill_color=deg, vorder=deg,
    ...               edge_color=ebet, eorder=eorder, edge_pen_width=ebet,
    ...               edge_control_points=control, # some curvy edges
    ...               output="graph-draw.pdf")
    <...>

    .. figure:: graph-draw.*
        :align: center

        SFDP force-directed layout of a Price network with 1500 nodes. The
        vertex size and color indicate the degree, and the edge color and width
        the edge betweeness centrality.

    """

    vprops = vprops.copy() if vprops is not None else {}
    eprops = eprops.copy() if eprops is not None else {}

    props, kwargs = parse_props("vertex", kwargs)
    vprops.update(props)
    props, kwargs = parse_props("edge", kwargs)
    eprops.update(props)

    if pos is None:
        if (g.num_vertices() > 2 and output is None and
            kwargs.get("update_layout", True)):
            L = np.sqrt(g.num_vertices())
            pos = random_layout(g, [L, L])
            if g.num_vertices() > 1000:
                if "multilevel" not in kwargs:
                    kwargs["multilevel"] = True
            if "layout_K" not in kwargs:
                kwargs["layout_K"] = _avg_edge_distance(g, pos) / 10
        else:
            pos = sfdp_layout(g)
    elif output is None:
        if "layout_K" not in kwargs:
            kwargs["layout_K"] = _avg_edge_distance(g, pos)
        if "update_layout" not in kwargs:
            kwargs["update_layout"] = False

    if output is None:
        return interactive_window(g, pos, vprops, eprops, vorder, eorder,
                                  nodesfirst, **kwargs)
    else:
        output = os.path.expanduser(output)
        base, ext = os.path.splitext(output)
        if ext == ".gz":
            out = gzip.GzipFile(output, "w")
            output = base
        elif ext == ".bz2":
            out = bz2.BZ2File(output, "w")
            output = base
        elif ext == ".zip":
            out = zipfile.ZipFile(output, "w")
            output = base
        else:
            out = open(output, "w")

        if fmt == "auto":
            fmt = os.path.splitext(output)[1].replace(".", "")
        if fmt == "pdf":
            srf = cairo.PDFSurface(out, output_size[0], output_size[1])
        elif fmt == "ps":
            srf = cairo.PSSurface(out, output_size[0], output_size[1])
        elif fmt == "svg":
            srf = cairo.SVGSurface(out, output_size[0], output_size[1])
        elif fmt == "png":
            srf = cairo.ImageSurface(cairo.FORMAT_ARGB32, output_size[0],
                                     output_size[1])
        else:
            raise ValueError("Invalid format type: " + fmt)

        cr = cairo.Context(srf)

        adjust_default_sizes(g, output_size, vprops, eprops)
        if fit_view:
            offset, zoom = fit_to_view(g, pos, output_size, vprops["size"],
                                       vprops["pen_width"], None,
                                       vprops.get("text", None),
                                       vprops.get("font_family",
                                                  _vdefaults["font_family"]),
                                       vprops.get("font_size",
                                                  _vdefaults["font_size"]),
                                       cr)
        else:
            offset, zoom = [0, 0], 1

        if "bg_color" in kwargs:
            bg_color = kwargs["bg_color"]
            del  kwargs["bg_color"]
            cr.set_source_rgba(bg_color[0], bg_color[1],
                               bg_color[2], bg_color[3])
            cr.paint()
        cr.translate(offset[0], offset[1])
        cr.scale(zoom, zoom)

        cairo_draw(g, pos, cr, vprops, eprops, vorder, eorder,
                   nodesfirst, **kwargs)
        del cr

        if output.endswith(".png"):
            srf.write_to_png(out)
        return pos


def adjust_default_sizes(g, geometry, vprops, eprops, force=False):
    if "size" not in vprops or force:
        A = geometry[0] * geometry[1]
        vprops["size"] = np.sqrt(A / g.num_vertices()) / 3.5

    if "pen_width" not in vprops or force:
        size = vprops["size"]
        if isinstance(vprops["size"], PropertyMap):
            size = vprops["size"].fa.mean()
        vprops["pen_width"] = size / 10
        if "pen_width" not in eprops or force:
            eprops["pen_width"] = size / 10
        if "marker_size" not in eprops or force:
            eprops["marker_size"] = size * 0.8


def scale_ink(scale, vprops, eprops):
    if "size" not in vprops:
        vprops["size"] = _vdefaults["size"]
    if "pen_width" not in vprops:
        vprops["pen_width"] = _vdefaults["pen_width"]
    if "font_size" not in vprops:
        vprops["font_size"] = _vdefaults["font_size"]
    if "pen_width" not in eprops:
        eprops["pen_width"] = _edefaults["pen_width"]
    if "marker_size" not in eprops:
        eprops["marker_size"] = _edefaults["marker_size"]

    for props in [vprops, eprops]:
        if isinstance(props["pen_width"], PropertyMap):
            props["pen_width"].fa *= scale
        else:
            props["pen_width"] *= scale
    if isinstance(vprops["size"], PropertyMap):
        vprops["size"].fa *= scale
    else:
        vprops["size"] *= scale
    if isinstance(vprops["font_size"], PropertyMap):
        vprops["font_size"].fa *= scale
    else:
        vprops["font_size"] *= scale
    if isinstance(eprops["marker_size"], PropertyMap):
        eprops["marker_size"].fa *= scale
    else:
        eprops["marker_size"] *= scale

def get_bb(g, pos, size, pen_width, size_scale=1, text=None, font_family=None,
           font_size=None, cr=None):
    size = size.fa[:g.num_vertices()] if isinstance(size, PropertyMap) else size
    pen_width = pen_width.fa if isinstance(pen_width, PropertyMap) else pen_width
    pos_x, pos_y = ungroup_vector_property(pos, [0, 1])
    if text is not None and text != "":
        if not isinstance(size, PropertyMap):
            uniform = (not isinstance(font_size, PropertyMap) and
                       not isinstance(font_family, PropertyMap))
            size = np.ones(len(pos_x.fa)) * size
        else:
            uniform = False
        for i, v in enumerate(g.vertices()):
            ff = font_family[v] if isinstance(font_family, PropertyMap) \
               else font_family
            cr.select_font_face(ff)
            fs = font_size[v] if isinstance(font_family, PropertyMap) \
               else font_size
            cr.set_font_size(fs)
            t = text[v] if isinstance(text, PropertyMap) else text
            if not isinstance(t, str):
                t = str(t)
            extents = cr.text_extents(t)
            s = max(extents[2], extents[3]) * 1.4
            size[i] = max(size[i] * size_scale, s) / size_scale
            if uniform:
                size[:] = size[i]
                break
    delta = (size * size_scale) / 2 + pen_width
    x_range = [pos_x.fa.min(), pos_x.fa.max()]
    y_range = [pos_y.fa.min(), pos_y.fa.max()]
    x_delta = [x_range[0] - (pos_x.fa - delta).min(),
               (pos_x.fa + delta).max() - x_range[1]]
    y_delta = [y_range[0] - (pos_y.fa - delta).min(),
               (pos_y.fa + delta).max() - y_range[1]]
    return x_range, y_range, x_delta, y_delta


def fit_to_view(g, pos, geometry, size, pen_width, M=None, text=None,
                font_family=None, font_size=None, cr=None):
    if M is not None:
        pos_x, pos_y = ungroup_vector_property(pos, [0, 1])
        P = np.zeros((2, len(pos_x.fa)))
        P[0, :] = pos_x.fa
        P[1, :] = pos_y.fa
        T = np.zeros((2, 2))
        O = np.zeros(2)
        T[0, 0], T[1, 0], T[0, 1], T[1, 1], O[0], O[1] = M
        P = np.dot(T, P)
        P[0] += O[0]
        P[1] += O[1]
        pos_x.fa = P[0, :]
        pos_y.fa = P[1, :]
        pos = group_vector_property([pos_x, pos_y])
    x_range, y_range, x_delta, y_delta = get_bb(g, pos, size, pen_width,
                                                1, text, font_family,
                                                font_size, cr)
    zoom_x = (geometry[0] - sum(x_delta)) / (x_range[1] - x_range[0])
    zoom_y = (geometry[1] - sum(y_delta)) / (y_range[1] - y_range[0])
    if np.isnan(zoom_x) or np.isinf(zoom_x) or zoom_x == 0:
        zoom_x = 1
    if np.isnan(zoom_y) or np.isinf(zoom_y) or zoom_y == 0:
        zoom_y = 1
    pad = 0.95
    zoom = min(zoom_x, zoom_y) * pad
    empty_x = (geometry[0] - sum(x_delta)) - (x_range[1] - x_range[0]) * zoom
    empty_y = (geometry[1] - sum(y_delta)) - (y_range[1] - y_range[0]) * zoom
    offset = [-x_range[0] * zoom + empty_x / 2 + x_delta[0],
              -y_range[0] * zoom + empty_y / 2 + y_delta[0]]
    return offset, zoom


def transform_scale(M, scale):
    p = M.transform_distance(scale / np.sqrt(2),
                             scale / np.sqrt(2))
    return np.sqrt(p[0] ** 2 + p[1] ** 2)


#
# The functions and classes below depend on GTK
# =============================================
#

try:
    from gi.repository import Gtk, Gdk, GdkPixbuf
    import gobject
    from gtk_draw import *
except (ImportError, RuntimeError), e:
    msg = "Error importing Gtk module: %s; GTK+ drawing will not work." % str(e)
    warnings.filterwarnings("always", msg, ImportWarning)
    warnings.warn(msg, ImportWarning)
