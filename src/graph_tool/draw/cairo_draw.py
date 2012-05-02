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
    warnings.warn("Error importing cairo. Graph drawing will not work.",
                  ImportWarning)

try:
    import matplotlib.cm
    import matplotlib.colors
except ImportError:
    warnings.warn("error importing matplotlib module. " + \
                  "Graph drawing will not work..", ImportWarning)

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
    warnings.warn("error importing cairo-based drawing library. " +
                  "Was graph-tool compiled with cairomm support?",
                  ImportWarning)

try:
    from gi.repository import Gtk, Gdk, GdkPixbuf
    import gobject
except ImportError:
    warnings.warn("Error importing Gtk module. Gtk drawing will " +
                  "not work.", ImportWarning)

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
    if lp.fa.max() == 0 and ll.fa.max() == 0:
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


def transform_scale(M, scale):
    p = M.transform_distance(scale / np.sqrt(2),
                             scale / np.sqrt(2))
    return np.sqrt(p[0] ** 2 + p[1] ** 2)


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


def point_in_poly(p, poly):
    i, c = 0, False
    j = len(poly) - 1
    while i < len(poly):
        if (((poly[i][1] > p[1]) != (poly[j][1] > p[1])) and
            (p[0] < (poly[j][0] - poly[i][0]) * (p[1] - poly[i][1]) /
             (poly[j][1] - poly[i][1]) + poly[i][0])):
            c = not c
        j = i
        i += 1
    return c


class VertexMatrix(object):
    def __init__(self, g, pos):
        self.g = g
        self.pos = pos
        self.m = None
        self.m_res = None
        self.update()

    def get_box(self, p, size=None):
        if size is None:
            return (int(round(p[0] / self.m_res)),
                    int(round(p[1] / self.m_res)))
        else:
            n = int(np.ceil(size / self.m_res))
            b = self.get_box(p)
            boxes = []
            for i in xrange(-n, n):
                for j in xrange(-n, n):
                    boxes.append((b[0] + i, b[1] + j))
            return boxes

    def update(self):
        pos_x, pos_y = ungroup_vector_property(self.pos, [0, 1])
        x_range = [pos_x.fa.min(), pos_x.fa.max()]
        y_range = [pos_y.fa.min(), pos_y.fa.max()]
        self.m_res = min(x_range[1] - x_range[0],
                         y_range[1] - y_range[0]) / np.sqrt(self.g.num_vertices())
        self.m_res *= np.sqrt(10)

        self.m = defaultdict(set)
        for v in self.g.vertices():
            i, j = self.get_box(self.pos[v])
            self.m[(i, j)].add(v)

    def update_vertex(self, v, new_pos):
        b = self.get_box(self.pos[v])
        self.m[b].remove(v)
        self.pos[v] = new_pos
        b = self.get_box(self.pos[v])
        self.m[b].add(v)

    def remove_vertex(self, v):
        b = self.get_box(self.pos[v])
        self.m[b].remove(v)

    def add_vertex(self, v):
        b = self.get_box(self.pos[v])
        self.m[b].add(v)

    def get_closest(self, pos):
        pos = np.array(pos)
        box = self.get_box(pos)
        dist = float("inf")
        clst = None
        for i in xrange(-1, 2):
            for j in xrange(-1, 2):
                b = (box[0] + i, box[1] + j)
                for v in self.m[b]:
                    ndist = ((pos - self.pos[v].a[:2]) ** 2).sum()
                    if ndist < dist:
                        dist = ndist
                        clst = v
        return clst

    def mark_polygon(self, points, selected):
        rect = [min([x[0] for x in points]), min([x[1] for x in points]),
                max([x[0] for x in points]), max([x[1] for x in points])]
        p1 = self.get_box(rect[:2])
        p2 = self.get_box(rect[2:])
        for i in xrange(p1[0], p2[0] + 1):
            for j in xrange(p1[1], p2[1] + 1):
                for v in self.m[(i, j)]:
                    p = self.pos[v]
                    if not point_in_poly(p, points):
                        continue
                    selected[v] = True


def apply_transforms(g, pos, m):
    m = tuple(m)
    g = GraphView(g, directed=True)
    libgraph_tool_draw.apply_transforms(g._Graph__graph, _prop("v", g, pos),
                                        m[0], m[1], m[2], m[3], m[4], m[5])


class GraphWidget(Gtk.DrawingArea):
    r"""Interactive GTK+ widget displaying a given graph.

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
    update_layout : bool (optional, default: ``True``)
        If ``True``, the layout will be updated dynamically.
    layout_K : float (optional, default: ``1.0``)
        Parameter ``K`` passed to :func:`~graph_tool.draw.sfdp_layout`.
    multilevel : bool (optional, default: ``False``)
        Parameter ``multilevel`` passed to :func:`~graph_tool.draw.sfdp_layout`.
    display_props : list of :class:`~graph_tool.PropertyMap` instances (optional, default: ``None``)
        List of properties to be displayed when the mouse passes over a vertex.
    display_props_size : float (optional, default: ``11``)
        Font size used to display the vertex properties.
    bg_color : str or sequence (optional, default: ``None``)
        Background color. The default is white.
    vertex_* : :class:`~graph_tool.PropertyMap` or arbitrary types (optional, default: ``None``)
        Parameters following the pattern ``vertex_<prop-name>`` specify the
        vertex property with name ``<prop-name>``, as an alternative to the
        ``vprops`` parameter.
    edge_* : :class:`~graph_tool.PropertyMap` or arbitrary types (optional, default: ``None``)
        Parameters following the pattern ``edge_<prop-name>`` specify the edge
        property with name ``<prop-name>``, as an alternative to the ``eprops``
        parameter.
    **kwargs
        Any extra parameters are passed to :func:`~graph_tool.draw.cairo_draw`.

    Notes
    -----

    The graph drawing can be panned by dragging with the middle mouse button
    pressed. The graph may be zoomed by scrolling with the mouse wheel, or
    equivalent (if the "shift" key is held, the vertex/edge sizes are scaled
    accordingly). The layout may be rotated by dragging while holding the
    "control" key. Pressing the "r" key centers and zooms the layout around the
    graph.  By pressing the "a" key, the current translation, scaling and
    rotation transformations are applied to the vertex positions themselves, and
    the transformation matrix is reset (if this is never done, the given
    position properties are never modified).

    Individual vertices may be selected by pressing the left mouse button. The
    currently selected vertex follows the mouse pointer. To stop the selection,
    the right mouse button must be pressed. Alternatively, a group of vertices
    may be selected by holding the "shift" button while the pointer is dragged
    while pressing the left button. The selected vertices may be moved by
    dragging the pointer with the left button pressed. They may be rotated by
    holding the "control" key and scrolling with the mouse. If the key "z" is
    pressed, the layout is zoomed to fit the selected vertices only.

    If the key "s" is pressed, the dynamic spring-block layout is
    activated. Vertices which are currently selected do are not updated.

    """

    def __init__(self, g, pos, vprops=None, eprops=None, vorder=None,
                 eorder=None, nodesfirst=False, update_layout=False,
                 layout_K=1., multilevel=False, display_props=None,
                 display_props_size=11, bg_color=None, **kwargs):
        Gtk.DrawingArea.__init__(self)

        vprops = {} if vprops is None else vprops
        eprops = {} if eprops is None else eprops

        props, kwargs = parse_props("vertex", kwargs)
        vprops.update(props)
        props, kwargs = parse_props("edge", kwargs)
        eprops.update(props)
        self.kwargs = kwargs

        self.g = g
        self.pos = pos
        self.vprops = vprops
        self.eprops = eprops
        self.vorder = vorder
        self.eorder = eorder
        self.nodesfirst = nodesfirst

        self.panning = None
        self.tmatrix = cairo.Matrix()  # position to surface
        self.smatrix = cairo.Matrix()  # surface to screen
        self.pointer = [0, 0]
        self.picked = False
        self.selected = g.new_vertex_property("bool")
        self.srect = None
        self.drag_begin = None
        self.moved_picked = False
        self.vertex_matrix = None

        self.display_prop = g.vertex_index if display_props is None \
                            else display_props
        self.display_prop_size = display_props_size

        self.geometry = None
        self.base = None
        self.background = None
        self.bg_color = bg_color if bg_color is not None else [1, 1, 1, 1]
        self.surface_callback = None

        self.layout_callback_id = None
        self.layout_K = layout_K
        self.layout_init_step = self.layout_K
        self.epsilon = 0.01 * self.layout_K
        self.multilevel_layout = multilevel

        if multilevel:
            self.cgs = coarse_graphs(g)
            u = self.cgs.next()
            self.cg, self.cpos, self.layout_K, self.cvcount, self.cecount = u
            self.ag = self.g
            self.apos = self.pos
            self.g = self.cg
            self.pos = self.cpos
            self.layout_step = self.layout_K
        else:
            self.cg = None
        if update_layout:
            self.reset_layout()

        # Event signals
        self.connect("motion_notify_event", self.motion_notify_event)
        self.connect("button_press_event", self.button_press_event)
        self.connect("button_release_event", self.button_release_event)
        self.connect("scroll_event", self.scroll_event)
        self.connect("key_press_event", self.key_press_event)
        self.connect("key_release_event", self.key_release_event)
        self.connect("destroy_event", self.cleanup)

        self.set_events(Gdk.EventMask.EXPOSURE_MASK
                        | Gdk.EventMask.LEAVE_NOTIFY_MASK
                        | Gdk.EventMask.BUTTON_PRESS_MASK
                        | Gdk.EventMask.BUTTON_RELEASE_MASK
                        | Gdk.EventMask.BUTTON_MOTION_MASK
                        | Gdk.EventMask.POINTER_MOTION_MASK
                        | Gdk.EventMask.POINTER_MOTION_HINT_MASK
                        | Gdk.EventMask.SCROLL_MASK
                        | Gdk.EventMask.KEY_PRESS_MASK
                        | Gdk.EventMask.KEY_RELEASE_MASK)

        self.set_property("can-focus", True)
        self.connect("draw", self.draw)

    def cleanup(self):
        """Cleanup callbacks."""
        if self.layout_callback_id is not None:
            ret = gobject.source_remove(self.layout_callback_id)
            if not ret:
                warnings.warn("error removing idle callback...")
            self.layout_callback_id = None

    def __del__(self):
        self.cleanup()

    # Layout update

    def reset_layout(self):
        """Reset the layout algorithm."""
        if self.layout_callback_id is not None:
            gobject.source_remove(self.layout_callback_id)
            self.layout_callback_id = None
        self.layout_step = self.layout_init_step
        self.layout_callback_id = gobject.idle_add(self.layout_callback)

    def layout_callback(self):
        """Perform one step of the layout algorithm."""
        if self.layout_callback_id is None:
            return False
        pos_temp = ungroup_vector_property(self.pos, [0, 1])
        sfdp_layout(self.g, K=self.layout_K,
                    max_iter=5, pos=self.pos,
                    pin=self.selected,
                    init_step=self.layout_step,
                    multilevel=False)
        self.layout_step *= 0.9
        if self.vertex_matrix is not None:
            self.vertex_matrix.update()
        self.regenerate_surface(lazy=False)
        self.queue_draw()
        ps = ungroup_vector_property(self.pos, [0, 1])
        delta = np.sqrt((pos_temp[0].fa - ps[0].fa) ** 2 +
                        (pos_temp[1].fa - ps[1].fa) ** 2).mean()
        if delta > self.epsilon:
            return True
        else:
            if self.multilevel_layout:
                try:
                    u = self.cgs.next()
                    self.cg, self.cpos, K, self.cvcount, self.cecount = u
                    self.layout_K *= 0.75
                    self.g = self.cg
                    self.pos = self.cpos
                    self.layout_step = max(self.layout_K,
                                           _avg_edge_distance(self.g,
                                                              self.pos) / 10)
                    if self.vertex_matrix is not None:
                        self.vertex_matrix = VertexMatrix(self.g, self.pos)
                    self.epsilon = 0.05 * self.layout_K * self.g.num_edges()
                    geometry = [self.get_allocated_width(),
                                self.get_allocated_height()]
                    adjust_default_sizes(self.g, geometry, self.vprops,
                                         self.eprops, force=True)
                    self.fit_to_window(ink=False)
                    self.regenerate_surface(lazy=False)
                except StopIteration:
                    self.g = self.ag
                    self.pos = self.apos
                    self.g.copy_property(self.cpos, self.pos)
                    if self.vertex_matrix is not None:
                        self.vertex_matrix = VertexMatrix(self.g, self.pos)
                    self.multilevel_layout = False
                    self.layout_init_step = max(self.layout_K,
                                                _avg_edge_distance(self.g,
                                                                   self.pos) /
                                                                   10)
                    self.epsilon = 0.01 * self.layout_K

                return True
            self.layout_callback_id = None
            return False

    # Actual drawing

    def regenerate_surface(self, lazy=True, timeout=350):
        r"""Redraw the graph surface. If lazy is True, the actual redrawing will
        be performed after the specified timeout."""
        if lazy:
            if self.surface_callback is not None:
                gobject.source_remove(self.surface_callback)
            f = lambda: self.regenerate_surface(lazy=False)
            self.surface_callback = gobject.timeout_add(timeout, f)
        else:
            geometry = [self.get_allocated_width() * 3,
                        self.get_allocated_height() * 3]

            m = cairo.Matrix()
            m.translate(self.get_allocated_width(),
                        self.get_allocated_height())
            self.smatrix = self.smatrix * m
            self.tmatrix = self.tmatrix * self.smatrix
            if (self.base is None or self.base.get_width() != geometry[0] or
                self.base.get_height() != geometry[1]):
                # self.base = cairo.ImageSurface(cairo.FORMAT_ARGB32,
                #                                *geometry)
                w = self.get_window()
                if w is None:
                    return False
                self.base = w.create_similar_surface(cairo.CONTENT_COLOR_ALPHA,
                                                     *geometry)
            cr = cairo.Context(self.base)
            cr.set_source_rgba(*self.bg_color)
            cr.paint()
            cr.set_matrix(self.tmatrix)
            cairo_draw(self.g, self.pos, cr, self.vprops, self.eprops,
                       self.vorder, self.eorder, self.nodesfirst, **self.kwargs)
            self.smatrix = cairo.Matrix()
            self.smatrix.translate(-self.get_allocated_width(),
                                   -self.get_allocated_height())
            if self.surface_callback is not None:
                gobject.source_remove(self.surface_callback)
                self.surface_callback = None
                self.queue_draw()
            return False

    def draw(self, da, cr):
        r"""Redraw the widget."""

        geometry = [self.get_allocated_width(),
                    self.get_allocated_height()]

        if self.geometry is None:
            adjust_default_sizes(self.g, geometry, self.vprops, self.eprops)
            self.fit_to_window(ink=False)
            self.regenerate_surface(lazy=False)
            self.geometry = geometry

        ul = self.pos_to_device((0, 0), surface=True)
        lr = self.pos_to_device((self.base.get_width(),
                                 self.base.get_height()),
                                surface=True)
        if (ul[0] > 0 or lr[0] < geometry[0] or
            ul[1] > 0 or lr[1] < geometry[1]):
            self.regenerate_surface()

        if self.background is None:
            # draw checkerboard
            self.background = cairo.ImageSurface(cairo.FORMAT_ARGB32, 14, 14)
            bcr = cairo.Context(self.background)
            bcr.rectangle(0, 0, 7, 7)
            bcr.set_source_rgb(102. / 256, 102. / 256, 102. / 256)
            bcr.fill()
            bcr.rectangle(7, 0, 7, 7)
            bcr.set_source_rgb(153. / 256, 153. / 256, 153. / 256)
            bcr.fill()
            bcr.rectangle(0, 7, 7, 7)
            bcr.set_source_rgb(153. / 256, 153. / 256, 153. / 256)
            bcr.fill()
            bcr.rectangle(7, 7, 7, 7)
            bcr.set_source_rgb(102. / 256, 102. / 256, 102. / 256)
            bcr.fill()
            del bcr
            self.background = cairo.SurfacePattern(self.background)
            self.background.set_extend(cairo.EXTEND_REPEAT)

        cr.set_source(self.background)
        cr.paint()

        cr.save()
        cr.set_matrix(self.smatrix)
        cr.set_source_surface(self.base)
        cr.paint()
        cr.restore()

        if self.picked is not None or self.picked is not False:
            vprops = {}
            vprops.update(self.vprops)
            vprops["halo"] = True
            vprops["color"] = [1., 1., 1., 0.]
            vprops["fill_color"] = [1., 1., 1., 0.]
            vprops["text_color"] = [1., 1., 1., 0.]

            eprops = {}
            eprops.update(self.eprops)
            eprops["color"] = [1., 1., 1., 0.]

            u = GraphView(self.g, vfilt=self.selected)

            cr.save()
            cr.set_matrix(self.tmatrix * self.smatrix)
            cairo_draw(u, self.pos, cr, vprops, eprops, self.vorder,
                       self.eorder, self.nodesfirst)
            cr.restore()

        if self.srect is not None:
            cr.move_to(self.srect[0], self.srect[1])
            cr.line_to(self.srect[0], self.srect[3])
            cr.line_to(self.srect[2], self.srect[3])
            cr.line_to(self.srect[2], self.srect[1])
            cr.line_to(self.srect[0], self.srect[1])
            cr.close_path()
            cr.set_source_rgba(0, 0, 1, 0.3)
            cr.fill()

        if self.surface_callback is not None:
            icon = self.render_icon(Gtk.STOCK_EXECUTE, Gtk.IconSize.BUTTON)
            Gdk.cairo_set_source_pixbuf(cr, icon, 10, 10)
            cr.paint()

        if (self.picked is not None and self.picked is not False and
            not isinstance(self.picked, PropertyMap)):
            if isinstance(self.display_prop, PropertyMap):
                txt = str(self.display_prop[self.picked])
            else:
                txt = ", ".join([str(x[self.picked])
                                 for x in self.display_prop])
            geometry = [self.get_allocated_width(),
                        self.get_allocated_height()]
            pos = [10, geometry[1] - 10]
            cr.set_font_size(self.display_prop_size)
            ext = cr.text_extents(txt)
            pad = 8
            cr.rectangle(pos[0] - pad / 2, pos[1] - ext[3] - pad / 2,
                         ext[2] + pad, ext[3] + pad)
            cr.set_source_rgba(1, 1, 1, 1.0)
            cr.fill()
            cr.move_to(pos[0], pos[1])
            cr.set_source_rgba(0, 0, 0, 1.0)
            cr.show_text(txt)

        return False

    # Position and transforms

    def pos_to_device(self, pos, dist=False, surface=False, cr=None):
        """Convert a position from the graph space to the widget space."""
        if cr is None:
            cr = Gdk.cairo_create(self.get_root_window())
            if surface:
                cr.set_matrix(self.smatrix)
            else:
                cr.set_matrix(self.tmatrix * self.smatrix)
        if dist:
            return cr.user_to_device_distance(pos[0], pos[1])
        else:
            return cr.user_to_device(pos[0], pos[1])

    def pos_from_device(self, pos, dist=False, surface=False, cr=None):
        """Convert a position from the widget space to the device space."""
        if cr is None:
            cr = Gdk.cairo_create(self.get_root_window())
            if surface:
                cr.set_matrix(self.smatrix)
            else:
                cr.set_matrix(self.tmatrix * self.smatrix)
        if dist:
            return cr.device_to_user_distance(pos[0], pos[1])
        else:
            return cr.device_to_user(pos[0], pos[1])

    def apply_transform(self):
        r"""Apply current transform matrix to vertex coordinates."""
        zoom = self.pos_from_device((1, 0), dist=True)[0]
        apply_transforms(self.g, self.pos, self.smatrix * self.tmatrix)
        self.tmatrix = cairo.Matrix()
        self.tmatrix.scale(zoom, zoom)
        self.smatrix = cairo.Matrix()
        apply_transforms(self.g, self.pos, self.smatrix * self.tmatrix)
        self.tmatrix = cairo.Matrix()
        self.tmatrix.scale(1. / zoom, 1. / zoom)
        if self.vertex_matrix is not None:
            self.vertex_matrix.update()
        self.fit_to_window()
        self.regenerate_surface()
        self.queue_draw()

    def fit_to_window(self, ink=False, g=None):
        r"""Fit graph to window."""
        geometry = [self.get_allocated_width(), self.get_allocated_height()]
        if g is None:
            g = self.g
        pos = g.own_property(self.pos)
        cr = self.get_window().cairo_create()
        offset, zoom = fit_to_view(g, pos, geometry,
                                   self.vprops.get("size", 0),
                                   self.vprops.get("pen_width", 0),
                                   self.tmatrix * self.smatrix,
                                   self.vprops.get("text", None),
                                   self.vprops.get("font_family",
                                                   _vdefaults["font_family"]),
                                   self.vprops.get("font_size",
                                                   _vdefaults["font_size"]),
                                   cr)
        m = cairo.Matrix()
        m.translate(offset[0], offset[1])
        m.scale(zoom, zoom)
        self.tmatrix = self.tmatrix * self.smatrix * m
        self.smatrix = cairo.Matrix()
        if ink:
            scale_ink(zoom, self.vprops, self.eprops)

    # Picking vertices

    def init_picked(self):
        r"""Init picked vertices."""
        self.selected.fa = False
        p = self.pos_from_device(self.pointer)
        if self.vertex_matrix is None:
            self.vertex_matrix = VertexMatrix(self.g, self.pos)
        self.picked = self.vertex_matrix.get_closest(p)
        if self.picked is not None:
            self.selected.a[int(self.picked)] = True

    # Key and pointer bindings

    def button_press_event(self, widget, event):
        r"""Handle button press."""
        x = event.x
        y = event.y
        state = event.state
        self.pointer = [x, y]

        if event.button == 1 and not state & Gdk.ModifierType.CONTROL_MASK:
            if state & Gdk.ModifierType.SHIFT_MASK:
                self.srect = [x, y, x, y]
            elif self.picked == False:
                self.init_picked()
                self.queue_draw()
            if self.drag_begin is None:
                self.drag_begin = [x, y]
            return True

        if (event.button == 2 or
            (event.button == 1 and state & Gdk.ModifierType.CONTROL_MASK)):
            self.panning = (event.x, event.y)
            return True

        if event.button == 3:
            if isinstance(self.picked, PropertyMap):
                self.picked = None
                self.selected.fa = False
                self.queue_draw()
            elif self.picked is not False:
                self.picked = False
                self.selected.fa = False
                self.queue_draw()
            return True

    def button_release_event(self, widget, event):
        r"""Handle button release."""
        state = event.state
        if event.button == 1:
            if self.srect is not None:
                if self.picked == False:
                    self.init_picked()
                if not isinstance(self.picked, PropertyMap):
                    self.picked = self.selected

                if state & Gdk.ModifierType.CONTROL_MASK:
                    old_picked = self.picked.fa.copy()
                    self.picked.fa = False

                p1 = [self.srect[0], self.srect[1]]
                p2 = [self.srect[2], self.srect[3]]
                poly = [p1, [p1[0], p2[1]], p2, [p2[0], p1[1]]]
                poly = [self.pos_from_device(x) for x in poly]

                self.vertex_matrix.mark_polygon(poly, self.picked)

                if state & Gdk.ModifierType.CONTROL_MASK:
                    self.picked.fa = old_picked - self.picked.fa & old_picked

                self.srect = None

                self.queue_draw()
            self.drag_begin = None

            if self.moved_picked:
                self.moved_picked = False
                self.regenerate_surface(timeout=100)
                self.queue_draw()

            return True

        if event.button == 2:
            self.panning = None
            self.queue_draw()
            return True

    def motion_notify_event(self, widget, event):
        r"""Handle pointer motion."""
        if event.is_hint:
            x, y, state = event.window.get_pointer()[1:]
        else:
            x = event.x
            y = event.y
            state = event.state
        self.pointer = [x, y]

        if (state & Gdk.ModifierType.BUTTON1_MASK and
            not state & Gdk.ModifierType.CONTROL_MASK):
            if state & Gdk.ModifierType.SHIFT_MASK:
                if self.srect is not None:
                    self.srect[2:] = self.pointer
                    self.queue_draw()
            elif (self.picked is not None and self.picked is not False
                  and self.srect is None):
                p = self.pos_from_device(self.pointer)
                if isinstance(self.picked, PropertyMap):
                    if self.drag_begin is not None:
                        c = self.pos_from_device(self.drag_begin)
                        u = GraphView(self.g, vfilt=self.picked)
                        delta = np.asarray(p) - np.asarray(c)
                        for v in u.vertices():
                            new_pos = self.pos[v].a + delta
                            self.vertex_matrix.update_vertex(self.g.vertex(int(v)),
                                                             new_pos)
                        self.drag_begin = self.pointer
                elif self.vertex_matrix is not None:
                    self.vertex_matrix.update_vertex(self.picked, p)
                self.moved_picked = True
                self.queue_draw()
        elif (state & Gdk.ModifierType.BUTTON2_MASK or
              (state & Gdk.ModifierType.BUTTON1_MASK and
               state & Gdk.ModifierType.CONTROL_MASK)):
            if self.panning is not None:
                offset = [x - self.panning[0],
                          y - self.panning[1]]
                m = cairo.Matrix()
                m.translate(offset[0], offset[1])
                self.smatrix = self.smatrix * m
            self.panning = (x, y)
            self.queue_draw()
        else:
            self.panning = None

            if self.picked is not False:
                p = self.pos_from_device(self.pointer)
                v = self.vertex_matrix.get_closest(p)
                if v is not None and not isinstance(self.picked, PropertyMap):
                    if self.picked is not None:
                        self.selected[self.picked] = False
                        if self.picked != v:
                            self.queue_draw()
                    self.picked = v
                    self.selected[v] = True
        return True

    def scroll_event(self, widget, event):
        r"""Handle scrolling."""
        state = event.state

        angle = 0
        zoom = 1.

        if event.direction == Gdk.ScrollDirection.UP:
            if state & Gdk.ModifierType.CONTROL_MASK:
                angle = 0.1
            else:
                zoom = 1. / 0.9
                if state & Gdk.ModifierType.SHIFT_MASK:
                    scale_ink(1. / 0.9, self.vprops, self.eprops)
        elif event.direction == Gdk.ScrollDirection.DOWN:
            if state & Gdk.ModifierType.CONTROL_MASK:
                angle = -0.1
            else:
                zoom = 0.9
                if state & Gdk.ModifierType.SHIFT_MASK:
                    scale_ink(0.9, self.vprops, self.eprops)

        # keep centered
        if zoom != 1:
            center = self.pointer
            cpos = self.pos_from_device(center, surface=True)

            m = cairo.Matrix()
            m.scale(zoom, zoom)
            self.smatrix = self.smatrix.multiply(m)

            ncpos = self.pos_from_device(center, surface=True)
            self.smatrix.translate(ncpos[0] - cpos[0],
                                   ncpos[1] - cpos[1])
            self.regenerate_surface()
        if angle != 0:
            if not isinstance(self.picked, PropertyMap):
                center = (self.pointer[0], self.pointer[1])
                m = cairo.Matrix()
                m.translate(center[0], center[1])
                m.rotate(angle)
                m.translate(-center[0], -center[1])
                self.smatrix = self.smatrix.multiply(m)
                self.regenerate_surface()
            else:
                center = self.pos_from_device(self.pointer)
                u = GraphView(self.g, vfilt=self.picked)

                if self.vertex_matrix is not None:
                    for v in u.vertices():
                        self.vertex_matrix.remove_vertex(self.g.vertex(int(v)))

                m = cairo.Matrix()
                m.rotate(angle)
                m.translate(-center[0], -center[1])

                apply_transforms(u, self.pos, m)

                m = cairo.Matrix()
                m.translate(center[0], center[1])
                apply_transforms(u, self.pos, m)

                if self.vertex_matrix is not None:
                    for v in u.vertices():
                        self.vertex_matrix.add_vertex(self.g.vertex(int(v)))
                self.moved_picked = True

        self.queue_draw()
        return True

    def key_press_event(self, widget, event):
        r"""Handle key press."""

        #print event.keyval
        if event.keyval == 114:
            self.fit_to_window()
            self.regenerate_surface(timeout=50)
            self.queue_draw()
        elif event.keyval == 115:
            self.reset_layout()
        elif event.keyval == 97:
            self.apply_transform()
        elif event.keyval == 112:
            if self.picked == False:
                self.init_picked()
            else:
                self.picked = False
                self.selected.fa = False
                self.vertex_matrix = None
                self.queue_draw()
        elif event.keyval == 0x7a:
            if isinstance(self.picked, PropertyMap):
                u = GraphView(self.g, vfilt=self.picked)
                self.fit_to_window(g=u)
                self.regenerate_surface(timeout=50)
                self.queue_draw()
        return True

    def key_release_event(self, widget, event):
        r"""Handle release event."""
        #print event.keyval
        if event.keyval == 65507:
            if self.moved_picked:
                self.moved_picked = False
                self.regenerate_surface(timeout=100)
                self.queue_draw()
        return True


class GraphWindow(Gtk.Window):
    r"""Interactive GTK+ window containing a :class:`~graph_tool.draw.GraphWidget`.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be drawn.
    pos : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Vector-valued vertex property map containing the x and y coordinates of
        the vertices. If not given, it will be computed using :func:`sfdp_layout`.
    geometry : tuple
        Widget geometry.
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
    update_layout : bool (optional, default: ``True``)
        If ``True``, the layout will be updated dynamically.
    **kwargs
        Any extra parameters are passed to :class:`~graph_tool.draw.GraphWidget` and
        :func:`~graph_tool.draw.cairo_draw`.
    """

    def __init__(self, g, pos, geometry, vprops=None, eprops=None, vorder=None,
                 eorder=None, nodesfirst=False, update_layout=False, **kwargs):
        Gtk.Window.__init__(self, title="graph-tool's interactive window™")
        icon = GdkPixbuf.Pixbuf.new_from_file('%s/graph-tool-logo.svg' %
                                              os.path.dirname(__file__))
        self.set_icon(icon)
        self.set_default_size(geometry[0], geometry[1])

        self.graph = GraphWidget(g, pos, vprops, eprops, vorder, eorder,
                                 nodesfirst, update_layout, **kwargs)
        self.add(self.graph)

    def __del__(self):
        self.graph.cleanup()


def interactive_window(g, pos=None, vprops=None, eprops=None, vorder=None,
                       eorder=None, nodesfirst=False, geometry=(500, 400),
                       update_layout=True, async=False, **kwargs):
    r"""
    Display an interactive GTK+ window containing the given graph.

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
    geometry : tuple (optional, default: ``(500, 400)``)
        Window geometry.
    update_layout : bool (optional, default: ``True``)
        If ``True``, the layout will be updated dynamically.
    async : bool (optional, default: ``False``)
        If ``True``, run asynchronously. (Requires :mod:`IPython`)
    **kwargs
        Any extra parameters are passed to :class:`~graph_tool.draw.GraphWindow`,
        :class:`~graph_tool.draw.GraphWidget` and :func:`~graph_tool.draw.cairo_draw`.

    Returns
    -------
    pos : :class:`~graph_tool.PropertyMap`
        Vector vertex property map with the x and y coordinates of the vertices.
    selected : :class:`~graph_tool.PropertyMap` (optional, only if ``output is None``)
        Boolean-valued vertex property map marking the vertices which were
        selected interactively.

    Notes
    -----

    See documentation of :class:`~graph_tool.draw.GraphWidget` for key bindings
    information.

    """
    if pos is None:
        if update_layout:
            pos = random_layout(g, [1, 1])
        else:
            pos = sfdp_layout(g)
    win = GraphWindow(g, pos, geometry, vprops, eprops, vorder, eorder,
                      nodesfirst, update_layout, **kwargs)
    win.show_all()
    if async:
        # just a placeholder for a proper main loop integration with gtk3 when
        # ipython implements it
        import IPython.lib.inputhook
        f = lambda: Gtk.main_iteration_do(False)
        IPython.lib.inputhook.set_inputhook(f)
    else:
        win.connect("delete_event", Gtk.main_quit)
        Gtk.main()
    return pos, win.graph.selected.copy()


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
    >>> seed(42)
    >>> g = gt.random_graph(1000, lambda: min(zipf(2.4), 40),
    ...                     lambda i, j: exp(abs(i - j)), directed=False)
    >>> # extract largest component
    >>> g = gt.GraphView(g, vfilt=gt.label_largest_component(g))
    >>> deg = g.degree_property_map("out")
    >>> deg.a = 2 * (sqrt(deg.a) * 0.5 + 0.4)
    >>> ebet = gt.betweenness(g)[1]
    >>> ebet.a /= ebet.a.max() / 10.
    >>> gt.graph_draw(g, vertex_size=deg, vertex_fill_color=deg, vorder=deg,
    ...               edge_color=ebet, eorder=ebet, edge_pen_width=ebet,
    ...               output="graph-draw.pdf")
    <...>

    .. figure:: graph-draw.*
        :align: center

        SFDP force-directed layout of a graph with a power-law degree
        distribution, and dissortative degree correlation. The vertex size and
        color indicate the degree, and the edge color and width the edge
        betweeness centrality.

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
