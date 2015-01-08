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

from .. import GraphView, PropertyMap, ungroup_vector_property,\
     group_vector_property, _prop
from .cairo_draw import *
from .cairo_draw import _vdefaults, _edefaults
from .. draw import sfdp_layout, random_layout, _avg_edge_distance, \
    coarse_graphs


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
            for i in range(-n, n):
                for j in range(-n, n):
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
        for i in range(-1, 2):
            for j in range(-1, 2):
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
        for i in range(p1[0], p2[0] + 1):
            for j in range(p1[1], p2[1] + 1):
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
    def __init__(self, g, pos, vprops=None, eprops=None, vorder=None,
                 eorder=None, nodesfirst=False, update_layout=False,
                 layout_K=1., multilevel=False, display_props=None,
                 display_props_size=11, bg_color=None, layout_callback=None,
                 key_press_callback=None, **kwargs):
        r"""Interactive GTK+ widget displaying a given graph.

        Parameters
        ----------
        g : :class:`~graph_tool.Graph`
            Graph to be drawn.
        pos : :class:`~graph_tool.PropertyMap`
            Vector-valued vertex property map containing the x and y coordinates of
            the vertices.
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
        layout_callback : function (optional, default: ``Node``)
            User-supplied callback to be called whenever the positions of the layout
            have changed. It needs to have the following signature:

            .. code-block::

               def callback(g, picked, pos, vprops, eprops):
                   ...

            where ``g`` is the graph being drawn, ``picked`` is either a single
            vertex or a boolean vertex property map representing the vertices
            currently selected, and ``vprops`` and ``eprops`` are dictionaries with
            the vertex and edge properties currently being used by the layout.
        key_press_callback : function (optional, default: ``Node``)

            User-supplied callback to be called whenever a key-press event has
            happened. It needs to have the following signature:

            .. code-block::

               def callback(g, keyval, picked, pos, vprops, eprops):
                   ...

            where ``g`` is the graph being drawn, ``keyval`` is the key id,
            ``picked`` is either a single vertex or a boolean vertex property map
            representing the vertices currently selected, and ``vprops`` and
            ``eprops`` are dictionaries with the vertex and edge properties
            currently being used by the layout.
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
        activated. Vertices which are currently selected are not updated.

        """

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
        self.base_geometry = None
        self.background = None
        self.bg_color = bg_color if bg_color is not None else [1, 1, 1, 1]
        self.surface_callback = None

        self.layout_callback_id = None
        self.layout_K = layout_K
        self.layout_init_step = self.layout_K
        self.epsilon = 0.01 * self.layout_K
        self.multilevel_layout = multilevel

        self.layout_user_callback = layout_callback
        self.key_press_user_callback = key_press_callback

        if multilevel:
            self.cgs = coarse_graphs(g)
            u = next(self.cgs)
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
        if self.layout_callback_id is None or self.g.num_vertices() == 0:
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

        if self.layout_user_callback is not None:
            self.layout_user_callback(self.g, self.picked, self.pos, self.vprops,
                                      self.eprops)

        if delta > self.epsilon:
            return True
        else:
            if self.multilevel_layout:
                try:
                    u = next(self.cgs)
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
            if (self.base is None or self.base_geometry[0] != geometry[0] or
                self.base_geometry[1] != geometry[1]):
                # self.base = cairo.ImageSurface(cairo.FORMAT_ARGB32,
                #                                *geometry)
                w = self.get_window()
                if w is None:
                    return False
                self.base = w.create_similar_surface(cairo.CONTENT_COLOR_ALPHA,
                                                     *geometry)
                self.base_geometry = geometry

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

        ul = self.pos_to_device((0, 0), surface=True, cr=cr)
        lr = self.pos_to_device(self.base_geometry, surface=True, cr=cr)

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
        ox, oy = self.get_window().get_position()
        if cr is None:
            cr = self.get_window().cairo_create()
            if surface:
                cr.set_matrix(self.smatrix)
            else:
                cr.set_matrix(self.tmatrix * self.smatrix)
        if dist:
            return cr.user_to_device_distance(pos[0], pos[1])
        else:
            x, y = cr.user_to_device(pos[0], pos[1])
            return (x - ox, y - oy)

    def pos_from_device(self, pos, dist=False, surface=False, cr=None):
        """Convert a position from the widget space to the device space."""
        ox, oy = self.get_window().get_position()
        if cr is None:
            cr = self.get_window().cairo_create()
            if surface:
                cr.set_matrix(self.smatrix)
            else:
                cr.set_matrix(self.tmatrix * self.smatrix)
        if dist:
            return cr.device_to_user_distance(pos[0], pos[1])
        else:
            return cr.device_to_user(pos[0] + ox, pos[1] + oy)

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
        ox, oy = self.get_window().get_position()
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
        m.translate(offset[0] + ox, offset[1] + oy)
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
                if self.layout_user_callback is not None:
                    self.layout_user_callback(self.g, self.picked, self.pos,
                                              self.vprops, self.eprops)
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

        if self.key_press_user_callback is not None:
            self.key_press_user_callback(self.g, event.keyval, self.picked,
                                         self.pos, self.vprops, self.eprops)

        if event.keyval == 65507:
            if self.moved_picked:
                self.moved_picked = False
                self.regenerate_surface(timeout=100)
                self.queue_draw()
        return True


class GraphWindow(Gtk.Window):
    def __init__(self, g, pos, geometry, vprops=None, eprops=None, vorder=None,
                 eorder=None, nodesfirst=False, update_layout=False, **kwargs):
        r"""Interactive GTK+ window containing a :class:`~graph_tool.draw.GraphWidget`.

        Parameters
        ----------
        g : :class:`~graph_tool.Graph`
            Graph to be drawn.
        pos : :class:`~graph_tool.PropertyMap`
            Vector-valued vertex property map containing the x and y coordinates of
            the vertices.
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


_window_list = []

def interactive_window(g, pos=None, vprops=None, eprops=None, vorder=None,
                       eorder=None, nodesfirst=False, geometry=(500, 400),
                       update_layout=True, async=False, no_main=False, **kwargs):
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
    no_main : bool (optional, default: ``False``)
        If ``True``, the GTK+ main loop will not be called.
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
    _window_list.append(win)
    if not no_main:
        if async:
            # just a placeholder for a proper main loop integration with gtk3 when
            # ipython implements it
            import IPython.lib.inputhook
            f = lambda: Gtk.main_iteration_do(False)
            IPython.lib.inputhook.set_inputhook(f)
        else:
            def destroy_callback(*args, **kwargs):
                global _window_list
                for w in _window_list:
                    w.destroy()
                Gtk.main_quit()
            win.connect("delete_event", destroy_callback)
            Gtk.main()
    return pos, win.graph.selected.copy()
