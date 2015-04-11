#! /usr/bin/env python
# -*- coding: utf-8 -*-

# This simple example on how to interactive animations using graph-tool. Here we
# show the BFS tree of length 3 for the currently selected vertex.

from graph_tool.all import *

# We need some Gtk functions
from gi.repository import Gtk, Gdk
import sys

offscreen = sys.argv[1] == "offscreen" if len(sys.argv) > 1 else False

# We will use the network of network scientists, and filter out the largest
# component
g = collection.data["netscience"]
g = GraphView(g, vfilt=label_largest_component(g), directed=False)
g = Graph(g, prune=True)

pos = g.vp["pos"]  # layout positions

ecolor = g.new_edge_property("vector<double>")
for e in g.edges():
    ecolor[e] = [0.6, 0.6, 0.6, 1]
vcolor = g.new_vertex_property("vector<double>")
for v in g.vertices():
    vcolor[v] = [0.6, 0.6, 0.6, 1]

win = GraphWindow(g, pos, geometry=(500, 400),
                  edge_color=ecolor,
                  vertex_fill_color=vcolor)

orange = [0.807843137254902, 0.3607843137254902, 0.0, 1.0]
old_src = None
count = 0
def update_bfs(widget, event):
    global old_src, g, count, win
    src = widget.picked
    if src is None:
        return True
    if isinstance(src, PropertyMap):
        src = [v for v in g.vertices() if src[v]]
        if len(src) == 0:
            return True
        src = src[0]
    if src == old_src:
        return True
    old_src = src
    pred = shortest_distance(g, src, max_dist=3, pred_map=True)[1]
    for e in g.edges():
        ecolor[e] = [0.6, 0.6, 0.6, 1]
    for v in g.vertices():
        vcolor[v] = [0.6, 0.6, 0.6, 1]
    for v in g.vertices():
        w = g.vertex(pred[v])
        if w < g.num_vertices():
            e = g.edge(w, v)
            if e is not None:
                ecolor[e] = orange
                vcolor[v] = vcolor[w] = orange
    widget.regenerate_surface(lazy=False)
    widget.queue_draw()

    if offscreen:
        window = widget.get_window()
        pixbuf = Gdk.pixbuf_get_from_window(window, 0, 0, 500, 400)
        pixbuf.savev(r'./frames/bfs%06d.png' % count, 'png', [], [])
        count += 1

# Bind the function above as a montion notify handler
win.graph.connect("motion_notify_event", update_bfs)

# We will give the user the ability to stop the program by closing the window.
win.connect("delete_event", Gtk.main_quit)

# Actually show the window, and start the main loop.
win.show_all()
Gtk.main()
