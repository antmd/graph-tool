#! /usr/bin/env python
# -*- coding: utf-8 -*-

# This simple example on how to do animations using graph-tool. Here we do a
# simple simulation of an S->I->R->S epidemic model, where each vertex can be in
# one of the following states: Susceptible (S), infected (I), recovered (R). A
# vertex in the S state becomes infected either spontaneously with a probability
# 'x' or because a neighbour is infected. An infected node becomes recovered
# with probability 'r', and a recovered vertex becomes again susceptible with
# probability 's'.

# DISCLAIMER: The following code is definitely not the most efficient approach
# if you want to simulate this dynamics for very large networks, and/or for very
# long times. The main purpose is simply to highlight the animation capabilities
# of graph-tool.

from graph_tool.all import *
from numpy.random import *
import sys, os, os.path
import cairo

seed(42)
seed_rng(42)

# We need some Gtk and gobject functions
from gi.repository import Gtk, Gdk, GdkPixbuf
import gi._gobject as gobject

# We will use the karate-club network
g = collection.data["karate"]

pos = g.vp["pos"]  # layout positions

# We will filter out vertices which are in the "Recovered" state, by masking
# them using a property map.
removed = g.new_vertex_property("bool")

# SIRS dynamics parameters:

x = 0.001    # spontaneous outbreak probability
r = 0.1      # I->R probability
s = 0.01     # R->S probability

# (Note that the S->I transition happens simultaneously for every vertex with a
#  probability equal to the fraction of non-recovered neighbours which are
#  infected.)

S = 0
I = 1
R = 2

# Initialize all vertices to the S state
state = g.new_vertex_property("int")
state.a = S

# Images used to draw the nodes. They need to be loaded as cairo surfaces.
Simg = cairo.ImageSurface.create_from_png(open("face-grin.png", "rb"))
Simg_fear = cairo.ImageSurface.create_from_png(open("face-surprise.png", "rb"))
Iimg = cairo.ImageSurface.create_from_png(open("zombie.png", "rb"))

vertex_sfcs = g.new_vertex_property("object")
for v in g.vertices():
    vertex_sfcs[v] = Simg

# Newly infected nodes will be highlighted in red
newly_infected = g.new_vertex_property("bool")

# If True, the frames will be dumped to disk as images.
offscreen = sys.argv[1] == "offscreen" if len(sys.argv) > 1 else False
max_count = 500
if offscreen and not os.path.exists("./frames"):
    os.mkdir("./frames")

# This creates a GTK+ window with the initial graph layout
if not offscreen:
    win = GraphWindow(g, pos, geometry=(500, 400),
                      vertex_size=42,
                      vertex_anchor=0,
                      edge_color=[0.6, 0.6, 0.6, 1],
                      vertex_surface=vertex_sfcs,
                      vertex_halo=newly_infected,
                      vertex_halo_color=[0.8, 0, 0, 0.6])
else:
    count = 0
    win = Gtk.OffscreenWindow()
    win.set_default_size(500, 400)
    win.graph = GraphWidget(g, pos,
                            vertex_size=42,
                            vertex_anchor=0,
                            edge_color=[0.6, 0.6, 0.6, 1],
                            vertex_surface=vertex_sfcs,
                            vertex_halo=newly_infected,
                            vertex_halo_color=[0.8, 0, 0, 0.6])
    win.add(win.graph)


# This function will be called repeatedly by the GTK+ main loop, and we use it
# to update the state according to the SIRS dynamics.
def update_state():
    newly_infected.a = False
    removed.a = False

    # visit the nodes in random order
    vs = list(g.vertices())
    shuffle(vs)
    for v in vs:
        if state[v] == I:
            if random() < r:
                state[v] = R
        elif state[v] == S:
            if random() < x:
                state[v] = I
            else:
                ns = list(v.out_neighbours())
                if len(ns) > 0:
                    w = ns[randint(0, len(ns))]  # choose a random neighbour
                    if state[w] == I:
                        state[v] = I
                        newly_infected[v] = True
        elif random() < s:
            state[v] = S
        if state[v] == R:
            removed[v] = True

        if state[v] == S:
            if I in [state[w] for w in v.out_neighbours()]:
                vertex_sfcs[v] = Simg_fear
            else:
                vertex_sfcs[v] = Simg
        else:
            vertex_sfcs[v] = Iimg

    # Filter out the recovered vertices
    g.set_vertex_filter(removed, inverted=True)

    # The following will force the re-drawing of the graph, and issue a
    # re-drawing of the GTK window.
    win.graph.regenerate_surface(lazy=False)
    win.graph.queue_draw()

    # if doing an offscreen animation, dump frame to disk
    if offscreen:
        global count
        pixbuf = win.get_pixbuf()
        pixbuf.savev(r'./frames/zombies%06d.png' % count, 'png', [], [])
        if count > max_count:
            sys.exit(0)
        count += 1

    # We need to return True so that the main loop will call this function more
    # than once.
    return True


# Bind the function above as an 'idle' callback.
cid = gobject.idle_add(update_state)

# We will give the user the ability to stop the program by closing the window.
win.connect("delete_event", Gtk.main_quit)

# Actually show the window, and start the main loop.
win.show_all()
Gtk.main()
