#! /usr/bin/env python

# We probably will need some things from several places
import sys, os
from pylab import * # for plotting
from numpy.random import * # for random sampling
seed(42)

# We need to import the graph_tool module itself
from graph_tool.all import *

# let's construct a Price network (the one that existed before Barabasi). It is
# a directed network, with preferential attachment. The algorithm below is a
# very naive, and a bit slow, but quite simple.

# We start with an empty, directed graph
g = Graph()

# We want also to keep the age information for each vertex and edge. For that
# let's create some property maps
v_age = g.new_vertex_property("int")
e_age = g.new_edge_property("int")

# The final size of the network
N = 100000

# We have to start with one vertex
v = g.add_vertex()
v_age[v] = 0

# we will keep a list of the vertices. The number of times a vertex is in this
# list will give the probability of it being selected.
vlist = [v]

# let's now add the new edges and vertices
for i in xrange(1, N):
    # create our new vertex
    v = g.add_vertex()
    v_age[v] = i

    # we need to sample a new vertex to be the target, based on its in-degree +
    # 1. For that, we simply randomly sample it from vlist.
    i = randint(0, len(vlist))
    target = vlist[i]

    # add edge
    e = g.add_edge(v, target)
    e_age[e] = i

    # put v and target in the list
    vlist.append(target)
    vlist.append(v)

# now we have a graph!
# let's calculate its in-degree distribution and plot it to a file

in_hist = vertex_hist(g, "in")

clf()
errorbar(in_hist[1], in_hist[0], fmt="o", yerr=sqrt(in_hist[0]), label="in")
gca().set_yscale("log")
gca().set_xscale("log")
legend(loc="best")
xlabel("$k_{in}$")
ylabel("$P(k_{in})$")
savefig("deg-hist.png")

# let's do a random walk on the graph and print the age of the vertices we find,
# just for fun.

v = g.vertex(randint(0, g.num_vertices()))
for i in xrange(0, 100):
    print "vertex:", v, "in-degree:", v.in_degree(), "out-degree:",\
          v.out_degree(), "age:", v_age[v]

    if v.out_degree() == 0:
        print "Nowhere else to go... We found the main hub!"
        break

    n_list = []
    for w in v.out_neighbours():
        n_list.append(w)
    v = n_list[randint(0, len(n_list))]

# let's save our graph for posterity We want to save the age properties as
# well... To do this, they must become "internal" properties, as such:

g.vertex_properties["age"] = v_age
g.edge_properties["age"] = e_age

# now we can save it
g.save("price.xml.gz")

