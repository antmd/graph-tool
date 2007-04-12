## Process this file with automake to produce Makefile.in

AUTOMAKE_OPTIONS = subdir-objects

SUBDIRS = graph .

bin_SCRIPTS = graph-tool

boost_workarounddir = boost-workaround
boost_workaround_DATA = \
    boost-workaround/boost/graph/filtered_graph.hpp \
    boost-workaround/boost/graph/fruchterman_reingold.hpp \
    boost-workaround/boost/graph/graphml.hpp


EXTRA_DIST = $(bin_SCRIPTS) $(boost_workaround_DATA)