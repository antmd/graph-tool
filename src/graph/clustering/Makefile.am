## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = $(MOD_CPPFLAGS)

AM_CFLAGS = $(AM_CXXFLAGS)

libgraph_tool_clusteringdir = $(MOD_DIR)/clustering

libgraph_tool_clustering_LTLIBRARIES = libgraph_tool_clustering.la

libgraph_tool_clustering_la_includedir = $(MOD_DIR)/include/clustering

libgraph_tool_clustering_la_LIBADD = $(MOD_LIBADD)

libgraph_tool_clustering_la_LDFLAGS = $(MOD_LDFLAGS)

libgraph_tool_clustering_la_SOURCES = \
    graph_clustering.cc \
    graph_extended_clustering.cc \
    graph_motifs.cc

libgraph_tool_clustering_la_include_HEADERS = \
    graph_clustering.hh \
    graph_extended_clustering.hh \
    graph_motifs.hh

