## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = $(MOD_CPPFLAGS)

AM_CFLAGS=$(AM_CXXFLAGS)

libgraph_tool_centralitydir = $(MOD_DIR)/centrality

libgraph_tool_centrality_LTLIBRARIES = libgraph_tool_centrality.la

libgraph_tool_centrality_la_includedir = $(MOD_DIR)/include/centrality

libgraph_tool_centrality_la_LIBADD = $(MOD_LIBADD)

libgraph_tool_centrality_la_LDFLAGS = $(MOD_LDFLAGS)

libgraph_tool_centrality_la_SOURCES = \
    graph_betweenness.cc \
    graph_centrality_bind.cc \
    graph_closeness.cc \
    graph_eigentrust.cc \
    graph_eigenvector.cc \
    graph_hits.cc \
    graph_katz.cc \
    graph_pagerank.cc \
    graph_trust_transitivity.cc

libgraph_tool_centrality_la_include_HEADERS = \
    graph_closeness.hh \
    graph_eigentrust.hh \
    graph_eigenvector.hh \
    graph_pagerank.hh \
    graph_hits.hh \
    graph_katz.hh \
    graph_trust_transitivity.hh \
    minmax.hh
