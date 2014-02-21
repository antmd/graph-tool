## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = $(MOD_CPPFLAGS)

AM_CFLAGS = $(AM_CXXFLAGS)

libgraph_tool_flowdir = $(MOD_DIR)/flow

libgraph_tool_flow_LTLIBRARIES = libgraph_tool_flow.la

libgraph_tool_flow_la_includedir = $(MOD_DIR)/include/flow

libgraph_tool_flow_la_LIBADD = $(MOD_LIBADD)

libgraph_tool_flow_la_LDFLAGS = $(MOD_LDFLAGS)

libgraph_tool_flow_la_SOURCES = \
    graph_edmonds_karp.cc \
    graph_push_relabel.cc \
    graph_kolmogorov.cc \
    graph_maximum_cardinality_matching.cc \
    graph_minimum_cut.cc \
    graph_flow_bind.cc

libgraph_tool_flow_la_include_HEADERS = \
    graph_augment.hh
