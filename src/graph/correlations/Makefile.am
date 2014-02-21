## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = $(MOD_CPPFLAGS)
AM_CFLAGS = $(AM_CXXFLAGS)

libgraph_tool_correlationsdir = $(MOD_DIR)/correlations

libgraph_tool_correlations_LTLIBRARIES = libgraph_tool_correlations.la

libgraph_tool_correlations_la_includedir = $(MOD_DIR)/include/correlations

libgraph_tool_correlations_la_LIBADD = $(MOD_LIBADD)

libgraph_tool_correlations_la_LDFLAGS = $(MOD_LDFLAGS)

libgraph_tool_correlations_la_SOURCES = \
    graph_assortativity.cc \
    graph_correlations.cc \
    graph_correlations_imp1.cc \
    graph_avg_correlations.cc \
    graph_avg_correlations_imp1.cc \
    graph_avg_correlations_combined.cc \
    graph_correlations_combined.cc \
    graph_correlations_bind.cc

libgraph_tool_correlations_la_include_HEADERS = \
    graph_assortativity.hh \
    graph_correlations.hh \
    graph_corr_hist.hh \
    graph_avg_correlations.hh

