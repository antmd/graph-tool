## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = $(MOD_CPPFLAGS)

AM_CFLAGS = $(AM_CXXFLAGS)

libgraph_tool_utildir = $(MOD_DIR)/util

libgraph_tool_util_LTLIBRARIES = libgraph_tool_util.la

libgraph_tool_util_la_includedir = $(MOD_DIR)/include/util

libgraph_tool_util_la_LIBADD = $(MOD_LIBADD)

libgraph_tool_util_la_LDFLAGS = $(MOD_LDFLAGS)

libgraph_tool_util_la_SOURCES = \
    graph_search.cc \
    graph_util_bind.cc

libgraph_tool_util_la_include_HEADERS = \
    graph_search.hh
