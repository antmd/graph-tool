## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = $(MOD_CPPFLAGS)

AM_CFLAGS = $(AM_CXXFLAGS)

libgraph_tool_layoutdir = $(MOD_DIR)/draw

libgraph_tool_layout_LTLIBRARIES = libgraph_tool_layout.la

libgraph_tool_layout_la_includedir = $(MOD_DIR)/include/layout

libgraph_tool_layout_la_LIBADD = $(MOD_LIBADD)

libgraph_tool_layout_la_LDFLAGS = $(MOD_LDFLAGS)

libgraph_tool_layout_la_SOURCES = \
    graph_arf.cc \
    graph_fruchterman_reingold.cc \
    graph_sfdp.cc \
    graph_radial.cc \
    graph_bind_layout.cc

libgraph_tool_layout_la_include_HEADERS = \
    graph_arf.hh \
    graph_sfdp.hh
