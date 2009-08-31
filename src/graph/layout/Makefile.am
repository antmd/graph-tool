## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = $(MOD_CPPFLAGS)

AM_CFLAGS = $(AM_CXXFLAGS)

libgraph_tool_layoutdir = $(pythondir)/graph_tool/draw

libgraph_tool_layout_LTLIBRARIES = libgraph_tool_layout.la

libgraph_tool_layout_la_includedir = $(pythondir)/graph_tool/include

libgraph_tool_layout_la_LIBADD = $(MOD_LIBADD)

libgraph_tool_layout_la_LDFLAGS = $(MOD_LDFLAGS)

libgraph_tool_layout_la_SOURCES = \
    graph_arf.cc \
    graph_bind_layout.cc

libgraph_tool_layout_la_include_HEADERS = \
    graph_arf.hh
