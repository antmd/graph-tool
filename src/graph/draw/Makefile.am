## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = $(MOD_CPPFLAGS) $(CAIROMM_CFLAGS)

AM_CFLAGS = $(AM_CXXFLAGS)

libgraph_tool_drawdir = $(MOD_DIR)/draw

libgraph_tool_draw_LTLIBRARIES = libgraph_tool_draw.la

libgraph_tool_draw_la_includedir = $(MOD_DIR)/include/draw

libgraph_tool_draw_la_LIBADD = $(MOD_LIBADD)

libgraph_tool_draw_la_LDFLAGS = $(MOD_LDFLAGS) $(CAIROMM_LIBS)

libgraph_tool_draw_la_SOURCES = \
    graph_cairo_draw.cc \
    graph_tree_cts.cc

libgraph_tool_draw_la_include_HEADERS = 
