## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = $(MOD_CPPFLAGS)

AM_CFLAGS=$(AM_CXXFLAGS)

libgraph_tool_searchdir = $(MOD_DIR)/search

libgraph_tool_search_LTLIBRARIES = libgraph_tool_search.la

libgraph_tool_search_la_includedir = $(MOD_DIR)/include/search

libgraph_tool_search_la_LIBADD = $(MOD_LIBADD)

libgraph_tool_search_la_LDFLAGS = $(MOD_LDFLAGS)

libgraph_tool_search_la_SOURCES = \
    graph_bfs.cc\
    graph_dfs.cc\
    graph_dijkstra.cc\
    graph_bellman_ford.cc\
    graph_astar.hh\
    graph_astar.cc\
    graph_astar_implicit.cc\
    graph_search_bind.cc

libgraph_tool_search_la_include_HEADERS = 
