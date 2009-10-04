## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = $(MOD_CPPFLAGS)

AM_CFLAGS=$(AM_CXXFLAGS)

libgraph_tool_communitydir = $(pythondir)/graph_tool/community

libgraph_tool_community_LTLIBRARIES = libgraph_tool_community.la

libgraph_tool_community_la_includedir = $(pythondir)/graph_tool/include

libgraph_tool_community_la_LIBADD = $(MOD_LIBADD)

libgraph_tool_community_la_LDFLAGS = $(MOD_LDFLAGS)

libgraph_tool_community_la_SOURCES = \
    graph_community.cc \
    graph_community_network.cc

libgraph_tool_community_la_include_HEADERS = \
    graph_community.hh \
    graph_community_network.hh
