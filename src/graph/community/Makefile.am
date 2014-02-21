## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = $(MOD_CPPFLAGS)

AM_CFLAGS=$(AM_CXXFLAGS)

libgraph_tool_communitydir = $(MOD_DIR)/community

libgraph_tool_community_LTLIBRARIES = libgraph_tool_community.la

libgraph_tool_community_la_includedir = $(MOD_DIR)/include/community

libgraph_tool_community_la_LIBADD = $(MOD_LIBADD)

libgraph_tool_community_la_LDFLAGS = $(MOD_LDFLAGS)

libgraph_tool_community_la_SOURCES = \
    graph_blockmodel.cc \
    graph_blockmodel_overlap.cc \
    graph_community.cc \
    graph_community_network.cc \
    graph_community_network_edges.cc \
    graph_community_network_vavg.cc \
    graph_community_network_eavg.cc \
    graph_community_network_eavg_imp1.cc \
    spence.cc

libgraph_tool_community_la_include_HEADERS = \
    graph_blockmodel.hh \
    graph_blockmodel_overlap.hh \
    graph_community.hh \
    graph_community_network.hh
