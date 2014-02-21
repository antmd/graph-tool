## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = $(MOD_CPPFLAGS) $(CAIROMM_CFLAGS)

AM_CFLAGS = $(AM_CXXFLAGS)

libgraph_tool_spectraldir = $(MOD_DIR)/spectral

libgraph_tool_spectral_LTLIBRARIES = libgraph_tool_spectral.la

libgraph_tool_spectral_la_includedir = $(MOD_DIR)/include/spectral

libgraph_tool_spectral_la_LIBADD = $(MOD_LIBADD)

libgraph_tool_spectral_la_LDFLAGS = $(MOD_LDFLAGS)

libgraph_tool_spectral_la_SOURCES = \
    graph_adjacency.cc \
    graph_incidence.cc \
    graph_laplacian.cc \
    graph_norm_laplacian.cc \
    graph_matrix.cc \
    graph_transition.cc

libgraph_tool_spectral_la_include_HEADERS = \
    graph_adjacency.hh \
    graph_incidence.hh \
    graph_laplacian.hh \
    graph_transition.hh
