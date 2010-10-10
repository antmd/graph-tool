## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = $(MOD_CPPFLAGS)

AM_CXXFLAGS = $(CXXFLAGS) $(CGAL_FLAGS)

libgraph_tool_generationdir = $(pythondir)/graph_tool/generation

libgraph_tool_generation_LTLIBRARIES = libgraph_tool_generation.la

libgraph_tool_generation_la_includedir = $(pythondir)/graph_tool/include

libgraph_tool_generation_la_LIBADD = $(MOD_LIBADD) $(CGAL_LIBADD)

libgraph_tool_generation_la_LDFLAGS = $(MOD_LDFLAGS)

libgraph_tool_generation_la_SOURCES = \
    graph_generation.cc \
    graph_rewiring.cc \
    graph_predecessor.cc \
    graph_line_graph.cc \
    graph_union.cc \
    graph_union_vprop.cc \
    graph_union_eprop.cc \
    graph_triangulation.cc \
    graph_lattice.cc \
    graph_geometric.cc \
    graph_price.cc


libgraph_tool_generation_la_include_HEADERS = \
    graph_generation.hh \
    graph_rewiring.hh \
    graph_predecessor.hh \
    graph_union.hh \
    graph_triangulation.hh \
    graph_lattice.hh \
    graph_geometric.hh \
    graph_price.hh \
    sampler.hh
