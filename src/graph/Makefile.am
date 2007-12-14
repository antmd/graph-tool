## Process this file with automake to produce Makefile.in

AM_CPPFLAGS =\
    -I. -I.. \
    -I../boost-workaround \
    -DHAVE_CONFIG_H

AM_CXXFLAGS =\
    -Wall \
    $(PYTHON_CPPFLAGS) \
    $(BOOST_CPPFLAGS)

AM_CFLAGS=$(AM_CXXFLAGS)

libgraph_tool_ladir = ${pythondir}
libgraph_tool_includedir = $(includedir)/graph_tool
install-data-hook:		
	chmod 755 $(includedir)/graph_tool


lib_LTLIBRARIES = libgraph_tool.la

libgraph_tool_la_SOURCES = \
    graph.hh\
    graph.cc\
    graph_python_interface.cc\
    graph_properties.cc\
    graph_correlations.cc\
    graph_edge_correlations.cc\
    graph_correlations_combined.cc\
    graph_correlations_neighbours.cc\
    graph_assortativity.cc\
    graph_clustering.cc\
    graph_extended_clustering.cc\
    graph_generation.cc\
    graph_distance.cc\
    graph_distance_sampled.cc\
    graph_reciprocity.cc\
    graph_minimum_spanning_tree.cc\
    graph_community.cc\
    graph_community_network.cc\
    graph_line_graph.cc\
    graph_betweenness.cc\
    graph_rewiring.cc\
    graph_layout.cc\
    graph_io.cc\
    graph_bind.cc\
    graphml.cpp\
    read_graphviz_spirit.cpp\
    ../boost-workaround/boost/graph/filtered_graph.hpp\
    ../boost-workaround/boost/graph/fruchterman_reingold.hpp\
    ../boost-workaround/boost/graph/graphml.hpp

libgraph_tool_include_HEADERS = \
    graph_adaptor.hh\
    graph.hh\
    graph_filtering.hh\
    graph_python_interface.hh\
    graph_selectors.hh\
    graph_properties.hh\
    shared_map.hh\
    histogram.hh\
    ../../config.h

libgraph_tool_la_LIBADD = \
    $(PYTHON_LDFLAGS) \
    $(BOOST_LDFLAGS) \
    $(OPENMP_LDFLAGS) \
    -lboost_python \
    -lboost_iostreams \
    -lexpat

# needed for typeinfo objects to work across DSO boundaries.
# see http://gcc.gnu.org/faq.html#dso
libgraph_tool_la_LDFLAGS = \
    -Wl,-E 
