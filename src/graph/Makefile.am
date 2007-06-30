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

libdir = ${pythondir}

lib_LTLIBRARIES = libgraph_tool.la

libgraph_tool_la_SOURCES = \
    graph_adaptor.hh\
    graph.hh\
    graph.cc\
    graph_filtering.hh\
    graph_python_filtering.hh\
    graph_selectors.hh\
    graph_properties.hh\
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
    graph_io.cc\
    graph_bind.cc\
    graphml.cpp\
    histogram.hh\
    shared_map.hh\
    read_graphviz_spirit.cpp

libgraph_tool_la_LIBADD = \
    $(PYTHON_LDFLAGS) \
    $(BOOST_LDFLAGS) \
    $(OPENMP_LDFLAGS) \
    -lboost_python \
    -lboost_iostreams \
    -lexpat

