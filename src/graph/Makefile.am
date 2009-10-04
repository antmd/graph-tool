## Process this file with automake to produce Makefile.in

SUBDIRS = generation stats clustering community util topology centrality correlations flow layout

AM_CPPFLAGS =\
    -I$(srcdir)/.. \
    -I$(srcdir)/../boost-workaround \
    -DHAVE_CONFIG_H

AM_CXXFLAGS =\
    -Wall \
    $(PYTHON_CPPFLAGS) \
    $(BOOST_CPPFLAGS)

AM_CFLAGS=$(AM_CXXFLAGS)

libgraph_tool_coredir = $(pythondir)/graph_tool
libgraph_tool_core_LTLIBRARIES = libgraph_tool_core.la

libgraph_tool_core_la_SOURCES = \
    graph_bind.cc \
    graph.cc \
    graph_copy.cc \
    graph_filtering.cc \
    graph_io.cc \
    graph_properties.cc \
    graph_properties_group.cc \
    graph_properties_ungroup.cc \
    graph_python_interface.cc \
    graph_python_interface_export.cc \
    graph_selectors.cc \
    graphml.cpp\
    read_graphviz_spirit.cpp\
    read_graphviz_new.cpp\
    ../boost-workaround/boost/graph/filtered_graph.hpp\
    ../boost-workaround/boost/graph/graphml.hpp\
    ../boost-workaround/boost/graph/betweenness_centrality.hpp\
    ../boost-workaround/boost/graph/copy.hpp\
    ../boost-workaround/boost/graph/isomorphism.hpp\
    ../boost-workaround/boost/graph/fruchterman_reingold.hpp\
    ../boost-workaround/boost/graph/kamada_kawai_spring_layout.hpp

libgraph_tool_core_la_includedir = $(pythondir)/graph_tool/include
libgraph_tool_core_la_include_HEADERS = \
    graph_adaptor.hh \
    graph_wrap.hh \
    graph_filtering.hh \
    graph.hh \
    graph_properties.hh \
    graph_properties_group.hh \
    graph_python_interface.hh \
    graph_selectors.hh \
    graph_util.hh \
    histogram.hh \
    mpl_nested_loop.hh \
    shared_map.hh \
    fast_vector_property_map.hh \
    numpy_bind.hh \
    ../../config.h

libgraph_tool_core_la_LIBADD = $(MOD_LIBADD)
libgraph_tool_core_la_LDFLAGS = $(MOD_LDFLAGS)

