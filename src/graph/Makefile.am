## Process this file with automake to produce Makefile.in

SUBDIRS = centrality clustering community correlations draw flow generation layout search spectral stats topology util

AM_CPPFLAGS =\
    -I$(srcdir)/.. \
    -I$(srcdir)/../boost-workaround \
    -DHAVE_CONFIG_H

AM_CXXFLAGS =\
    -Wall \
    $(PYTHON_CPPFLAGS) \
    $(BOOST_CPPFLAGS)

AM_CFLAGS=$(AM_CXXFLAGS)

libgraph_tool_coredir = $(MOD_DIR)
libgraph_tool_core_LTLIBRARIES = libgraph_tool_core.la

libgraph_tool_core_la_SOURCES = \
    graph.cc \
    graph_exceptions.cc \
    graph_bind.cc \
    graph_copy.cc \
    graph_filtering.cc \
    graph_io.cc \
    graph_openmp.cc \
    graph_properties.cc \
    graph_properties_copy.cc \
    graph_properties_group.cc \
    graph_properties_ungroup.cc \
    graph_python_interface.cc \
    graph_python_interface_export.cc \
    graph_selectors.cc \
    graphml.cpp \
    random.cc \
    read_graphviz_new.cpp

libgraph_tool_core_la_includedir = $(MOD_DIR)/include
libgraph_tool_core_la_include_HEADERS = \
    ../../config.h \
    fast_vector_property_map.hh \
    gml.hh \
    graph.hh \
    graph_adjacency.hh \
    graph_adaptor.hh \
    graph_exceptions.hh \
    graph_filtering.hh \
    graph_io_binary.hh \
    graph_properties.hh \
    graph_properties_group.hh \
    graph_python_interface.hh \
    graph_selectors.hh \
    graph_util.hh \
    histogram.hh \
    mpl_nested_loop.hh \
    numpy_bind.hh \
    numpy_bind_old.hh \
    random.hh \
    str_repr.hh \
    shared_map.hh \
    transform_iterator.hh

libgraph_tool_core_la_workarounddir = $(MOD_DIR)/include/boost-workaround/boost/graph/
libgraph_tool_core_la_workaround_HEADERS = \
    ../boost-workaround/boost/graph/betweenness_centrality.hpp \
    ../boost-workaround/boost/graph/detail/read_graphviz_new.hpp \
    ../boost-workaround/boost/graph/graphml.hpp \
    ../boost-workaround/boost/graph/graphviz.hpp \
    ../boost-workaround/boost/graph/isomorphism.hpp \
    ../boost-workaround/boost/graph/kamada_kawai_spring_layout.hpp \
    ../boost-workaround/boost/graph/metric_tsp_approx.hpp \
    ../boost-workaround/boost/graph/named_function_params-alt.hpp \
    ../boost-workaround/boost/graph/overloading.hpp \
    ../boost-workaround/boost/graph/push_relabel_max_flow.hpp \
    ../boost-workaround/boost/graph/copy_alt.hpp \
    ../boost-workaround/boost/graph/reverse_graph_alt.hpp \
    ../boost-workaround/boost/graph/stoer_wagner_min_cut.hpp

libgraph_tool_core_la_LIBADD = $(MOD_LIBADD)
libgraph_tool_core_la_LDFLAGS = $(MOD_LDFLAGS)

