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

libgraph_tool_correlationsdir = $(pythondir)/graph_tool

libgraph_tool_correlations_LTLIBRARIES = libgraph_tool_correlations.la

libgraph_tool_correlations_la_includedir = $(pythondir)/graph_tool/include

libgraph_tool_correlations_la_SOURCES = \
    graph_assortativity.cc \
    graph_correlations.cc \
    graph_correlations_combined.cc \
    graph_correlations_combined_imp1.cc \
    graph_correlations_combined_corr.cc \
    graph_correlations_combined_corr_imp1.cc \
    graph_correlations_imp1.cc \
    graph_correlations_imp2.cc \
    graph_correlations_imp3.cc \
    graph_correlations_neighbours.cc \
    graph_correlations_neighbours_imp1.cc \
    graph_correlations_neighbours_imp2.cc \
    graph_correlations_neighbours_imp3.cc \
    graph_correlations_neighbours_imp4.cc \
    graph_correlations_neighbours_imp5.cc \
    graph_correlations_neighbours_imp6.cc \
    graph_edge_correlations.cc \
    graph_edge_correlations_imp1.cc \
    graph_edge_correlations_imp2.cc \
    graph_edge_correlations_imp3.cc \
    graph_edge_correlations_imp4.cc \
    graph_edge_correlations_imp5.cc

libgraph_tool_correlations_la_include_HEADERS = \
    graph_assortativity.hh \
    graph_correlations_combined.hh \
    graph_correlations.hh \
    graph_correlations_neighbours.hh \
    graph_edge_correlations.hh

libgraph_tool_correlations_la_LIBADD = \
    $(PYTHON_LDFLAGS) \
    $(BOOST_LDFLAGS) \
    $(OPENMP_LDFLAGS) \
    -lboost_python \
    -lboost_iostreams \
    -lexpat

# needed for typeinfo objects to work across DSO boundaries.
# see http://gcc.gnu.org/faq.html#dso
libgraph_tool_correlations_la_LDFLAGS = \
    -module \
    -avoid-version \
    -export-dynamic \
    -no-undefined \
    -Wl,-E
