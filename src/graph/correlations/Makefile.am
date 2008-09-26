## Process this file with automake to produce Makefile.in

AM_CPPFLAGS =\
    -I. -I.. \
    -I../../boost-workaround \
    -DHAVE_CONFIG_H

AM_CFLAGS=$(AM_CXXFLAGS)

libgraph_tool_correlationsdir = $(pythondir)/graph_tool/correlations

libgraph_tool_correlations_LTLIBRARIES = libgraph_tool_correlations.la

libgraph_tool_correlations_la_includedir = $(pythondir)/graph_tool/include

libgraph_tool_correlations_la_SOURCES = \
    graph_assortativity.cc \
    graph_correlations.cc \
    graph_correlations_imp1.cc \
    graph_avg_correlations.cc \
    graph_avg_correlations_imp1.cc \
    graph_avg_correlations_combined.cc \
    graph_correlations_combined.cc \
    graph_correlations_bind.cc

libgraph_tool_correlations_la_include_HEADERS = \
    graph_assortativity.hh \
    graph_correlations.hh

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
