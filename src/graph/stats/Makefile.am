## Process this file with automake to produce Makefile.in

AM_CPPFLAGS =\
    -I. -I.. \
    -I../boost-workaround \
    -DHAVE_CONFIG_H

AM_CFLAGS=$(AM_CXXFLAGS)

libgraph_tool_statsdir = $(pythondir)/graph_tool/stats

libgraph_tool_stats_LTLIBRARIES = libgraph_tool_stats.la

libgraph_tool_stats_la_includedir = $(pythondir)/graph_tool/include

libgraph_tool_stats_la_SOURCES = \
    graph_components.cc \
    graph_histograms.cc \
    graph_parallel.cc \
    graph_stats_bind.cc

libgraph_tool_stats_la_include_HEADERS = \
    graph_components.hh \
    graph_parallel.hh \
    graph_histograms.hh

libgraph_tool_stats_la_LIBADD = \
    $(PYTHON_LDFLAGS) \
    $(BOOST_LDFLAGS) \
    $(OPENMP_LDFLAGS) \
    -lboost_python \
    -lboost_iostreams \
    -lexpat

# needed for typeinfo objects to work across DSO boundaries.
# see http://gcc.gnu.org/faq.html#dso
libgraph_tool_stats_la_LDFLAGS = \
    -module \
    -avoid-version \
    -export-dynamic \
    -no-undefined \
    -Wl,-E
