## Process this file with automake to produce Makefile.in

AM_CPPFLAGS =\
    -I. -I.. \
    -I../boost-workaround \
    -DHAVE_CONFIG_H

AM_CFLAGS=$(AM_CXXFLAGS)

libgraph_tool_clusteringdir = $(pythondir)/graph_tool/clustering

libgraph_tool_clustering_LTLIBRARIES = libgraph_tool_clustering.la

libgraph_tool_clustering_la_includedir = $(pythondir)/graph_tool/include

libgraph_tool_clustering_la_SOURCES = \
    graph_clustering.cc \
    graph_extended_clustering.cc

libgraph_tool_clustering_la_include_HEADERS = \
    graph_clustering.hh \
    graph_extended_clustering.hh

libgraph_tool_clustering_la_LIBADD = \
    $(PYTHON_LDFLAGS) \
    $(BOOST_LDFLAGS) \
    $(OPENMP_LDFLAGS) \
    -lboost_python \
    -lboost_iostreams \
    -lexpat

# needed for typeinfo objects to work across DSO boundaries.
# see http://gcc.gnu.org/faq.html#dso
libgraph_tool_clustering_la_LDFLAGS = \
    -module \
    -avoid-version \
    -export-dynamic \
    -no-undefined \
    -Wl,-E
