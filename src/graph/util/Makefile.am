## Process this file with automake to produce Makefile.in

AM_CPPFLAGS =\
    -I$(srcdir) -I$(srcdir)/.. \
    -I$(srcdir)/../boost-workaround \
    -DHAVE_CONFIG_H

AM_CFLAGS=$(AM_CXXFLAGS)

libgraph_tool_utildir = $(pythondir)/graph_tool/util

libgraph_tool_util_LTLIBRARIES = libgraph_tool_util.la

libgraph_tool_util_la_includedir = $(pythondir)/graph_tool/include

libgraph_tool_util_la_SOURCES = \
    graph_search.cc \
    graph_util_bind.cc

libgraph_tool_util_la_include_HEADERS = \
    graph_search.hh

libgraph_tool_util_la_LIBADD = \
    $(PYTHON_LDFLAGS) \
    $(BOOST_LDFLAGS) \
    $(OPENMP_LDFLAGS) \
    -lboost_python \
    -lboost_iostreams \
    -lexpat

# needed for typeinfo objects to work across DSO boundaries.
# see http://gcc.gnu.org/faq.html#dso
libgraph_tool_util_la_LDFLAGS = \
    -module \
    -avoid-version \
    -export-dynamic \
    -no-undefined \
    -Wl,-E
