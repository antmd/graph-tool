## Process this file with automake to produce Makefile.in

AM_CPPFLAGS =\
    -I$(srcdir) -I$(srcdir)/.. \
    -I$(srcdir)/../boost-workaround \
    -DHAVE_CONFIG_H

AM_CFLAGS=$(AM_CXXFLAGS)

libgraph_tool_miscdir = $(pythondir)/graph_tool/misc

libgraph_tool_misc_LTLIBRARIES = libgraph_tool_misc.la

libgraph_tool_misc_la_includedir = $(pythondir)/graph_tool/include

libgraph_tool_misc_la_SOURCES = \
    graph_rewiring.cc \
    graph_misc.cc

libgraph_tool_misc_la_include_HEADERS = \
    graph_rewiring.hh

libgraph_tool_misc_la_LIBADD = \
    $(PYTHON_LDFLAGS) \
    $(BOOST_LDFLAGS) \
    $(OPENMP_LDFLAGS) \
    -lboost_python \
    -lboost_iostreams \
    -lexpat

# needed for typeinfo objects to work across DSO boundaries.
# see http://gcc.gnu.org/faq.html#dso
libgraph_tool_misc_la_LDFLAGS = \
    -module \
    -avoid-version \
    -export-dynamic \
    -no-undefined \
    -Wl,-E
