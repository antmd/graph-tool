## Process this file with automake to produce Makefile.in

#SUBDIRS = correlations

AM_CPPFLAGS =\
    -I. -I.. \
    -I../boost-workaround \
    -I $(pythondir)/numpy/core/include/numpy/ \
    -DHAVE_CONFIG_H

AM_CXXFLAGS =\
    -Wall \
    $(PYTHON_CPPFLAGS) \
    $(BOOST_CPPFLAGS)

AM_CFLAGS=$(AM_CXXFLAGS)

libgraph_tool_coredir = $(pythondir)/graph_tool

libgraph_tool_core_LTLIBRARIES = libgraph_tool_core.la

libgraph_tool_core_la_includedir = $(pythondir)/graph_tool/include

libgraph_tool_core_la_SOURCES = \
    graph_bind.cc \
    graph.cc \
    graph_copy.cc \
    graph_filtering.cc \
    graph_io.cc \
    graph_properties.cc \
    graph_python_interface.cc \
    graph_python_interface_export.cc \
    graph_selectors.cc \
    graphml.cpp\
    read_graphviz_spirit.cpp\
    ../boost-workaround/boost/graph/filtered_graph.hpp\
    ../boost-workaround/boost/graph/graphml.hpp

libgraph_tool_core_la_include_HEADERS = \
    graph_adaptor.hh \
    graph_filtering.hh \
    graph.hh \
    graph_properties.hh \
    graph_python_interface.hh \
    graph_selectors.hh \
    graph_util.hh \
    histogram.hh \
    mpl_nested_loop.hh \
    shared_map.hh \
    numpy_bind.hh \
    ../../config.h

# leave out precompiled headers for now, since it doesn't seem to help much
# BUILT_SOURCES = \
#     graph_filtering.hh.gch

# ## Header precompilation
# ## FIXME: need a better way to convince libtool to let us do this.
# $(libgraph_tool_la_include_HEADERS):

# $(BUILT_SOURCES):
# 	$(CXXCOMPILE) $(CXXFLAGS) $(AM_CXXFLAGS) $(AM_CPPFLAGS) -fPIC -DPIC -x c++-header `basename $@ .gch`
# mostlyclean-local:
# 	-rm -f *.gch

libgraph_tool_core_la_LIBADD = \
    $(PYTHON_LDFLAGS) \
    $(BOOST_LDFLAGS) \
    $(OPENMP_LDFLAGS) \
    -lboost_python \
    -lboost_iostreams \
    -lexpat

# needed for typeinfo objects to work across DSO boundaries.
# see http://gcc.gnu.org/faq.html#dso
libgraph_tool_core_la_LDFLAGS = \
    -module \
    -avoid-version \
    -export-dynamic \
    -no-undefined \
    -Wl,-E

