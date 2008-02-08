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

libgraph_tooldir = $(pythondir)/graph_tool

libgraph_tool_LTLIBRARIES = libgraph_tool.la

libgraph_tool_la_includedir = $(pythondir)/graph_tool/include

libgraph_tool_la_SOURCES = \
    graph_assortativity.cc \
    graph_betweenness.cc \
    graph_bind.cc \
    graph.cc \
    graph_clustering.cc \
    graph_community.cc \
    graph_community_network.cc \
    graph_correlations.cc \
    graph_correlations_combined.cc \
    graph_correlations_combined_corr.cc \
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
    graph_distance.cc \
    graph_distance_sampled.cc \
    graph_edge_correlations.cc \
    graph_edge_correlations_imp1.cc \
    graph_edge_correlations_imp2.cc \
    graph_edge_correlations_imp3.cc \
    graph_edge_correlations_imp4.cc \
    graph_edge_correlations_imp5.cc \
    graph_extended_clustering.cc \
    graph_filtering.cc \
    graph_generation.cc \
    graph_io.cc \
    graph_layout.cc \
    graph_line_graph.cc \
    graph_minimum_spanning_tree.cc \
    graph_properties.cc \
    graph_python_interface.cc \
    graph_python_interface_export.cc \
    graph_reciprocity.cc \
    graph_rewiring.cc \
    graph_selectors.cc \
    graphml.cpp\
    read_graphviz_spirit.cpp\
    ../boost-workaround/boost/graph/filtered_graph.hpp\
    ../boost-workaround/boost/graph/fruchterman_reingold.hpp\
    ../boost-workaround/boost/graph/graphml.hpp

libgraph_tool_la_include_HEADERS = \
    graph_adaptor.hh \
    graph_assortativity.hh \
    graph_clustering.hh \
    graph_community.hh \
    graph_community_network.hh \
    graph_correlations_combined.hh \
    graph_correlations.hh \
    graph_correlations_neighbours.hh \
    graph_distance.hh \
    graph_distance_sampled.hh \
    graph_edge_correlations.hh \
    graph_extended_clustering.hh \
    graph_filtering.hh \
    graph.hh \
    graph_properties.hh \
    graph_python_interface.hh \
    graph_rewiring.hh \
    graph_selectors.hh \
    graph_util.hh \
    histogram.hh \
    mpl_nested_loop.hh \
    shared_map.hh \
    ../../config.h

BUILT_SOURCES = \
    graph_filtering.hh.gch

## Header precompilation
## FIXME: need a better way to convince libtool to let us do this.
$(libgraph_tool_la_include_HEADERS):

$(BUILT_SOURCES):
	$(CXXCOMPILE) $(CXXFLAGS) $(AM_CXXFLAGS) $(AM_CPPFLAGS) -fPIC -DPIC -x c++-header `basename $@ .gch`
mostlyclean-local:
	-rm -f *.gch

libgraph_tool_la_LIBADD = \
    $(PYTHON_LDFLAGS) \
    $(BOOST_LDFLAGS) \
    $(OPENMP_LDFLAGS) \
    -lboost_python \
    -lboost_iostreams \
    -lexpat

# needed for typeinfo objects to work across DSO boundaries.
# see http://gcc.gnu.org/faq.html#dso
libgraph_tool_la_LDFLAGS = \
    -module \
    -avoid-version \
    -export-dynamic \
    -no-undefined \
    -Wl,-E

