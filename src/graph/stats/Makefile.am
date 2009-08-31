## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = $(MOD_CPPFLAGS)

AM_CXXFLAGS = $(CXXFLAGS)

libgraph_tool_statsdir = $(pythondir)/graph_tool/stats

libgraph_tool_stats_LTLIBRARIES = libgraph_tool_stats.la

libgraph_tool_stats_la_includedir = $(pythondir)/graph_tool/include

libgraph_tool_stats_la_SOURCES = \
    graph_histograms.cc \
    graph_average.cc \
    graph_parallel.cc \
    graph_distance.cc \
    graph_distance_sampled.cc \
    graph_stats_bind.cc


libgraph_tool_stats_la_include_HEADERS = \
    graph_parallel.hh \
    graph_histograms.hh \
    graph_average.hh \
    graph_distance_sampled.hh \
    graph_distance.hh

libgraph_tool_stats_la_LIBADD = $(MOD_LIBADD)

libgraph_tool_stats_la_LDFLAGS = $(MOD_LDFLAGS)
