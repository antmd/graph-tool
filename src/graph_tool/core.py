#! /usr/bin/env python
# graph_tool.py -- a general graph manipulation python module
#
# Copyright (C) 2007 Tiago de Paula Peixoto <tiago@forked.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
# RTLD_GLOBAL needs to be set in dlopen() if we want typeinfo and friends to
# work properly across DSO boundaries. See http://gcc.gnu.org/faq.html#dso

# The "except" is because the dl module raises a system error on ia64 and x86_64
# systems because "int" and addresses are different sizes.
try:
    from dl import RTLD_LAZY, RTLD_NOW, RTLD_GLOBAL
except ImportError:
    RTLD_LAZY = 1
    RTLD_NOW = 2
    RTLD_GLOBAL = 256
_orig_dlopen_flags = sys.getdlopenflags()

sys.setdlopenflags(RTLD_LAZY|RTLD_GLOBAL)
import libgraph_tool_core
sys.setdlopenflags(_orig_dlopen_flags) # reset it to normal case to avoid
                                       # unnecessary symbol collision
__version__ = libgraph_tool_core.mod_info().version

import os, os.path, re, struct, fcntl, termios, gzip, bz2, string,\
       textwrap, time, signal, traceback, shutil, time, math, inspect, \
       functools, types

################################################################################
# Utility functions
################################################################################

def _degree(name):
    """Retrieve the degree type from string"""
    deg = name
    if name == "in-degree" or name == "in":
        deg = libgraph_tool_core.Degree.In
    if name == "out-degree" or name == "out":
        deg = libgraph_tool_core.Degree.Out
    if name == "total-degree" or name == "total":
        deg = libgraph_tool_core.Degree.Total
    return deg

def _parse_range(range):
    """Parse a range in the form 'lower[*] upper[*]' (where an '*' indicates
    a closed boundary) or 'comp val', where 'comp' is [==|=|!=|>=|<=|>|<].
    """
    try:
        rcomp = re.compile(r"(>|<|>=|<=|==|=|!=)\s*([^=<>!]+)")
        rrange = re.compile(r"([^\*]+)(\*?)\s+([^\*]+)(\*?)")
        inverted = False
        if rcomp.match(range):
            comp = rcomp.match(range).group(1)
            thres = rcomp.match(range).group(2)
            thres = float(thres)
            if "<" in comp:
                ran = (float('-Inf'), thres)
                if "=" in comp:
                    inc = (False, True)
                else:
                    inc = (False, False)
            elif ">" in comp:
                ran = (thres, float('Inf'))
                if "=" in comp:
                    inc = (True, False)
                else:
                    inc = (False, False)
            elif comp == "==" or comp == "=" or comp == "!=":
                ran = (thres, thres)
                inc = (True, True)
                if comp == "!=":
                    inverted = True
        elif rrange.match(range):
            m = rrange.match(range)
            r1 = float(m.group(1))
            i1 = "*" in m.group(2)
            r2 = float(m.group(3))
            i2 = "*" in m.group(4)
            ran = (r1, r2)
            inc = (i1, i2)
        else:
            raise ValueError("invalid value for range: " + range)
    except (ValueError, TypeError), e:
        raise ValueError("invalid value for range: %s: %s " % (range, str(e)))
    return (ran, inc, inverted)

################################################################################
# Decorators
# Some useful function decorators which will be used below
################################################################################

def _wraps(func):
    """This decorator works like the functools.wraps meta-decorator, but
    also preserves the function's argument signature. This uses eval, and is
    thus a bit of a hack, but there no better way I know of to do this."""
    def decorate(f):
        argspec = inspect.getargspec(func)
        args_call = inspect.formatargspec(argspec[0])
        argspec = inspect.formatargspec(argspec[0], defaults=argspec[3])
        argspec = argspec.lstrip("(").rstrip(")")
        wrap = eval("lambda %s: f%s" % (argspec, args_call), locals())
        return functools.wraps(func)(wrap)
    return decorate

def _attrs(**kwds):
    """Decorator which adds arbitrary attributes to methods"""
    def decorate(f):
        for k in kwds:
            setattr(f, k, kwds[k])
        return f
    return decorate

def _handle_exceptions(func):
    """Decorator which will catch and properly propagate exceptions raised
    by GraphInterface in all the methods"""
    @_wraps(func)
    def wrap(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except (IOError, RuntimeError), e:
            libgraph_tool_core.raise_error(str(e))
    return wrap

def _limit_args(allowed_vals):
    """Decorator which will limit the values of given arguments to a specified
    list of allowed values, and raise TypeError exception if the given value
    doesn't belong. 'allowed_vals' is a dict containing the allowed value list
    for each limited function argument."""
    def decorate(func):
        @_wraps(func)
        def wrap(*args, **kwargs):
            arg_names = inspect.getargspec(func)[0]
            arguments = zip(arg_names, args)
            arguments += [(k, kwargs[k]) for k in kwargs.keys()]
            for a in arguments:
                if allowed_vals.has_key(a[0]):
                    if a[1] not in allowed_vals[a[0]]:
                        raise TypeError("value for '%s' must be one of: %s" % \
                                         (a[0], ", ".join(allowed_vals[a[0]])))
            return func(*args, **kwargs)
        return wrap
    return decorate

def _lazy_load(func):
    """Decorator which will call the 'load' method before executing the
    decorated function, thus implementing 'lazy loading'."""
    @_wraps(func)
    def wrap(*args, **kwargs):
        self = args[0]
        if self.lazy_filename != None and self.lazy_format != None:
            self.load(self.lazy_filename, self.lazy_format)
            self.lazy_filename = None
            self.lazy_format = None
        return func(*args, **kwargs)
    return wrap

################################################################################
# Graph class
# The main graph interface
################################################################################

from libgraph_tool_core import GraphError, Vector_bool, Vector_int32_t, \
     Vector_int64_t, Vector_double, Vector_long_double

class Graph(object):
    """The main graph type which encapsulates the GraphInterface"""

    def __init__(self, g = None):
        if g == None:
            self.__graph = libgraph_tool_core.GraphInterface()
            self.lazy_filename = None
            self.lazy_format = None
        else:
            self.__graph = libgraph_tool_core.GraphInterface(g.__graph)
            self.lazy_filename = g.lazy_filename
            self.lazy_format = g.lazy_format

    @_handle_exceptions
    @_lazy_load
    def copy(self):
        """Returns a deep copy of self"""
        new_graph = Graph()
        new_graph.__graph = libgraph_tool_core.GraphInterface(self.__graph)
        return new_graph

    # Graph access

    @_handle_exceptions
    @_lazy_load
    def vertices(self):
        """Return iterator over the vertices"""
        return self.__graph.Vertices()

    @_handle_exceptions
    @_lazy_load
    def vertex(self, i):
        """Return the i-th vertex from the graph"""
        return self.__graph.Vertex(i)

    @_handle_exceptions
    @_lazy_load
    def edges(self):
        """Return iterator over the edges"""
        return self.__graph.Edges()

    @_handle_exceptions
    @_lazy_load
    def add_vertex(self):
        """Add a new vertex to the graph, and return it"""
        return self.__graph.AddVertex()

    @_handle_exceptions
    @_lazy_load
    def remove_vertex(self, vertex, reindex_edges=True):
        """Remove a vertex from the graph"""
        k = vertex.in_degree() + vertex.out_degree()
        self.__graph.RemoveVertex(vertex)
        if k > 0 and reindex_edges:
            self.__graph.ReIndexEdges()

    @_handle_exceptions
    @_lazy_load
    def add_edge(self, source, target):
        """Add a new edge from 'source' to 'target' to the graph, and return
        it"""
        return self.__graph.AddEdge(source, target)

    @_handle_exceptions
    @_lazy_load
    def remove_edge(self, edge, reindex=True):
        """Remove a edge from the graph"""
        self.__graph.RemoveEdge(edge)
        if reindex:
            self.__graph.ReIndexEdges()

    @_handle_exceptions
    @_lazy_load
    def reindex_edges(self):
        """re-index all edges from 0 to g.num_edges() - 1"""
        self.__graph.ReIndexEdges()

    @_handle_exceptions
    @_lazy_load
    def clear(self):
        """Remove all vertices and edges from the graph (as well as its
        properties)"""
        self.__graph.Clear()

    @_handle_exceptions
    @_lazy_load
    def __get_vertex_properties(self):
        return PropertyDict(self, self.__graph.GetVertexProperties(),
                            lambda g, key: g.get_vertex_property(key),
                            None,
                            lambda g, key: g.remove_vertex_property(key),
                            lambda g, key: g.__graph.GetVertexProperties()[key]\
                            .value_type())
    vertex_properties = property(__get_vertex_properties,
                                 doc="Dictionary of vertex properties")

    @_handle_exceptions
    @_lazy_load
    def __get_edge_properties(self):
        return PropertyDict(self, self.__graph.GetEdgeProperties(),
                            lambda g, key: g.get_edge_property(key),
                            None,
                            lambda g, key: g.remove_edge_property(key),
                            lambda g, key: g.__graph.GetEdgeProperties()[key]\
                            .value_type())
    edge_properties = property(__get_edge_properties,
                                 doc="Dictionary of edge properties")

    @_handle_exceptions
    @_lazy_load
    def __get_graph_properties(self):
        valdict = dict((k, v[self.__graph]) \
                       for k,v in self.__graph.GetGraphProperties().iteritems())
        return PropertyDict(self, valdict,
                            lambda g, key: g.get_graph_property(key),
                            lambda g, key, val: g.set_graph_property(key, val),
                            lambda g, key: g.remove_graph_property(key),
                            lambda g, key: g.__graph.GetGraphProperties()[key].value_type())

    @_handle_exceptions
    @_lazy_load
    def __set_graph_properties(self, val):
        props = self.__graph.GetGraphProperties()
        for k,v in val.iteritems():
            if not props.has_key(k):
                if v.__class__ == str:
                    typename = "string"
                elif v.__class__ == float:
                    typename = "double"
                elif v.__class__ == int:
                    typename = "int"
                elif v.__class__ == bool:
                    typename = "boolean"
                else:
                    raise TypeError("invalid value type for graph " + \
                                    "property '%s'" % k)
                self.add_graph_property(k, typename)
                props = self.__graph.GetGraphProperties()
            props[k][self.__graph] = v

    graph_properties = property(__get_graph_properties, __set_graph_properties,
                                doc="Dictionary of graph properties")

    # Basic options (File i/o and such)
    __groups = ["Basic Options"]

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    def load(self, filename, format="auto"):
        """Load graph from 'filename'. The format is guessed from the file name,
        or can be specified by 'format', which can be either 'xml' or 'dot'."""

        if format == 'auto':
            self.__graph.ReadFromFile(filename)
        else:
            self.__graph.ReadFromFile(filename, format)
        self.lazy_filename = filename
        self.lazy_format = format

    @_handle_exceptions
    def lazy_load(self, filename, format="auto"):
        """Calls load() with the same arguments when the graph is first
        accessed."""

        self.lazy_filename = filename
        self.lazy_format = format

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def save(self, filename, format="auto"):
        """Save graph to file. The format is guessed from the 'file' name, or
        can be specified by 'format', which can be either 'xml' or 'dot'."""

        if format == 'auto':
            self.__graph.WriteToFile(filename)
        else:
            self.__graph.WriteToFile(filename, format)

    # Graph filtering
    __groups.append("Filtering")

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def directed(self):
        """Treat graph as directed (default)."""
        self.__graph.SetDirected(True)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def undirected(self):
        """Treat graph as undirected."""
        self.__graph.SetDirected(False)

    @_handle_exceptions
    @_lazy_load
    def set_directed(self, is_directed):
        """Set the directedness of the graph"""
        self.__graph.SetDirected(is_directed)

    @_handle_exceptions
    @_lazy_load
    def is_directed(self):
        """Get the directedness of the graph"""
        return self.__graph.GetDirected()

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def reversed(self):
        """Reverse the direction of the edges."""
        self.__graph.SetReversed(not self.__graph.GetReversed())

    @_handle_exceptions
    @_lazy_load
    def set_reversed(self, is_reversed):
        """Reverse the direction of the edges, if 'reversed' is True, or
        maintain the original direction otherwise."""
        self.__graph.SetReversed(is_reversed)

    @_handle_exceptions
    @_lazy_load
    def is_reversed(self):
        """Return 'True' if the edges are reversed, and 'False' otherwise."""
        return self.__graph.GetReversed()

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def set_vertex_filter(self, property, inverted=False):
        """Choose vertex boolean property"""
        self.__graph.SetVertexFilterProperty(property, inverted)

    @_handle_exceptions
    @_lazy_load
    def get_vertex_filter(self):
        """Get the vertex filter property and whether or not it is inverted, or
        None if filter is not active"""
        if self.__graph.IsVertexFilterActive():
            return self.__graph.GetVertexFilterProperty()

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def set_edge_filter(self, property, inverted=False):
        """Choose edge boolean property"""
        self.__graph.SetEdgeFilterProperty(property, inverted)

    @_handle_exceptions
    @_lazy_load
    def get_edge_filter(self, property):
        """Get the edge filter property and whether or not it is inverted, or
        None if filter is not active"""
        if self.__Graph.IsEdgeFilterActive():
            return self.__Graph.GetEdgeFilterProperty()

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def exclude_vertex_range(self, property, range):
        """Choose vertex property and range to exclude."""
        (ran, inc, inverted) = _parse_range(range)
        self.__graph.SetVertexFilterRange(property, ran, inc, not inverted)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def keep_vertex_range(self, property, range):
        """Choose vertex property and range to keep (and exclude the rest)."""
        (ran, inc, inverted) = _parse_range(range)
        self.__graph.SetVertexFilterRange(property, ran, inc, inverted)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def reset_vertex_filter(self):
        """Remove edge filter."""
        self.__graph.SetVertexFilterProperty('',False)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def exclude_edge_range(self, property, edge_range):
        """Choose edge property and range to exclude."""
        (ran, inc, inverted) = _parse_range(edge_range)
        self.__graph.SetEdgeFilterRange(property, ran, inc, not inverted)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def keep_edge_range(self, property, range):
        """Choose edge property and range to keep (and exclude the rest)."""
        (ran, inc, inverted) = _parse_range(range)
        self.__graph.SetEdgeFilterRange(property, ran, inc, inverted)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    def reset_edge_filter(self):
        """Remove edge filter"""
        self.__graph.SetEdgeFilterProperty('',False)

    # Graph modification
    __groups.append("Graph Modification")

    @_handle_exceptions
    @_lazy_load
    def add_vertex_property(self, name, type):
        """Add a new (uninitialized) vertex property map of type 'type', and
        return it"""
        self.__graph.AddVertexProperty(name, type)
        return self.vertex_properties[name]

    @_handle_exceptions
    @_lazy_load
    def add_edge_property(self, name, type):
        """Add a new (uninitialized) edge property map of type 'type', and
        return it"""
        self.__graph.AddEdgeProperty(name, type)
        return self.edge_properties[name]

    @_handle_exceptions
    @_lazy_load
    def add_graph_property(self, name, type, val=None):
        """Add a new (uninitialized) graph property map of type 'type'"""
        self.__graph.AddGraphProperty(name, type)
        if val != None:
            self.set_graph_property(name, val)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def set_graph_property(self, property, val):
        """Set the selected graph property"""
        try:
            self.__graph.GetGraphProperties()[property][self.__graph] = val
        except KeyError:
            raise
        except:
            type_name = self.__graph.GetGraphProperties()[property].get_type()
            raise ValueError("wrong value for type '%s': '%s'" % \
                             (type_name, str(val)))

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def get_graph_property(self, property):
        """Get the selected graph property"""
        return self.graph_properties[property]

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def remove_vertex_property(self, property):
        """Remove the selected vertex property"""
        self.__graph.RemoveVertexProperty(property)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def remove_edge_property(self, property):
        """Remove the selected edge property"""
        self.__graph.RemoveEdgeProperty(property)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def remove_graph_property(self, property):
        """Remove the selected graph property"""
        self.__graph.RemoveGraphProperty(property)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def insert_vertex_index(self, property):
        """Insert vertex index as property"""
        self.__graph.InsertVertexIndexProperty(property)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def insert_edge_index(self, property):
        """Insert edge index as property"""
        self.__graph.InsertEdgeIndexProperty(property)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def list_properties(self):
        """List all properties"""
        w = max([len(x) for x in self.graph_properties.keys() + \
                 self.vertex_properties.keys() + \
                 self.edge_properties.keys()]) + 4
        w = w if w > 14 else 14
        for k,v in self.graph_properties.iteritems():
            print "%%-%ds (graph)   (type: %%s, val: %%s)" % w % \
                  (k, self.graph_properties.get_type(k), str(v))
        for k,v in self.vertex_properties.iteritems():
            print "%%-%ds (vertex)  (type: %%s)" % w % (k, v.get_type())
        for k,v in self.edge_properties.iteritems():
            print "%%-%ds (edge)    (type: %%s)" % w % (k, v.get_type())

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def purge_vertices(self):
        """Remove all vertices of the graph which are currently being filtered
        out, and return to the unfiltered state"""
        self.__graph.PurgeVertices()
        self.__graph.SetVertexFilterProperty('')

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def purge_edges(self):
        """Remove all edges of the graph which are currently being filtered out,
        and return to the unfiltered state"""
        self.__graph.PurgeEdges()
        self.__graph.SetEdgeFilterProperty('')

    # Basic graph statistics
    __groups.append("Basic Statistics")

    @_attrs(opt_group=__groups[-1], has_output=True)
    @_handle_exceptions
    @_lazy_load
    def num_vertices(self):
        """Get the number of vertices."""
        return self.__graph.GetNumberOfVertices()

    @_attrs(opt_group=__groups[-1], has_output=True)
    @_handle_exceptions
    @_lazy_load
    def num_edges(self):
        """Get the number of edges"""
        return self.__graph.GetNumberOfEdges()

class PropertyDict(dict):
    """Wrapper for the dict of vertex, graph or edge properties, which sets the
    value on the property map when changed in the dict."""
    def __init__(self, g, old, get_func, set_func, del_func, type_func):
        dict.__init__(self)
        dict.update(self, old)
        self.g = g
        self.get_func = get_func
        self.set_func = set_func
        self.del_func = del_func
    def __setitem__(self, key, val):
        if self.set_func != None:
            self.set_func(self.g, key, val)
        else:
            raise KeyError("Property dict cannot be set")
        dict.__setitem__(self, key, val)
    def __delitem__(self, key):
        self.del_func(self.g, key)
        dict.__delitem__(self, key)
    def get_type(self, key):
        """Return the type of the internal property"""
        return self.type_func(self.g, key)
