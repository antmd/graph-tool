#! /usr/bin/env python
# graph_tool -- a general graph manipulation python module
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
import libgraph_tool_core as libcore
sys.setdlopenflags(_orig_dlopen_flags) # reset it to normal case to avoid
                                       # unnecessary symbol collision
__version__ = libcore.mod_info().version

import io # sets up libcore io routines

import os, os.path, re, struct, fcntl, termios, gzip, bz2, string,\
       textwrap, time, signal, traceback, shutil, time, math, inspect, \
       functools, types
from StringIO import StringIO
from decorators import _wraps, _require, _attrs, _handle_exceptions, _limit_args

################################################################################
# Utility functions
################################################################################

_require("name", str)
def _degree(name):
    """Retrieve the degree type from string"""
    deg = name
    if name == "in-degree" or name == "in":
        deg = libcore.Degree.In
    if name == "out-degree" or name == "out":
        deg = libcore.Degree.Out
    if name == "total-degree" or name == "total":
        deg = libcore.Degree.Total
    return deg

def _prop(t, g, prop):
    """Returns either a property map, or an internal vertex property map with a
    given name"""
    if type(prop) == str:
        try:
            pmap = g.properties[(t,prop)]
        except KeyError:
            raise GraphError("no internal %s property named: %d" %\
                             ("vertex" if t == "v" else \
                              ("edge" if t == "e" else "graph"),prop))
    else:
        pmap = prop
    return pmap._PropertyMap__map

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
# Property Maps
################################################################################

class PropertyMap(object):
    """Property Map class"""
    def __init__(self, pmap, g, key_type, key_trans = None):
        self.__map = pmap
        self.__g = g
        self.__key_type = key_type
        self.__key_trans = key_trans if key_trans != None else lambda k: k
        self.__register_map()

    def __register_map(self):
        self.__g._Graph__known_properties.append((self.key_type(), self.__map))
    def __unregister_map(self):
        i = self.__g._Graph__known_properties.index((self.key_type(), self.__map))
        del self.__g._Graph__known_properties[i]

    def __del__(self):
        self.__unregister_map()

    def __getitem__(self, k):
        return self.__map[self.__key_trans(k)]

    def __setitem__(self, k, v):
        self.__map[self.__key_trans(k)] = v

    def get_graph(self):
        """Get the graph to which the map refers"""
        return self.__g

    def key_type(self):
        """The key type of the map. Either 'g', 'v' or 'e'"""
        return self.__key_type

    def value_type(self):
        """The value type of the map"""
        return self.__map.value_type()

class PropertyDict(dict):
    """Wrapper for the dict of vertex, graph or edge properties, which sets the
    value on the property map when changed in the dict."""
    def __init__(self, g, old, get_func, set_func, del_func):
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

################################################################################
# Graph class
# The main graph interface
################################################################################

from libgraph_tool_core import Vertex, Edge, GraphError,\
     Vector_bool, Vector_int32_t,  Vector_int64_t, Vector_double,\
     Vector_long_double, Vector_string, new_vertex_property, new_edge_property,\
     new_graph_property

class Graph(object):
    """The main graph type"""

    def __init__(self, g = None):
        self.__properties = {}
        self.__known_properties = []
        self.__filter_state = {"reversed": False,
                               "edge_filter": (None,False),
                               "vertex_filter": (None,False),
                               "directed": True,
                               "reversed": False}
        self.__stashed_filter_state = []

        if g == None:
            self.__graph = libcore.GraphInterface()
        else:
            self.__graph = libcore.GraphInterface(g.__graph)
            for k,v in g.__properties.iteritems():
                new_p = self.new_property(v.key_type(), v.value_type())
                self.copy_property(v, new_p, g)
                self.properties[(v.key_type(),k)] = new_p
            self.__stashed_filter_state = [g.__filter_state]
            self.pop_filter()

    @_handle_exceptions
    def copy(self):
        """Returns a deep copy of self"""
        new_graph = Graph(self)
        return new_graph

    # Graph access
    # ============

    @_handle_exceptions
    def vertices(self):
        """Return iterator over the vertices"""
        return self.__graph.Vertices()

    @_handle_exceptions
    def vertex(self, i):
        """Return the i-th vertex from the graph"""
        return self.__graph.Vertex(int(i))

    @_handle_exceptions
    def edges(self):
        """Return iterator over the edges"""
        return self.__graph.Edges()

    @_handle_exceptions
    def add_vertex(self):
        """Add a new vertex to the graph, and return it"""
        return self.__graph.AddVertex()

    @_handle_exceptions
    def remove_vertex(self, vertex, reindex_edges=True):
        """Remove a vertex from the graph"""
        k = vertex.in_degree() + vertex.out_degree()
        index = self.vertex_index[vertex]
        for pmap in self.__known_properties:
            if pmap[0] == "v":
                self.__graph.ShiftVertexProperty(pmap[1].get_map(), index)
        self.__graph.RemoveVertex(vertex)

    @_handle_exceptions
    def add_edge(self, source, target):
        """Add a new edge from 'source' to 'target' to the graph, and return
        it"""
        return self.__graph.AddEdge(source, target)

    @_handle_exceptions
    def remove_edge(self, edge):
        """Remove an edge from the graph"""
        self.__graph.RemoveEdge(edge)

    @_handle_exceptions
    def clear(self):
        """Remove all vertices and edges from the graph"""
        self.__graph.Clear()

    @_handle_exceptions
    def clear_edges(self):
        """Remove all edges from the graph"""
        self.__graph.ClearEdges()

    # Internal property maps
    # ======================

    # all properties
    @_handle_exceptions
    def __get_properties(self):
        return PropertyDict(self, self.__properties,
                            lambda g,k: g.__properties[k],
                            lambda g,k,v: g.__set_property(k[0],k[1],v),
                            lambda g,k: g.__del_property(k[0],k[1]))

    @_handle_exceptions
    @_limit_args({"t":["v", "e", "g"]})
    @_require("k", str)
    @_require("v", PropertyMap)
    def __set_property(self, t, k, v):
        if t != v.key_type():
            raise ValueError("wrong key type for property map")
        self.__properties[(t,k)] = v

    @_handle_exceptions
    @_limit_args({"t":["v", "e", "g"]})
    @_require("k", str)
    def __del_property(self, t, k):
        del self.__properties[(t,k)]

    properties = property(__get_properties,
                          doc="Dictionary of internal properties")

    def __get_specific_properties(self, t):
        props = dict([(k[1],v) for k,v in self.__properties.iteritems() \
                      if k[0] == t ])
        return props

    # vertex properties
    @_handle_exceptions
    def __get_vertex_properties(self):
        return PropertyDict(self, self.__get_specific_properties("v"),
                            lambda g,k: g.__properties[("v",k)],
                            lambda g,k,v: g.__set_property("v",k,v),
                            lambda g,k: g.__del_property("v",k))
    vertex_properties = property(__get_vertex_properties,
                                 doc="Dictionary of vertex properties")
    # edge properties
    @_handle_exceptions
    def __get_edge_properties(self):
        return PropertyDict(self, self.__get_specific_properties("e"),
                            lambda g,k: g.__properties[("e",k)],
                            lambda g,k,v: g.__set_property("e",k,v),
                            lambda g,k: g.__del_property("e",k))
    edge_properties = property(__get_edge_properties,
                                 doc="Dictionary of edge properties")

    # graph properties
    @_handle_exceptions
    def __get_graph_properties(self):
        return PropertyDict(self, self.__get_specific_properties("g"),
                            lambda g,k: g.__properties[("g",k)],
                            lambda g,k,v: g.__set_property("g",k,v),
                            lambda g,k: g.__del_property("g",k))
    graph_properties = property(__get_graph_properties,
                                 doc="Dictionary of graph properties")

    @_handle_exceptions
    def list_properties(self):
        """List all internal properties"""
        if len(self.__properties) == 0:
            return
        w = max([len(x[0]) for x in self.__properties.keys()]) + 4
        w = w if w > 14 else 14

        for k,v in self.__properties.iteritems():
            if k[0] == "g":
                print "%%-%ds (graph)   (type: %%s, val: %%s)" % w % \
                      (k[1], v.value_type(), str(v[self]))
        for k,v in self.__properties.iteritems():
            if k[0] == "v":
                print "%%-%ds (vertex)  (type: %%s)" % w % (k[1],
                                                            v.value_type())
        for k,v in self.__properties.iteritems():
            if k[0] == "e":
                print "%%-%ds (edge)    (type: %%s)" % w % (k[1],
                                                            v.value_type())

    # index properties

    def _get_vertex_index(self):
        return PropertyMap(libcore.get_vertex_index(self.__graph), self, "v")
    vertex_index = property(_get_vertex_index, doc="Vertex index map")

    def _get_edge_index(self):
        return PropertyMap(libcore.get_edge_index(self.__graph), self, "v")
    edge_index = property(_get_edge_index, doc="Edge index map")

    # Property map creation

    @_handle_exceptions
    def new_property(self, key_type, type):
        """Create a new (uninitialized) vertex property map of key type
        'key_type', value type 'type', and return it"""
        if key_type == "v" or key_type == "vertex":
            return self.new_vertex_property(type)
        if key_type == "e" or key_type == "edge":
            return self.new_edge_property(type)
        if key_type == "g" or key_type == "graph":
            return self.new_graph_property(type)
        raise GraphError("unknown key type: " + key_type)

    @_handle_exceptions
    def new_vertex_property(self, type):
        """Create a new (uninitialized) vertex property map of type 'type', and
        return it"""
        return PropertyMap(new_vertex_property(type,
                                               self.__graph.GetVertexIndex()),
                           self, "v")

    @_handle_exceptions
    def new_edge_property(self, type):
        """Create a new (uninitialized) edge property map of type 'type', and
        return it"""
        return PropertyMap(new_edge_property(type, self.__graph.GetEdgeIndex()),
                           self, "e")

    @_handle_exceptions
    def new_graph_property(self, type):
        """Create a new (uninitialized) graph property map of type 'type', and
        return it"""
        return PropertyMap(new_graph_property(type,
                                              self.__graph.GetGraphIndex()),
                           self, "g", lambda k: k.__graph)

    # property map copying
    @_handle_exceptions
    @_require("src", PropertyMap)
    @_require("tgt", PropertyMap)
    def copy_property(self, src, tgt, g=None):
        """Copy contents of src property to tgt property. Parameter g specifices
        the (identical) source graph to copy properties from (defaults to self)
        """
        if src.key_type() != tgt.key_type():
            raise GraphError("source and target properties must have the same" +
                             " key type")
        if g == None:
            g = self
        if g != self:
            g.stash_filter()
        self.stash_filter()
        if src.key_type() == "v":
            self.__graph.CopyVertexProperty(g.__graph, src._PropertyMap__map,
                                            tgt._PropertyMap__map)
        elif src.key_type() == "e":
            self.__graph.CopyEdgeProperty(g.__graph, src._PropertyMap__map,
                                          tgt._PropertyMap__map)
        else:
            tgt[self] = src[g]
        self.pop_filter()
        if g != self:
            g.pop_filter()

    # I/O operations
    # ==============

    @_handle_exceptions
    def load(self, filename, format="auto"):
        """Load graph from 'filename' (which can also be a file-like
        object). The format is guessed from the file name, or can be specified
        by 'format', which can be either 'xml' or 'dot'."""
        if format == 'auto' and isinstance(filename, str):
            if filename.endswith(".xml") or filename.endswith(".xml.gz") or \
                   filename.endswith(".xml.bz2"):
                format = "xml"
            elif filename.endswith(".dot") or filename.endswith(".dot.gz") or \
                     filename.endswith(".dot.bz2"):
                format = "dot"
            else:
                libcore.raise_error\
                    ("cannot determine file format of: " + filename )
        elif format == "auto":
            format = "xml"
        if isinstance(filename, str):
            props = self.__graph.ReadFromFile(filename, None, format)
        else:
            props = self.__graph.ReadFromFile("", filename, format)
        for name, prop in props[0].iteritems():
            self.vertex_properties[name] = PropertyMap(prop, self, "v")
        for name, prop in props[1].iteritems():
            self.edge_properties[name] = PropertyMap(prop, self, "e")
        for name, prop in props[2].iteritems():
            self.graph_properties[name] = PropertyMap(prop, self, "g",
                                                      lambda k: k.__graph)

    @_handle_exceptions
    def save(self, filename, format="auto"):
        """Save graph to file. The format is guessed from the 'file' name, or
        can be specified by 'format', which can be either 'xml' or 'dot'."""
        if format == 'auto' and isinstance(filename, str):
            if filename.endswith(".xml") or filename.endswith(".xml.gz") or \
                   filename.endswith(".xml.bz2"):
                format = "xml"
            elif filename.endswith(".dot") or filename.endswith(".dot.gz") or \
                     filename.endswith(".dot.bz2"):
                format = "dot"
            else:
                libcore.raise_error\
                    ("cannot determine file format of: " + filename )
        elif format == "auto":
            format = "xml"
        props = [(name[1], prop._PropertyMap__map) for name,prop in \
                 self.__properties.iteritems()]
        if isinstance(filename, str):
            self.__graph.WriteToFile(filename, None, format, props)
        else:
            self.__graph.WriteToFile("", filename, format, props)

    # Directedness
    # ============

    @_handle_exceptions
    def directed(self):
        """Treat graph as directed (default)."""
        self.__graph.SetDirected(True)

    @_handle_exceptions
    def undirected(self):
        """Treat graph as undirected."""
        self.__graph.SetDirected(False)

    @_handle_exceptions
    def set_directed(self, is_directed):
        """Set the directedness of the graph"""
        self.__graph.SetDirected(is_directed)

    @_handle_exceptions
    def is_directed(self):
        """Get the directedness of the graph"""
        return self.__graph.GetDirected()

    # Reversedness
    # ============

    @_handle_exceptions
    def reversed(self):
        """Reverse the direction of the edges."""
        self.__graph.SetReversed(not self.__graph.GetReversed())

    @_handle_exceptions
    def set_reversed(self, is_reversed):
        """Reverse the direction of the edges, if 'reversed' is True, or
        maintain the original direction otherwise."""
        self.__graph.SetReversed(is_reversed)

    @_handle_exceptions
    def is_reversed(self):
        """Return 'True' if the edges are reversed, and 'False' otherwise."""
        return self.__graph.GetReversed()

    # Filtering
    # =========

    @_handle_exceptions
    def set_vertex_filter(self, property, inverted=False):
        """Choose vertex boolean filter property"""
        self.__graph.SetVertexFilterProperty(property, inverted)
        self.__filter_state["vertex_filter"] = (property, inverted)

    @_handle_exceptions
    def get_vertex_filter(self):
        """Get the vertex filter property and whether or not it is inverted"""
        return self.__filter_state["vertex_filter"]

    @_handle_exceptions
    def set_edge_filter(self, property, inverted=False):
        """Choose edge boolean property"""
        self.__graph.SetEdgeFilterProperty(property, inverted)
        self.__filter_state["edge_filter"] = (property, inverted)

    @_handle_exceptions
    def get_edge_filter(self):
        """Get the edge filter property and whether or not it is inverted"""
        return self.__filter_state["edge_filter"]

    @_handle_exceptions
    def purge_vertices(self):
        """Remove all vertices of the graph which are currently being filtered
        out, and return to the unfiltered state"""
        self.__graph.PurgeVertices()
        self.__graph.SetVertexFilterProperty('')

    @_handle_exceptions
    def purge_edges(self):
        """Remove all edges of the graph which are currently being filtered out,
        and return to the unfiltered state"""
        self.__graph.PurgeEdges()
        self.__graph.SetEdgeFilterProperty('')

    @_handle_exceptions
    def stash_filter(self):
        """Stash current filter state and recover unfiltered graph"""
        self.__stashed_filter_state.append(self.__filter_state)
        self.set_vertex_filter("")
        self.set_edge_filter("")
        self.directed()
        self.set_reversed(False)

    @_handle_exceptions
    def pop_filter(self):
        """Pop last stashed filter state"""
        state = self.__stashed_filter_state.pop()
        self.set_vertex_filter(state["vertex_filter"][0],
                               state["vertex_filter"][1])
        self.set_edge_filter(state["edge_filter"][0],
                             state["edge_filter"][1])
        self.set_directed(state["directed"])
        self.set_reversed(state["reversed"])


    # Basic graph statistics
    # ======================

    @_handle_exceptions
    def num_vertices(self):
        """Get the number of vertices."""
        return self.__graph.GetNumberOfVertices()

    @_handle_exceptions
    def num_edges(self):
        """Get the number of edges"""
        return self.__graph.GetNumberOfEdges()

    def underlying_graph(self):
        """Retrieve a GraphInterface from a Graph (for internal use only)"""
        return self.__graph


    # Pickling support
    # ================

    def __getstate__(self):
        state = dict()
        state["vertex_filt"] = self.get_vertex_filter()
        self.reset_vertex_filter()
        state["edge_filt"] = self.get_edge_filter()
        self.reset_edge_filter()
        sio = StringIO()
        self.save(sio, "xml")
        state["blob"] = sio.getvalue()
        return state
    def __setstate__(self, state):
        self.__init__()
        blob = state["blob"]
        if blob != "":
            sio = StringIO(blob)
            self.load(sio, "xml")
        if state["vertex_filt"] != None:
            self.set_vertex_filter(state["vertex_filt"][0],
                                   state["vertex_filt"][1])
        if state["edge_filt"] != None:
            self.set_edge_filter(state["edge_filt"][0],
                                 state["edge_filt"][1])

def value_types():
    """Return a list of possible properties value types"""
    return libcore.get_property_types()
