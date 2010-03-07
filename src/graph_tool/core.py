#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# graph_tool -- a general graph manipulation python module
#
# Copyright (C) 2007-2010 Tiago de Paula Peixoto <tiago@forked.de>
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

from dl_import import *
dl_import("import libgraph_tool_core as libcore")
import libgraph_tool_core as libcore   # for pylint
__version__ = libcore.mod_info().version

import io  # sets up libcore io routines

import sys
import os, os.path, re, struct, fcntl, termios, gzip, bz2, string,\
       textwrap, time, signal, traceback, shutil, time, math, inspect, \
       functools, types, weakref, copy
from StringIO import StringIO
from decorators import _wraps, _require, _attrs, _limit_args

################################################################################
# Utility functions
################################################################################


def _prop(t, g, prop):
    """Return either a property map, or an internal property map with a given
    name."""
    if type(prop) == str:
        try:
            pmap = g.properties[(t, prop)]
        except KeyError:
            raise KeyError("no internal %s property named: %s" %\
                           ("vertex" if t == "v" else \
                            ("edge" if t == "e" else "graph"), prop))
    else:
        pmap = prop
    if pmap == None:
        return libcore.any()
    else:
        if t != prop.key_type():
            names = {'e': 'edge', 'v': 'vertex', 'g': 'graph'}
            raise ValueError("Expected '%s' property map, got '%s'" %
                             (names[t], names[prop.key_type()]))
        return pmap._PropertyMap__map.get_map()


def _degree(g, name):
    """Retrieve the degree type from string, or returns the corresponding
    property map."""
    deg = name
    if name == "in-degree" or name == "in":
        deg = libcore.Degree.In
    elif name == "out-degree" or name == "out":
        deg = libcore.Degree.Out
    elif name == "total-degree" or name == "total":
        deg = libcore.Degree.Total
    else:
        deg = _prop("v", g, deg)
    return deg


def _type_alias(type_name):
    alias = {"int8_t": "bool",
             "boolean": "bool",
             "int": "int32_t",
             "long": "int32_t",
             "long long": "int64_t",
             "object": "python::object"}
    if type_name in value_types():
        return type_name
    if type_name in alias:
        return alias[type_name]
    ma = re.compile(r"vector<(.*)>").match(type_name)
    if ma:
        t = ma.group(1)
        if t in alias:
            return "vector<%s>" % alias[t]
    raise ValueError("invalid property value type: " + type_name)


def show_config():
    """Show ``graph_tool`` build configuration."""
    info = libcore.mod_info()
    print "version:", info.version
    print "gcc version:", info.gcc_version
    print "compilation flags:", info.cxxflags
    print "install prefix:", info.install_prefix
    print "python dir:", info.python_dir
    print "graph filtering:", libcore.graph_filtering_enabled()
    print "openmp:", libcore.openmp_enabled()
    print "uname:", " ".join(os.uname())

################################################################################
# Property Maps
################################################################################


class PropertyMap(object):
    """Property map class."""
    def __init__(self, pmap, g, key_type, key_trans=None):
        self.__map = pmap
        self.__g = weakref.ref(g)
        self.__key_type = key_type
        self.__key_trans = key_trans if key_trans != None else lambda k: k
        self.__register_map()

    def __register_map(self):
        if self.__g() == None:
            return
        self.__g()._Graph__known_properties.append((self.key_type(),
                                                    weakref.ref(self.__map)))

    def __unregister_map(self):
        if self.__g() == None:
            return
        i = self.__g()._Graph__known_properties.index((self.key_type(),
                                                       weakref.ref(self.__map)))
        del self.__g()._Graph__known_properties[i]

    def __del__(self):
        self.__unregister_map()

    def __getitem__(self, k):
        return self.__map[self.__key_trans(k)]

    def __setitem__(self, k, v):
        self.__map[self.__key_trans(k)] = v

    def __repr__(self):
        # provide some more useful information
        if self.key_type() == "e":
            k = "Edge"
        elif self.key_type() == "v":
            k = "Vertex"
        else:
            k = "Graph"
        g = self.__g()
        if g == None:
            g = "a non-existent graph"
        else:
            g = "Graph 0x%x" % id(g)
        return ("<PropertyMap object with key type '%s' and value type '%s',"
                + " for %s, at 0x%x>") % (k, self.value_type(), g, id(self))

    def get_graph(self):
        """Get the graph class to which the map refers."""
        return self.__g()

    def key_type(self):
        """Return the key type of the map. Either 'g', 'v' or 'e'."""
        return self.__key_type

    def value_type(self):
        """Return the value type of the map."""
        return self.__map.value_type()

    def get_array(self):
        """Get an array with property values.

        .. WARNING::

           The returned array does not own the data, which belongs to the
           property map. Therefore, the returned array cannot have a longer
           lifespan than the property map itself! Furthermore, if the graph
           changes, it may leave the pointer to the data in the array dangling!
           Do *not* store the array if the graph is to be modified, or the
           original property map deleted; *store a copy instead*!
        """
        self.__g().stash_filter(edge=True, vertex=True)
        if self.__key_type == 'v':
            n = self.__g().num_vertices()
        elif self.__key_type == 'e':
            n = self.__g().num_edges()
        else:
            n = 1
        self.__g().pop_filter(edge=True, vertex=True)
        return self.__map.get_array(n)

    def __get_array(self):
        if self.get_array() != None:
            return self.get_array()[:]
        else:
            return None

    def __set_array(self, v):
        self.get_array()[:] = v

    a = property(__get_array, __set_array,
                 doc=r"""Shortcut to the :meth:`~PropertyMap.get_array` method
                 as a property. A view to the array is returned, instead of the
                 array, for convenience.""")

    def is_writable(self):
        """Return True if the property is writable."""
        return self.__map.is_writable()


def _check_prop_writable(prop, name=None):
    if not prop.is_writable():
        raise ValueError("property map%s is not writable." %\
                         ((" '%s'" % name) if name != None else ""))


def _check_prop_scalar(prop, name=None, floating=False):
    scalars = ["bool", "int32_t", "int64_t", "unsigned long",
               "double", "long double"]
    if floating:
        scalars = ["double", "long double"]

    if prop.value_type() not in scalars:
        raise ValueError("property map%s is not of scalar%s type." %\
                         (((" '%s'" % name) if name != None else ""),
                          (" floating" if floating else "")))


def _check_prop_vector(prop, name=None, scalar=True, floating=False):
    scalars = ["bool", "int32_t", "int64_t", "unsigned long",
               "double", "long double"]
    if not scalar:
        scalars += ["string"]
    if floating:
        scalars = ["double", "long double"]
    vals = ["vector<%s>" % v for v in scalars]
    if prop.value_type() not in vals:
        raise ValueError("property map%s is not of vector%s type." %\
                         (((" '%s'" % name) if name != None else ""),
                          (" floating" if floating else "")))


def group_vector_property(g, props, value_type=None, vprop=None, pos=None):
    """Group list of properties ``props`` into a vector property map of the same
    type.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to which the property maps belong.
    props : list of :class:`~graph_tool.PropertyMap`
        Properties to be grouped.
    value_type : string (optional, default: None)
        If supplied, defines the value type of the grouped property.
    vprop : :class:`~graph_tool.PropertyMap` (optional, default: None)
        If supplied, the properties are grouped into this property map.
    pos : list of ints (optional, default: None)
        If supplied, should contain a list of indexes where each corresponding
        element of ``props`` should be inserted.

    Returns
    -------
    vprop : :class:`~graph_tool.PropertyMap`
       A vector property map with the grouped values of each property map in
       ``props``.
    """
    vtypes = set()
    keys = set()
    for i, p in enumerate(props):
        if "vector" in p.value_type():
            raise ValueError("property map 'props[%d]' is a vector property." %
                             i)
        vtypes.add(p.value_type())
        keys.add(p.key_type())
    if len(keys) > 1:
        raise ValueError("'props' must be of the same key type.")
    k = keys.pop()

    if vprop == None:
        if value_type == None and len(vtypes) == 1:
            value_type = vtypes.pop()

        if value_type != None:
            value_type = "vector<%s>" % value_type
            if k == 'v':
                vprop = g.new_vertex_property(value_type)
            elif k == 'e':
                vprop = g.new_edge_property(value_type)
            else:
                vprop = g.new_graph_property(value_type)
        else:
            ValueError("Can't automatically determine property map value" +
                       " type. Please provide the 'value_type' parameter.")
    _check_prop_vector(vprop, name="vprop", scalar=False)

    for i, p in enumerate(props):
        if k != "g":
            g.stash_filter(directed=True, reversed=True)
            g.set_directed(True)
            g.set_reversed(False)
            libcore.group_vector_property(g._Graph__graph,
                                          _prop(k, g, vprop),
                                          _prop(k, g, p),
                                          i if pos == None else pos[i],
                                          k == 'e')
            g.pop_filter(directed=True, reversed=True)
        else:
            vprop[g][i if pos == None else pos[i]] = p[g]
    return vprop


def ungroup_vector_property(g, vprop, pos, props=None):
    """Ungroup vector property map ``vprop`` into a list of individual property
    maps.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to which the property map belong.
    vprop : :class:`~graph_tool.PropertyMap`
        Vector property map to be ungrouped.
    pos : list of ints (optional, default: None)
        If supplied, should contain a list of indexes where each corresponding
        element of ``vprop`` should be inserted into the ``props`` list.
    props : list of :class:`~graph_tool.PropertyMap`  (optional, default: None)
        If supplied, should contain a list of property maps to which ``vprop``
        should be ungroupped.

    Returns
    -------
    props : list of :class:`~graph_tool.PropertyMap`
       A list of property maps with the ungrouped values of ``vprop``.
    """

    _check_prop_vector(vprop, name="vprop", scalar=False)
    k = vprop.key_type()
    value_type = vprop.value_type().split("<")[1].split(">")[0]
    if props == None:
        if k == 'v':
            props = [g.new_vertex_property(value_type) for i in pos]
        elif k == 'e':
            props = [g.new_edge_property(value_type) for i in pos]
        else:
            props = [g.new_graph_property(value_type) for i in pos]

    for i, p in enumerate(pos):
        if props[i].key_type() != k:
            raise ValueError("'props' must be of the same key type as 'vprop'.")

        if k != 'g':
            g.stash_filter(directed=True, reversed=True)
            g.set_directed(True)
            g.set_reversed(False)
            libcore.ungroup_vector_property(g._Graph__graph,
                                            _prop(k, g, vprop),
                                            _prop(k, g, props[i]),
                                            p, k == 'e')
            g.pop_filter(directed=True, reversed=True)
        else:
            if len(vprop[g]) <= pos[i]:
                vprop[g].resize(pos[i] + 1)
            props[i][g] = vprop[g][pos[i]]
    return props


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

from libgraph_tool_core import Vertex, Edge, Vector_bool, Vector_int32_t, \
     Vector_int64_t, Vector_double, Vector_long_double, Vector_string, \
     new_vertex_property, new_edge_property, new_graph_property


class Graph(object):
    """This class encapsulates either a directed multigraph (default or if
    ``directed=True``) or an undirected multigraph (if ``directed=False``), with
    optional internal edge, vertex or graph properties.

    It is implemented as an adjacency list, where both vertex and edge lists are
    C++ STL vectors.
    """

    def __init__(self, g=None, directed=True):
        """Construct a graph. If ``g`` is specified, the graph (and its internal
        properties) will be copied. The ``directed`` parameter specifies whether
        the graph should be directed or undirected."""
        self.__properties = {}
        self.__known_properties = []
        self.__filter_state = {"reversed": False,
                               "edge_filter": (None, False),
                               "vertex_filter": (None, False),
                               "directed": True}
        self.__stashed_filter_state = []

        if g == None:
            self.__graph = libcore.GraphInterface()
            self.set_directed(directed)
        else:
            self.__graph = libcore.GraphInterface(g.__graph)
            for k, v in g.__properties.iteritems():
                new_p = self.new_property(v.key_type(), v.value_type())
                self.copy_property(v, new_p, g)
                self.properties[k] = new_p

            self.__stashed_filter_state = [self.get_filter_state()]

            v_filt, v_rev = g.__filter_state["vertex_filter"]
            if v_filt != None:
                if v_filt not in g.vertex_properties.values():
                    new_filt = self.new_vertex_property("bool")
                    self.copy_property(v_filt, new_filt)

                else:
                    for k, v in g.vertex_properties.iteritems():
                        if v == v_filt:
                            new_filt = self.vertex_properties[k]
                self.__stashed_filter_state[0]["vertex_filter"] = (new_filt,
                                                                   v_rev)
            e_filt, e_rev = g.__filter_state["edge_filter"]
            if e_filt != None:
                if e_filt not in g.edge_properties.values():
                    new_filt = self.new_edge_property("bool")
                    self.copy_property(e_filt, new_filt)

                else:
                    for k, v in g.edge_properties.iteritems():
                        if v == e_filt:
                            new_filt = self.edge_properties[k]
                self.__stashed_filter_state[0]["edge_filter"] = (new_filt,
                                                                 e_rev)
            self.pop_filter()
        # internal index maps
        self.__vertex_index = \
                 PropertyMap(libcore.get_vertex_index(self.__graph), self, "v")
        self.__edge_index = \
                 PropertyMap(libcore.get_edge_index(self.__graph), self, "e")

    def copy(self):
        """Return a deep copy of self. All internal property maps are also
        copied."""
        return Graph(self)

    def __repr__(self):
        # provide some more useful information
        d = "directed" if self.is_directed() else "undirected"
        fr = ", reversed" if self.is_reversed() and self.is_directed() else ""
        f = ""
        if self.get_edge_filter()[0] != None:
            f += ", edges filtered by %s" % (str(self.get_edge_filter()))
        if self.get_vertex_filter()[0] != None:
            f += ", vertices filtered by %s" % (str(self.get_vertex_filter()))
        n = self.num_vertices()
        e = self.num_edges()
        return "<Graph object, %s%s, with %d %s and %d edge%s%s at 0x%x>"\
               % (d, fr, n, "vertex" if n == 1 else "vertices", e,
                  "" if e == 1 else "s", f, id(self))

    # Graph access
    # ============

    def vertices(self):
        """Return an iterator over the vertices.

        Examples
        --------
        >>> g = gt.Graph()
        >>> vlist = g.add_vertex(5)
        >>> vlist2 = []
        >>> for v in g.vertices():
        ...     vlist2.append(v)
        ...
        >>> assert(vlist == vlist2)

        """
        return libcore.get_vertices(weakref.ref(self.__graph))

    def vertex(self, i, use_index=False):
        """Return the i-th vertex from the graph. If use_index=True, the vertex
        with index i is returned (which can differ from the i-th vertex in case
        of filtered graphs)."""
        if use_index:
            self.stash_filter(vertex=True)
        try:
            v = libcore.get_vertex(weakref.ref(self.__graph), int(i))
        finally:
            if use_index:
                self.pop_filter(vertex=True)
        return v

    def edges(self):
        """Return an iterator over the edges."""
        return libcore.get_edges(weakref.ref(self.__graph))

    def add_vertex(self, n=1):
        """Add a vertex to the graph, and return it. If ``n > 1``, ``n``
        vertices are inserted and a list is returned."""
        vlist = [libcore.add_vertex(weakref.ref(self.__graph)) \
                 for i in xrange(0, n)]
        if n == 1:
            return vlist[0]
        return vlist

    def remove_vertex(self, vertex):
        """Remove a vertex from the graph."""
        index = self.vertex_index[vertex]
        for pmap in self.__known_properties:
            if pmap[0] == "v" and pmap[1]() != None and \
                   pmap[1]() != self.__vertex_index._PropertyMap__map:
                self.__graph.ShiftVertexProperty(pmap[1]().get_map(), index)
        self.clear_vertex(vertex)
        libcore.remove_vertex(self.__graph, vertex)

    def remove_vertex_if(self, predicate):
        """Remove all the vertices from the graph for which ``predicate(v)``
        evaluates to ``True``. """
        N = self.num_vertices()
        for i in xrange(0, N):
            v = self.vertex(N - i - 1)
            if predicate(v):
                self.remove_vertex(v)

    def clear_vertex(self, vertex):
        """Remove all in and out-edges from the given vertex."""
        del_es = set()
        for e in vertex.all_edges():
            del_es.add(e)
        for e in del_es:
            self.remove_edge(e)

    def add_edge(self, source, target):
        """Add a new edge from ``source`` to ``target`` to the graph, and return
        it."""
        return libcore.add_edge(weakref.ref(self.__graph), source, target)

    def remove_edge(self, edge):
        """Remove an edge from the graph."""
        return libcore.remove_edge(self.__graph, edge)

    def remove_edge_if(self, predicate):
        """Remove all edges from the graph, for which ``predicate(e)`` evaluates
        to ``True``."""
        for v in self.vertices():
            del_es = []
            for e in v.out_edges():
                if predicate(e):
                    del_es.append(e)
            for e in del_es:
                self.remove_edge(e)

    def clear(self):
        """Remove all vertices and edges from the graph."""
        self.__graph.Clear()

    def clear_edges(self):
        """Remove all edges from the graph."""
        self.__graph.ClearEdges()

    # Internal property maps
    # ======================

    # all properties
    def __get_properties(self):
        return PropertyDict(self, self.__properties,
                            lambda g, k: g.__properties[k],
                            lambda g, k, v: g.__set_property(k[0], k[1], v),
                            lambda g, k: g.__del_property(k[0], k[1]))

    @_limit_args({"t": ["v", "e", "g"]})
    @_require("k", str)
    @_require("v", PropertyMap)
    def __set_property(self, t, k, v):
        if t != v.key_type():
            raise ValueError("wrong key type for property map")
        self.__properties[(t, k)] = v

    @_limit_args({"t": ["v", "e", "g"]})
    @_require("k", str)
    def __del_property(self, t, k):
        del self.__properties[(t, k)]

    properties = property(__get_properties,
                          doc=
    """Dictionary of internal properties. Keys must always be a tuple, where the
    first element if a string from the set {'v', 'e', 'g'}, representing a
    vertex, edge or graph property, and the second element is the name of the
    property map.

    Examples
    --------
    >>> g = gt.Graph()
    >>> g.properties[("e", "foo")] = g.new_edge_property("vector<double>")
    >>> del g.properties[("e", "foo")]
    """)

    def __get_specific_properties(self, t):
        props = dict([(k[1], v) for k,v in self.__properties.iteritems() \
                      if k[0] == t ])
        return props

    # vertex properties
    def __get_vertex_properties(self):
        return PropertyDict(self, self.__get_specific_properties("v"),
                            lambda g, k: g.__properties[("v", k)],
                            lambda g, k, v: g.__set_property("v", k, v),
                            lambda g, k: g.__del_property("v", k))
    vertex_properties = property(__get_vertex_properties,
                                 doc="Dictionary of vertex properties")

    # edge properties
    def __get_edge_properties(self):
        return PropertyDict(self, self.__get_specific_properties("e"),
                            lambda g, k: g.__properties[("e", k)],
                            lambda g, k, v: g.__set_property("e", k, v),
                            lambda g, k: g.__del_property("e", k))
    edge_properties = property(__get_edge_properties,
                                 doc="Dictionary of edge properties")

    # graph properties
    def __get_graph_properties(self):
        return PropertyDict(self, self.__get_specific_properties("g"),
                            lambda g, k: g.__properties[("g", k)],
                            lambda g, k, v: g.__set_property("g", k, v),
                            lambda g, k: g.__del_property("g", k))
    graph_properties = property(__get_graph_properties,
                                 doc="Dictionary of graph properties")

    def list_properties(self):
        """List all internal properties.

        Examples
        --------
        >>> g = gt.Graph()
        >>> g.properties[("e", "foo")] = g.new_edge_property("vector<double>")
        >>> g.vertex_properties["foo"] = g.new_vertex_property("double")
        >>> g.vertex_properties["bar"] = g.new_vertex_property("python::object")
        >>> g.graph_properties["gnat"] = g.new_graph_property("string", "hi there!")
        >>> g.list_properties()
        gnat           (graph)   (type: string, val: hi there!)
        bar            (vertex)  (type: python::object)
        foo            (vertex)  (type: double)
        foo            (edge)    (type: vector<double>)
        """

        if len(self.__properties) == 0:
            return
        w = max([len(x[0]) for x in self.__properties.keys()]) + 4
        w = w if w > 14 else 14

        for k, v in self.__properties.iteritems():
            if k[0] == "g":
                print "%%-%ds (graph)   (type: %%s, val: %%s)" % w % \
                      (k[1], v.value_type(), str(v[self]))
        for k, v in self.__properties.iteritems():
            if k[0] == "v":
                print "%%-%ds (vertex)  (type: %%s)" % w % (k[1],
                                                            v.value_type())
        for k, v in self.__properties.iteritems():
            if k[0] == "e":
                print "%%-%ds (edge)    (type: %%s)" % w % (k[1],
                                                            v.value_type())

    # index properties

    def _get_vertex_index(self):
        return self.__vertex_index
    vertex_index = property(_get_vertex_index,
                            doc="Vertex index map. This map is immutable.")

    def _get_edge_index(self):
        return self.__edge_index
    edge_index = property(_get_edge_index, doc="Edge index map.")

    # Property map creation

    def new_property(self, key_type, value_type):
        """Create a new (uninitialized) vertex property map of key type
        ``key_type`` (``v``, ``e`` or ``g``), value type ``value_type``, and
        return it.
        """
        if key_type == "v" or key_type == "vertex":
            return self.new_vertex_property(value_type)
        if key_type == "e" or key_type == "edge":
            return self.new_edge_property(value_type)
        if key_type == "g" or key_type == "graph":
            return self.new_graph_property(value_type)
        raise ValueError("unknown key type: " + key_type)

    def new_vertex_property(self, value_type):
        """Create a new (uninitialized) vertex property map of type
        ``value_type``, and return it."""
        return PropertyMap(new_vertex_property(_type_alias(value_type),
                                               self.__graph.GetVertexIndex()),
                           self, "v")

    def new_edge_property(self, value_type):
        """Create a new (uninitialized) edge property map of type
        ``value_type``, and return it."""
        return PropertyMap(new_edge_property(_type_alias(value_type),
                                             self.__graph.GetEdgeIndex()),
                           self, "e")

    def new_graph_property(self, value_type, val=None):
        """Create a new graph property map of type ``value_type``, and return
        it. If ``val`` is not None, the property is initialized to its value."""
        prop = PropertyMap(new_graph_property(_type_alias(value_type),
                                              self.__graph.GetGraphIndex()),
                           self, "g", lambda k: k.__graph)
        if val != None:
            prop[self] = val
        return prop

    # property map copying
    @_require("src", PropertyMap)
    @_require("tgt", (PropertyMap, type(None)))
    def copy_property(self, src, tgt=None, value_type=None, g=None):
        """Copy contents of ``src`` property to ``tgt`` property. If ``tgt`` is
        None, then a new property map of the same type (or with the type given
        by the optional ``value_type`` parameter) is created, and returned. The
        optional parameter g specifies the (identical) source graph to copy
        properties from (defaults to self).
        """
        if tgt == None:
            tgt = self.new_property(src.key_type(),
                                    (src.value_type()
                                     if value_type == None else value_type))
            ret = tgt
        else:
            ret = None

        if src.key_type() != tgt.key_type():
            raise ValueError("source and target properties must have the same" +
                             " key type")
        if g == None:
            g = self
        if g != self:
            g.stash_filter()
        self.stash_filter()
        if src.key_type() == "v":
            self.__graph.CopyVertexProperty(g.__graph, _prop("v", g, src),
                                            _prop("v", g, tgt))
        elif src.key_type() == "e":
            self.__graph.CopyEdgeProperty(g.__graph, _prop("e", g, src),
                                            _prop("e", g, tgt))
        else:
            tgt[self] = src[g]
        self.pop_filter()
        if g != self:
            g.pop_filter()
        return ret

    # degree property map
    @_limit_args({"deg": ["in", "out", "total"]})
    def degree_property_map(self, deg):
        """Create and return a vertex property map containing the degree type
        given by ``deg``."""
        return PropertyMap(self.__graph.DegreeMap(deg), self, "v")

    # I/O operations
    # ==============
    def load(self, file_name, file_format="auto"):
        """Load graph from ``file_name`` (which can be either a string or a
        file-like object). The format is guessed from ``file_name``, or can be
        specified by ``file_format``, which can be either "xml" or "dot". """

        if type(file_name) == str:
            file_name = os.path.expanduser(file_name)
        if file_format == 'auto' and isinstance(file_name, str):
            if file_name.endswith(".xml") or file_name.endswith(".xml.gz") or \
                   file_name.endswith(".xml.bz2"):
                file_format = "xml"
            elif file_name.endswith(".dot") or file_name.endswith(".dot.gz") or \
                     file_name.endswith(".dot.bz2"):
                file_format = "dot"
            else:
                raise ValueError("cannot determine file format of: " + file_name)
        elif file_format == "auto":
            file_format = "xml"
        if isinstance(file_name, str):
            props = self.__graph.ReadFromFile(file_name, None, file_format)
        else:
            props = self.__graph.ReadFromFile("", file_name, format)
        for name, prop in props[0].iteritems():
            self.vertex_properties[name] = PropertyMap(prop, self, "v")
        for name, prop in props[1].iteritems():
            self.edge_properties[name] = PropertyMap(prop, self, "e")
        for name, prop in props[2].iteritems():
            self.graph_properties[name] = PropertyMap(prop, self, "g",
                                                      lambda k: k.__graph)

    def save(self, file_name, file_format="auto"):
        """Save graph to ``file_name`` (which can be either a string or a
        file-like object). The format is guessed from the ``file_name``, or can
        be specified by ``file_format``, which can be either "xml" or "dot". """

        if type(file_name) == str:
            file_name = os.path.expanduser(file_name)
        if file_format == 'auto' and isinstance(file_name, str):
            if file_name.endswith(".xml") or file_name.endswith(".xml.gz") or \
                   file_name.endswith(".xml.bz2"):
                file_format = "xml"
            elif file_name.endswith(".dot") or file_name.endswith(".dot.gz") or \
                     file_name.endswith(".dot.bz2"):
                file_format = "dot"
            else:
                raise ValueError("cannot determine file file_format of: " + file_name)
        elif file_format == "auto":
            file_format = "xml"
        props = [(name[1], prop._PropertyMap__map) for name, prop in \
                 self.__properties.iteritems()]
        if isinstance(file_name, str):
            self.__graph.WriteToFile(file_name, None, file_format, props)
        else:
            self.__graph.WriteToFile("", file_name, file_format, props)

    # Directedness
    # ============

    def set_directed(self, is_directed):
        """Set the directedness of the graph."""
        self.__graph.SetDirected(is_directed)

    def is_directed(self):
        """Get the directedness of the graph."""
        return self.__graph.GetDirected()

    # Reversedness
    # ============

    def set_reversed(self, is_reversed):
        """Reverse the direction of the edges, if ``reversed`` is ``True``, or
        maintain the original direction otherwise."""
        self.__graph.SetReversed(is_reversed)

    def is_reversed(self):
        """Return ``True`` if the edges are reversed, and ``False`` otherwise.
        """
        return self.__graph.GetReversed()

    # Filtering
    # =========

    def set_vertex_filter(self, prop, inverted=False):
        """Choose vertex boolean filter property. Only the vertices with value
        different than zero are kept in the filtered graph. If the ``inverted``
        option is supplied with value ``True``, only the vertices with value
        zero are kept. If the supplied property is ``None``, any previous
        filtering is removed."""

        self.__graph.SetVertexFilterProperty(_prop("v", self, prop),
                                             inverted)
        self.__filter_state["vertex_filter"] = (prop, inverted)

    def get_vertex_filter(self):
        """Return a tuple with the vertex filter property and bool value
        indicating whether or not it is inverted."""
        return self.__filter_state["vertex_filter"]

    def set_edge_filter(self, prop, inverted=False):
        """Choose edge boolean filter property. Only the edges with value
        different than zero are kept in the filtered graph. If the ``inverted``
        option is supplied with value ``True``, only the edges with value zero
        are kept. If the supplied property is ``None``, any previous filtering
        is removed."""
        self.__graph.SetEdgeFilterProperty(_prop("e", self, prop), inverted)
        self.__filter_state["edge_filter"] = (prop, inverted)

    def get_edge_filter(self):
        """Return a tuple with the edge filter property and bool value
        indicating whether or not it is inverted."""
        return self.__filter_state["edge_filter"]

    def purge_vertices(self):
        """Remove all vertices of the graph which are currently being filtered
        out, and return it to the unfiltered state."""
        self.__graph.PurgeVertices()
        self.set_vertex_filter(None)

    def purge_edges(self):
        """Remove all edges of the graph which are currently being filtered out,
        and return it to the unfiltered state."""
        self.__graph.PurgeEdges()
        self.set_edge_filter(None)

    def stash_filter(self, edge=False, vertex=False, directed=False,
                     reversed=False, all=True):
        """Stash current filter state and set the graph to its unfiltered
        state. The optional keyword arguments specify which type of filter
        should be stashed."""
        if edge or vertex or directed or reversed:
            all = False
        self.__stashed_filter_state.append(self.get_filter_state())
        if libcore.graph_filtering_enabled():
            if vertex or all:
                self.set_vertex_filter(None)
            if edge or all:
                self.set_edge_filter(None)
        if directed or all:
            self.set_directed(True)
        if reversed or all:
            self.set_reversed(False)

    def pop_filter(self, edge=False, vertex=False, directed=False,
                   reversed=False, all=True):
        """Pop last stashed filter state. The optional keyword arguments specify
        which type of filter should be recovered."""
        if edge or vertex or directed or reversed:
            all = False
        state = self.__stashed_filter_state.pop()
        if libcore.graph_filtering_enabled():
            if vertex or all:
                self.set_vertex_filter(state["vertex_filter"][0],
                                       state["vertex_filter"][1])
            if edge or all:
                self.set_edge_filter(state["edge_filter"][0],
                                     state["edge_filter"][1])
        if directed or all:
            self.set_directed(state["directed"])
        if reversed or all:
            self.set_reversed(state["reversed"])

    def get_filter_state(self):
        """Return a copy of the filter state of the graph."""
        self.__filter_state["directed"] = self.is_directed()
        self.__filter_state["reversed"] = self.is_reversed()
        return copy.copy(self.__filter_state)

    def set_filter_state(self, state):
        """Set the filter state of the graph."""
        if libcore.graph_filtering_enabled():
            self.set_vertex_filter(state["vertex_filter"][0],
                                   state["vertex_filter"][1])
            self.set_edge_filter(state["edge_filter"][0],
                                 state["edge_filter"][1])
        self.set_directed(state["directed"])
        self.set_reversed(state["reversed"])

    # Basic graph statistics
    # ======================

    def num_vertices(self):
        """Get the number of vertices."""
        return self.__graph.GetNumberOfVertices()

    def num_edges(self):
        """Get the number of edges."""
        return self.__graph.GetNumberOfEdges()

    # Pickling support
    # ================

    def __getstate__(self):
        state = dict()
        sio = StringIO()
        if libcore.graph_filtering_enabled():
            if self.get_vertex_filter()[0] != None:
                self.vertex_properties["_Graph__pickle__vfilter"] = \
                    self.get_vertex_filter()[0]
                state["vfilter"] = self.get_vertex_filter()[1]
            if self.get_edge_filter()[0] != None:
                self.edge_properties["_Graph__pickle__efilter"] = \
                    self.get_edge_filter()[0]
                state["efilter"] = self.get_edge_filter()[1]
        self.save(sio, "xml")
        state["blob"] = sio.getvalue()
        return state

    def __setstate__(self, state):
        self.__init__()
        blob = state["blob"]
        if blob != "":
            sio = StringIO(blob)
            self.load(sio, "xml")
        if "vfilt" in state:
            vprop = self.vertex_properties["_Graph__pickle__vfilter"]
            self.set_vertex_filter(vprop, state["vfilt"])
        if "efilt" in state:
            eprop = self.edge_properties["_Graph__pickle__efilter"]
            self.set_edge_filter(eprop, state["efilt"])


def load_graph(file_name, file_format="auto"):
    """Load a graph from ``file_name`` (which can be either a string or a
        file-like object). The format is guessed from ``file_name``, or can be
        specified by ``file_format``, which can be either "xml" or "dot"."""
    g = Graph()
    g.load(file_name, file_format=file_format)
    return g


def value_types():
    """Return a list of possible properties value types."""
    return libcore.get_property_types()

# Vertex and Edge Types
# =====================
from libgraph_tool_core import Vertex, Edge


def _out_neighbours(self):
    """Return an iterator over the out-neighbours."""
    for e in self.out_edges():
        yield e.target()
Vertex.out_neighbours = _out_neighbours


def _in_neighbours(self):
    """Return an iterator over the in-neighbours."""
    for e in self.in_edges():
        yield e.source()
Vertex.in_neighbours = _in_neighbours


def _all_edges(self):
    """Return an iterator over all edges (both in or out)."""
    for e in self.out_edges():
        yield e
    for e in self.in_edges():
        yield e
Vertex.all_edges = _all_edges


def _all_neighbours(self):
    """Return an iterator over all neighbours (both in or out)."""
    for v in self.out_neighbours():
        yield v
    for v in self.in_neighbours():
        yield v
Vertex.all_neighbours = _all_neighbours


def _vertex_repr(self):
    if not self.is_valid():
        return "<invalid Vertex object at 0x%x>" % (id(self))
    return "<Vertex object with index '%d' at 0x%x>" % (int(self), id(self))
Vertex.__repr__ = _vertex_repr


def _edge_iter(self):
    """Iterate over the source and target"""
    for v in [self.source(), self.target()]:
        yield v


def _edge_repr(self):
    if not self.is_valid():
        return "<invalid Edge object at 0x%x>" % (id(self))

    return ("<Edge object with source '%d' and target '%d'" +
            " at 0x%x>") % (int(self.source()), int(self.target()), id(self))


# There are several edge classes... me must cycle through them all to modify
# them.
def init_edge_classes():
    for directed in [True, False]:
        for e_reversed in [True, False]:
            for e_filtered in [True, False]:
                for v_filtered in [True, False]:
                    g = Graph(directed=directed)
                    g.set_reversed(e_reversed)
                    v = g.add_vertex()
                    g.add_edge(v, v)
                    if e_filtered:
                        e_filter = g.new_edge_property("bool")
                        e_filter.a = [1]
                        g.set_edge_filter(e_filter)
                    if v_filtered:
                        v_filter = g.new_vertex_property("bool")
                        v_filter.a = [1]
                        g.set_vertex_filter(v_filter)
                    e = g.edges().next()
                    e.__class__.__repr__ = _edge_repr
                    e.__class__.__iter__ = _edge_iter

init_edge_classes()
