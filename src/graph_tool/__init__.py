#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# graph_tool -- a general graph manipulation python module
#
# Copyright (C) 2007-2011 Tiago de Paula Peixoto <tiago@skewed.de>
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

"""
graph_tool - efficient graph analysis and manipulation
======================================================

Summary
-------

.. autosummary::
   :nosignatures:

   Graph
   GraphView
   Vertex
   Edge
   PropertyMap
   PropertyArray
   load_graph
   group_vector_property
   ungroup_vector_property
   value_types
   show_config


This module provides:

   1. A :class:`~graph_tool.Graph` class for graph representation and manipulation
   2. Property maps for Vertex, Edge or Graph.
   3. Fast algorithms implemented in C++.

How to use the documentation
----------------------------

Documentation is available in two forms: docstrings provided
with the code, and the full documentation available in
`the graph-tool homepage <http://graph-tool.skewed.de>`_.

We recommend exploring the docstrings using `IPython
<http://ipython.scipy.org>`_, an advanced Python shell with TAB-completion and
introspection capabilities.

The docstring examples assume that ``graph_tool.all`` has been imported as
``gt``::

   >>> import graph_tool.all as gt

Code snippets are indicated by three greater-than signs::

   >>> x = x + 1

Use the built-in ``help`` function to view a function's docstring::

   >>> help(gt.Graph)

Contents
--------
"""

__author__ = "Tiago de Paula Peixoto <tiago@skewed.de>"
__copyright__ = "Copyright 2007-2011 Tiago de Paula Peixoto"
__license__ = "GPL version 3 or above"
__URL__ = "http://graph-tool.skewed.de"

# import numpy and scipy before everything to avoid weird segmentation faults
# depending on the order things are imported.

import numpy
import numpy.ma
import scipy
import scipy.stats

from dl_import import *
dl_import("import libgraph_tool_core as libcore")
import libgraph_tool_core as libcore   # for pylint
__version__ = libcore.mod_info().version

import io  # sets up libcore io routines

import sys
import os
import re
import gzip
import bz2
import weakref
import copy

from StringIO import StringIO
from decorators import _wraps, _require, _attrs, _limit_args
from inspect import ismethod

__all__ = ["Graph", "GraphView", "Vertex", "Edge", "Vector_bool",
           "Vector_int32_t", "Vector_int64_t", "Vector_double",
           "Vector_long_double", "Vector_string", "value_types", "load_graph",
           "PropertyMap", "group_vector_property", "ungroup_vector_property",
           "show_config", "PropertyArray", "__author__", "__copyright__",
           "__URL__", "__version__"]

# this is rather pointless, but it works around a sphinx bug
graph_tool = sys.modules[__name__]

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
             "long": "int64_t",
             "long long": "int64_t",
             "object": "python::object",
             "float": "double"}
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


def _python_type(type_name):
    type_name = _type_alias(type_name)
    if "vector" in type_name:
        ma = re.compile(r"vector<(.*)>").match(type_name)
        t = ma.group(1)
        return list, _python_type(t)
    if "int" in type_name:
        return int
    if type_name == "bool":
        return bool
    if "double" in type_name:
        return float
    if "string" in type_name:
        return str
    return object


def _gt_type(obj):
    t = type(obj)
    if t is numpy.longlong or t is numpy.uint64:
        return "long long"
    if t is int or issubclass(t, numpy.int):
        return "int"
    if t is numpy.float128:
        return "long double"
    if t is float or issubclass(t, numpy.float):
        return "double"
    if t is str:
        return "string"
    if t is bool:
        return "bool"
    if issubclass(t, list) or issubclass(t, numpy.ndarray):
        return "vector<%s>" % _gt_type(obj[0])
    return "object"


def _convert(prop, val):
    # attempt to convert to a compatible python type. This is useful,
    # for instance, when dealing with numpy types.
    vtype = _python_type(prop.value_type())
    if type(vtype) is tuple:
        return [vtype[1](x) for x in val]
    return vtype(val)


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


class PropertyArray(numpy.ndarray):
    """This is a :class:`~numpy.ndarray` subclass which keeps a reference of its :class:`~graph_tool.PropertyMap` owner, and detects if the underlying data has been invalidated."""

    __array_priority__ = -10

    def _get_pmap(self):
        return self._prop_map

    def _set_pmap(self, value):
        self._prop_map = value

    prop_map = property(_get_pmap, _set_pmap,
                        doc=":class:`~graph_tool.PropertyMap` owner instance.")

    def __new__(cls, input_array, prop_map):
        obj = numpy.asarray(input_array).view(cls)
        obj.prop_map = prop_map

        # check if data really belongs to property map
        if (prop_map._get_data().__array_interface__['data'][0] !=
            obj._get_base_data()):
            obj.prop_map = None
            # do a copy
            obj = numpy.asarray(obj)

        return obj

    def _get_base(self):
        base = self
        while base.base is not None:
            base = base.base
        return base

    def _get_base_data(self):
        return self._get_base().__array_interface__['data'][0]

    def _check_data(self):
        if self.prop_map is None:
            return

        data = self.prop_map._get_data()

        if (data is None or
            data.__array_interface__['data'][0] != self._get_base_data()):
            raise ValueError(("The graph correspondig to the underlying" +
                              " property map %s has changed. The" +
                              " PropertyArray at 0x%x is no longer valid!") %
                             (repr(self.prop_map), id(self)))

    def __array_finalize__(self, obj):
        if type(obj) is PropertyArray:
            obj._check_data()

        if obj is not None:
            # inherit prop_map only if the data is the same
            if (type(obj) is PropertyArray and
                self._get_base_data() == obj._get_base_data()):
                self.prop_map = getattr(obj, 'prop_map', None)
            else:
                self.prop_map = None
        self._check_data()

    def __array_prepare__(self, out_arr, context=None):
        self._check_data()
        return numpy.ndarray.__array_prepare__(self, out_arr, context)

    def __array_wrap__(self, out_arr, context=None):
        #demote to ndarray
        obj = numpy.ndarray.__array_wrap__(self, out_arr, context)
        return numpy.asarray(obj)

    # Overload members and operators to add data checking

    def _wrap_method(method):
        method = getattr(numpy.ndarray, method)

        def checked_method(self, *args, **kwargs):
            self._check_data()
            return method(self, *args, **kwargs)

        if ismethod(method):
            checked_method = _wraps(method)(checked_method)
        checked_method.__doc__ = getattr(method, "__doc__", None)
        return checked_method

    for method in ['all', 'any', 'argmax', 'argmin', 'argsort', 'astype',
                   'byteswap', 'choose', 'clip', 'compress', 'conj',
                   'conjugate', 'copy', 'cumprod', 'cumsum', 'diagonal', 'dot',
                   'dump', 'dumps', 'fill', 'flat', 'flatten', 'getfield',
                   'imag', 'item', 'itemset', 'itemsize', 'max', 'mean', 'min',
                   'newbyteorder', 'nonzero', 'prod', 'ptp', 'put', 'ravel',
                   'real', 'repeat', 'reshape', 'resize', 'round',
                   'searchsorted', 'setfield', 'setflags', 'sort', 'squeeze',
                   'std', 'sum', 'swapaxes', 'take', 'tofile', 'tolist',
                   'tostring', 'trace', 'transpose', 'var', 'view',
                   '__getitem__']:
        if hasattr(numpy.ndarray, method):
            locals()[method] = _wrap_method(method)


class PropertyMap(object):
    """This class provides a mapping from vertices, edges or whole graphs to arbitrary properties.

    See :ref:`sec_property_maps` for more details.

    The possible property value types are listed below.

    .. table::

        =======================     ======================
         Type name                  Alias
        =======================     ======================
        ``bool``                    ``uint8_t``
        ``int32_t``                 ``int``
        ``int64_t``                 ``long``, ``long long``
        ``double``                  ``float``
        ``long double``
        ``string``
        ``vector<bool>``            ``vector<uint8_t>``
        ``vector<int32_t>``         ``vector<int>``
        ``vector<int64_t>``         ``vector<long>``, ``vector<long long>``
        ``vector<double>``          ``vector<float>``
        ``vector<long double>``
        ``vector<string>``
        ``python::object``          ``object``
        =======================     ======================
    """
    def __init__(self, pmap, g, key_type, key_trans=None):
        self.__map = pmap
        self.__g = weakref.ref(g)
        self.__base_g = lambda: None
        try:
            if isinstance(g, GraphView):
                self.__base_g = weakref.ref(g.base)  # keep reference to the
                                                     # base graph, in case the
                                                     # graph view is deleted.
        except NameError:
            pass  # ignore if GraphView is yet undefined
        self.__key_type = key_type
        self.__key_trans = key_trans if key_trans is not None else lambda k: k
        self.__register_map()

    def __register_map(self):
        for g in [self.__g(), self.__base_g()]:
            if g is not None:
                g._Graph__known_properties.append((self.key_type(),
                                                   weakref.ref(self.__map)))

    def __unregister_map(self):
        for g in [self.__g(), self.__base_g()]:
            if g is not None:
                i = g._Graph__known_properties.index((self.key_type(),
                                                      weakref.ref(self.__map)))
                del g._Graph__known_properties[i]

    def __del__(self):
        self.__unregister_map()

    def __getitem__(self, k):
        return self.__map[self.__key_trans(k)]

    def __setitem__(self, k, v):
        key = self.__key_trans(k)
        try:
            self.__map[key] = v
        except TypeError:
            self.__map[key] = _convert(self, v)

    def __repr__(self):
        # provide some more useful information
        if self.key_type() == "e":
            k = "Edge"
        elif self.key_type() == "v":
            k = "Vertex"
        else:
            k = "Graph"
        g = self.get_graph()
        if g == None:
            g = "a non-existent graph"
        else:
            g = "Graph 0x%x" % id(g)
        return ("<PropertyMap object with key type '%s' and value type '%s',"
                + " for %s, at 0x%x>") % (k, self.value_type(), g, id(self))

    def copy(self, value_type=None):
        """Return a copy of the property map. If ``value_type`` is specified,
        the value type is converted to the chosen type."""
        return self.get_graph().copy_property(self, value_type=value_type)

    def get_graph(self):
        """Get the graph class to which the map refers."""
        g = self.__g()
        if g is None:
            g = self.__base_g()
        return g

    def key_type(self):
        """Return the key type of the map. Either 'g', 'v' or 'e'."""
        return self.__key_type

    def value_type(self):
        """Return the value type of the map."""
        return self.__map.value_type()

    def python_value_type(self):
        """Return the python-compatible value type of the map."""
        return _python_type(self.__map.value_type())

    def get_array(self):
        """Get a :class:`~graph_tool.PropertyArray` with the property values.

        .. note::

           An array is returned *only if* the value type of the property map is
           a scalar. For vector, string or object types, ``None`` is returned
           instead.

        .. warning::

           The returned array does not own the data, which belongs to the
           property map. Therefore, if the graph changes, the array may become
           *invalid* and any operation on it will fail with a
           :class:`ValueError` exception. Do **not** store the array if
           the graph is to be modified; store a **copy** instead.
        """
        a = self._get_data()
        if a is None:
            return None
        return PropertyArray(a, prop_map=self)

    def _get_data(self):
        g = self.get_graph()
        if g is None:
            return None
        g.stash_filter(edge=True, vertex=True)
        if self.__key_type == 'v':
            n = g.num_vertices()
        elif self.__key_type == 'e':
            n = g._Graph__graph.GetMaxEdgeIndex() + 1
        else:
            n = 1
        g.pop_filter(edge=True, vertex=True)
        a = self.__map.get_array(n)
        return a

    def __set_array(self, v):
        a = self.get_array()
        if a is None:
            return
        a[:] = v

    a = property(get_array, __set_array,
                 doc=r"""Shortcut to the :meth:`~PropertyMap.get_array` method
                 as an attribute. This makes assignments more convenient, e.g.:

                 >>> g = gt.Graph()
                 >>> g.add_vertex(10)
                 [...]
                 >>> prop = g.new_vertex_property("double")
                 >>> prop.a = np.random.random(10)           # Assignment from array
                 """)

    def __get_set_f_array(self, v=None, get=True):
        g = self.get_graph()
        if g is None:
            return None
        a = self.get_array()
        filt = [None]
        if self.__key_type == 'v':
            filt = g.get_vertex_filter()
        elif self.__key_type == 'e':
            filt = g.get_edge_filter()
        if get:
            if a is None:
                return a
            if filt[0] is None:
                return a
            return a[filt[0].a == (not filt[1])]
        else:
            if a is None:
                return
            if filt[0] is None:
                a[:] = v
            else:
                a[filt[0].a == (not filt[1])] = v

    fa = property(__get_set_f_array,
                  lambda self, v: self.__get_set_f_array(v, False),
                  doc=r"""The same as the :attr:`~PropertyMap.a` attribute, but
                  instead an *indexed* array is returned, which contains only
                  entries for vertices/edges which are not filtered out. If
                  there are no filters in place, the array is not indexed, and
                  is identical to the :attr:`~PropertyMap.a` attribute.

                  Note that because advanced indexing is triggered, a **copy**
                  of the array is returned, not a view, as for the
                  :attr:`~PropertyMap.a` attribute. Nevertheless, the assignment
                  of values to the *whole* array at once works as expected.""")

    def __get_set_m_array(self, v=None, get=True):
        g = self.get_graph()
        if g is None:
            return None
        a = self.get_array()
        filt = [None]
        if self.__key_type == 'v':
            filt = g.get_vertex_filter()
        elif self.__key_type == 'e':
            filt = g.get_edge_filter()
        if filt[0] is None or a is None:
            if get:
                return a
            else:
                return
        ma = numpy.ma.array(a, mask=(filt[0].a == False) if not filt[1] else (filt[0].a == True))
        if get:
            return ma
        else:
            ma[:] = v

    ma = property(__get_set_m_array,
                  lambda self, v: self.__get_set_m_array(v, False),
                  doc=r"""The same as the :attr:`~PropertyMap.a` attribute, but
                  instead a :class:`~numpy.ma.MaskedArray` object is returned,
                  which contains only entries for vertices/edges which are not
                  filtered out. If there are no filters in place, a regular
                  :class:`~graph_tool.PropertyArray` is returned, which is
                  identical to the :attr:`~PropertyMap.a` attribute.""")

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


def group_vector_property(props, value_type=None, vprop=None, pos=None):
    """Group list of properties ``props`` into a vector property map of the same type.

    Parameters
    ----------
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

    Examples
    --------
    >>> from numpy.random import seed, randint
    >>> from numpy import array
    >>> seed(42)
    >>> g = gt.random_graph(100, lambda: (3, 3))
    >>> props = [g.new_vertex_property("int") for i in xrange(3)]
    >>> for i in xrange(3):
    ...    props[i].a = randint(0, 100, g.num_vertices())
    >>> gprop = gt.group_vector_property(props)
    >>> print gprop[g.vertex(0)].a
    [71 40 96]
    >>> print array([p[g.vertex(0)] for p in props])
    [71 40 96]
    """
    g = props[0].get_graph()
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
            libcore.group_vector_property(g._Graph__graph, _prop(k, g, vprop),
                                          _prop(k, g, p),
                                          i if pos == None else pos[i],
                                          k == 'e')
            g.pop_filter(directed=True, reversed=True)
        else:
            vprop[g][i if pos is None else pos[i]] = p[g]
    return vprop


def ungroup_vector_property(vprop, pos, props=None):
    """Ungroup vector property map ``vprop`` into a list of individual property maps.

    Parameters 
    ----------
    vprop : :class:`~graph_tool.PropertyMap`
        Vector property map to be ungrouped.
    pos : list of ints
        A list of indexes corresponding to where each element of ``vprop``
        should be inserted into the ungrouped list.
    props : list of :class:`~graph_tool.PropertyMap`  (optional, default: None)
        If supplied, should contain a list of property maps to which ``vprop``
        should be ungroupped.

    Returns
    -------
    props : list of :class:`~graph_tool.PropertyMap`
       A list of property maps with the ungrouped values of ``vprop``.

    Examples
    --------
    >>> from numpy.random import seed, randint
    >>> from numpy import array
    >>> seed(42)
    >>> g = gt.random_graph(100, lambda: (3, 3))
    >>> prop = g.new_vertex_property("vector<int>")
    >>> for v in g.vertices():
    ...    prop[v] = randint(0, 100, 3)
    >>> uprops = gt.ungroup_vector_property(prop, [0, 1, 2])
    >>> print prop[g.vertex(0)].a
    [71 60 20]
    >>> print array([p[g.vertex(0)] for p in uprops])
    [71 60 20]
    """

    g = vprop.get_graph()
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
    """Generic multigraph class.

    This class encapsulates either a directed multigraph (default or if
    ``directed=True``) or an undirected multigraph (if ``directed=False``),
    with optional internal edge, vertex or graph properties.

    If ``g`` is specified, the graph (and its internal properties) will be
    copied.

    If ``prune`` is set to True, and ``g`` is specified, only the filtered graph
    will be copied, and the new graph object will not be filtered. Optionally, a
    tuple of three booleans can be passed as value to ``prune``, to specify a
    different behavior to vertex, edge, and reversal filters, respectively.

    The graph is implemented as an `adjacency list`_, where both vertex and edge
    lists are C++ STL vectors.

    .. _adjacency list: http://en.wikipedia.org/wiki/Adjacency_list

    """

    def __init__(self, g=None, directed=True, prune=False):
        self.__properties = {}
        self.__known_properties = []
        self.__filter_state = {"reversed": False,
                               "edge_filter": (None, False),
                               "vertex_filter": (None, False),
                               "directed": True}
        self.__stashed_filter_state = []

        if g is None:
            self.__graph = libcore.GraphInterface()
            self.set_directed(directed)
        else:
            if isinstance(prune, bool):
                vprune = eprune = rprune = prune
            else:
                vprune, eprune, rprune = prune
            if not (vprune or eprune or rprune):
                g.stash_filter(vertex=vprune, edge=vprune, reversed=rprune)
            self.__graph = libcore.GraphInterface(g.__graph, False)
            if not (vprune or eprune or rprune):
                g.pop_filter(vertex=vprune, edge=vprune, reversed=rprune)

            for k, v in g.__properties.iteritems():
                new_p = self.new_property(v.key_type(), v.value_type())
                self.copy_property(v, new_p, g=g)
                self.properties[k] = new_p

            self.__stashed_filter_state = [self.get_filter_state()]

            if not vprune:
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
            if not eprune:
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
            if not rprune:
                self.__stashed_filter_state[0]["reversed"] = g.is_reversed()

            # directedness is always a filter
            self.__stashed_filter_state[0]["directed"] = g.is_directed()

            self.pop_filter()

            if vprune or eprune:
                self.reindex_edges()

        # internal index maps
        self.__vertex_index = \
                 PropertyMap(libcore.get_vertex_index(self.__graph), self, "v")
        self.__edge_index = \
                 PropertyMap(libcore.get_edge_index(self.__graph), self, "e")

        # modification permissions
        self.__perms = {"add_edge": True, "del_edge": True,
                        "add_vertex": True, "del_vertex": True}

    def copy(self):
        """Return a deep copy of self. All :ref:`internal property maps <sec_internal_props>`
        are also copied."""
        return Graph(self)

    def __repr__(self):
        # provide more useful information
        d = "directed" if self.is_directed() else "undirected"
        fr = ", reversed" if self.is_reversed() and self.is_directed() else ""
        f = ""
        if self.get_edge_filter()[0] is not None:
            f += ", edges filtered by %s" % (str(self.get_edge_filter()))
        if self.get_vertex_filter()[0] is not None:
            f += ", vertices filtered by %s" % (str(self.get_vertex_filter()))
        n = self.num_vertices()
        e = self.num_edges()
        return "<%s object, %s%s, with %d %s and %d edge%s%s at 0x%x>"\
               % (type(self).__name__, d, fr, n,
                  "vertex" if n == 1 else "vertices", e, "" if e == 1 else "s",
                  f, id(self))

    # Graph access
    # ============

    def __check_perms(self, ptype):
        if not self.__perms[ptype]:
            raise RuntimeError("the graph cannot be modified at this point!")

    def vertices(self):
        """Return an :meth:`iterator <iterator.__iter__>` over the vertices.

        .. note::

           The order of the vertices traversed by the iterator **always**
           corresponds to the vertex index ordering, as given by the
           :attr:`~graph_tool.Graph.vertex_index` property map.

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

    def vertex(self, i, use_index=True):
        """Return the vertex with index ``i``. If ``use_index=False``, the
        ``i``-th vertex is returned (which can differ from the vertex with index
        ``i`` in case of filtered graphs). """
        if use_index:
            self.stash_filter(vertex=True)
        try:
            v = libcore.get_vertex(weakref.ref(self.__graph), int(i))
        finally:
            if use_index:
                self.pop_filter(vertex=True)
        return v

    def edge(self, s, t, all_edges=False):
        """Return the edge from vertex ``s`` to ``t``, if it exists. If
        ``all_edges=True`` then a list is returned with all the parallel edges
        from ``s`` to ``t``, otherwise only one edge is returned.

        This operation will take :math:`O(k(s))` time, where :math:`k(s)` is the
        out-degree of vertex :math:`s`.
        """
        s = self.vertex(int(s))
        t = self.vertex(int(t))
        edges = []
        for e in s.out_edges():
            if e.target() == t:
                if not all_edges:
                    return e
                edges.append(e)
        if all_edges:
            return edges
        return None

    def edges(self):
        """Return an :meth:`iterator <iterator.__iter__>` over the edges.

        .. note::

           The order of the edges traversed by the iterator **does not**
           necessarily correspond to the edge index ordering, as given by the
           :attr:`~graph_tool.Graph.edge_index` property map. This will only
           happen after :meth:`~graph_tool.Graph.reindex_edges` is called, or in
           certain situations such as just after a graph is loaded from a
           file. However, further manipulation of the graph may destroy the
           ordering.

        """
        return libcore.get_edges(weakref.ref(self.__graph))

    def add_vertex(self, n=1):
        """Add a vertex to the graph, and return it. If ``n > 1``, ``n``
        vertices are inserted and a list is returned."""
        self.__check_perms("add_vertex")
        vlist = []
        vfilt = self.get_vertex_filter()
        for i in xrange(n):
            v = libcore.add_vertex(weakref.ref(self.__graph))
            if vfilt[0] is not None:
                vfilt[0][v] = not vfilt[1]
            vlist.append(v)
        if n == 1:
            return vlist[0]
        return vlist

    def remove_vertex(self, vertex):
        """Remove a vertex from the graph."""
        self.__check_perms("del_vertex")
        index = self.vertex_index[vertex]
        for pmap in self.__known_properties:
            if pmap[0] == "v" and pmap[1]() != None and pmap[1]().is_writable():
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
        self.__check_perms("add_edge")
        e = libcore.add_edge(weakref.ref(self.__graph), source, target)
        efilt = self.get_edge_filter()
        if efilt[0] is not None:
            efilt[0][e] = not efilt[1]
        return e

    def remove_edge(self, edge):
        """Remove an edge from the graph."""
        self.__check_perms("del_edge")
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
        self.__check_perms("del_vertex")
        self.__check_perms("del_edge")
        self.__graph.Clear()

    def clear_edges(self):
        """Remove all edges from the graph."""
        self.__check_perms("del_edge")
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
    vertex, edge or graph property, respectively, and the second element is the
    name of the property map.

    Examples
    --------
    >>> g = gt.Graph()
    >>> g.properties[("e", "foo")] = g.new_edge_property("vector<double>")
    >>> del g.properties[("e", "foo")]
    """)

    def __get_specific_properties(self, t):
        props = dict([(k[1], v) for k, v in self.__properties.iteritems() \
                      if k[0] == t])
        return props

    # vertex properties
    def __get_vertex_properties(self):
        return PropertyDict(self, self.__get_specific_properties("v"),
                            lambda g, k: g.__properties[("v", k)],
                            lambda g, k, v: g.__set_property("v", k, v),
                            lambda g, k: g.__del_property("v", k))
    vertex_properties = property(__get_vertex_properties,
                                 doc="Dictionary of internal vertex properties. The keys are the property names.")

    # edge properties
    def __get_edge_properties(self):
        return PropertyDict(self, self.__get_specific_properties("e"),
                            lambda g, k: g.__properties[("e", k)],
                            lambda g, k, v: g.__set_property("e", k, v),
                            lambda g, k: g.__del_property("e", k))
    edge_properties = property(__get_edge_properties,
                                 doc="Dictionary of internal edge properties. The keys are the property names.")

    # graph properties
    def __get_graph_properties(self):
        return PropertyDict(self, self.__get_specific_properties("g"),
                            lambda g, k: g.__properties[("g", k)],
                            lambda g, k, v: g.__set_property("g", k, v),
                            lambda g, k: g.__del_property("g", k))
    graph_properties = property(__get_graph_properties,
                                 doc="Dictionary of internal graph properties. The keys are the property names.")

    def own_property(self, prop):
        """'Own' property map 'prop', which may belong to another graph."""
        return PropertyMap(prop._PropertyMap__map, self, prop.key_type())

    def list_properties(self):
        """Print a list of all internal properties.

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
                            doc="""Vertex index map.

                            It maps for each vertex in the graph an unique
                            integer in the range [0, :meth:`~graph_tool.Graph.num_vertices` - 1].

                            .. note::

                                This is a special instance of a :class:`~graph_tool.PropertyMap`
                                class, which is **immutable**, and cannot be
                                accessed as an array.""")

    def _get_edge_index(self):
        return self.__edge_index
    edge_index = property(_get_edge_index, doc="""Edge index map.

                            It maps for each edge in the graph an unique
                            integer.

                            .. warning::

                                Differently from :attr:`~graph_tool.Graph.vertex_index`,
                                this is a **regular** instance of a :class:`~graph_tool.PropertyMap`
                                class, and is therefore **mutable**!

                                Additionally, the indexes may not necessarily
                                lie in the range [0, :meth:`~graph_tool.Graph.num_edges` - 1].
                                However this will always happen whenever no
                                edges are deleted from the graph.

                                The internal consistency expected by most
                                algorithms and the proper functioning of
                                property maps assume that the indexes are unique
                                and constant, which is guaranteed by the
                                library.  Therefore it is recommended **never**
                                to modify these values, unless you know what you
                                are doing.""")

    def _get_max_edge_index(self):
        return self.__graph.GetMaxEdgeIndex()
    max_edge_index = property(_get_max_edge_index,
                              doc="The maximum value of the edge index map.")

    def reindex_edges(self):
        """
        Reset the edge indexes so that they lie in the [0, :meth:`~graph_tool.Graph.num_edges` - 1]
        range. The index ordering will be compatible with the sequence returned
        by the :meth:`~graph_tool.Graph.edges` function.

        .. WARNING::

           Calling this function will invalidate all existing edge property
           maps, if the index ordering is modified! The property maps will still
           be usable, but their contents will still be tied to the old indexes,
           and thus may become scrambled.
        """
        self.__graph.ReIndexEdges()

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
        if tgt is None:
            tgt = self.new_property(src.key_type(),
                                    (src.value_type()
                                     if value_type == None else value_type))
            ret = tgt
        else:
            ret = None

        if src.key_type() != tgt.key_type():
            raise ValueError("source and target properties must have the same" +
                             " key type")
        if g is None:
            g = self
        if g is not self:
            g.stash_filter(directed=True, reversed=True)
        self.stash_filter(directed=True, reversed=True)
        if src.key_type() == "v":
            self.__graph.CopyVertexProperty(g.__graph, _prop("v", g, src),
                                            _prop("v", self, tgt))
        elif src.key_type() == "e":
            self.__graph.CopyEdgeProperty(g.__graph, _prop("e", g, src),
                                            _prop("e", self, tgt))
        else:
            tgt[self] = src[g]
        self.pop_filter(directed=True, reversed=True)
        if g is not self:
            g.pop_filter(directed=True, reversed=True)
        return ret

    # degree property map
    @_limit_args({"deg": ["in", "out", "total"]})
    def degree_property_map(self, deg):
        """Create and return a vertex property map containing the degree type
        given by ``deg``."""
        return PropertyMap(self.__graph.DegreeMap(deg), self, "v")

    # I/O operations
    # ==============
    def __get_file_format(self, file_name):
        fmt = None
        for f in ["xml", "dot", "gml"]:
            names = ["." + f, ".%s.gz" % f, ".%s.bz2" % f]
            for name in names:
                if file_name.endswith(name):
                    fmt = f
                    break
        if fmt is None:
            raise ValueError("cannot determine file format of: " + file_name)
        return fmt

    def load(self, file_name, fmt="auto"):
        """Load graph from ``file_name`` (which can be either a string or a
        file-like object). The format is guessed from ``file_name``, or can be
        specified by ``fmt``, which can be either "xml", "dot" or "gml". """

        if type(file_name) == str:
            file_name = os.path.expanduser(file_name)
        if fmt == 'auto' and isinstance(file_name, str):
            fmt = self.__get_file_format(file_name)
        elif fmt == "auto":
            fmt = "xml"
        if isinstance(file_name, str):
            props = self.__graph.ReadFromFile(file_name, None, fmt)
        else:
            props = self.__graph.ReadFromFile("", file_name, fmt)
        for name, prop in props[0].iteritems():
            self.vertex_properties[name] = PropertyMap(prop, self, "v")
        for name, prop in props[1].iteritems():
            self.edge_properties[name] = PropertyMap(prop, self, "e")
        for name, prop in props[2].iteritems():
            self.graph_properties[name] = PropertyMap(prop, self, "g",
                                                      lambda k: k.__graph)

    def save(self, file_name, fmt="auto"):
        """Save graph to ``file_name`` (which can be either a string or a
        file-like object). The format is guessed from the ``file_name``, or can
        be specified by ``fmt``, which can be either "xml", "dot" or "gml". """

        if type(file_name) == str:
            file_name = os.path.expanduser(file_name)
        if fmt == 'auto' and isinstance(file_name, str):
            fmt = self.__get_file_format(file_name)
        elif fmt == "auto":
            fmt = "xml"
        props = [(name[1], prop._PropertyMap__map) for name, prop in \
                 self.__properties.iteritems()]
        if isinstance(file_name, str):
            self.__graph.WriteToFile(file_name, None, fmt, props)
        else:
            self.__graph.WriteToFile("", file_name, fmt, props)

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
        """Reverse the direction of the edges, if ``is_reversed`` is ``True``,
        or maintain the original direction otherwise."""
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
        """Get the number of vertices.

        .. note::

            If the vertices are being filtered, this operation is
            :math:`O(N)`. Otherwise it is :math:`O(1)`.
        """
        return self.__graph.GetNumberOfVertices()

    def num_edges(self):
        """Get the number of edges.

        .. note::

            If the edges are being filtered, this operation is
            :math:`O(E)`. Otherwise it is :math:`O(1)`.
        """
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


def load_graph(file_name, fmt="auto"):
    """
    Load a graph from ``file_name`` (which can be either a string or a file-like object).

    The format is guessed from ``file_name``, or can be specified by
    ``fmt``, which can be either "xml", "dot" or "gml".
    """
    g = Graph()
    g.load(file_name, fmt)
    return g


def value_types():
    """Return a list of possible properties value types."""
    return libcore.get_property_types()

# Vertex and Edge Types
# =====================
from libgraph_tool_core import Vertex, Edge

Vertex.__doc__ = """Vertex descriptor.

This class represents a vertex in a :class:`~graph_tool.Graph` instance.

:class:`~graph_tool.Vertex` instances are hashable, and are convertible to
integers, corresponding to its index (see :attr:`~graph_tool.Graph.vertex_index`).
"""


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

_edge_doc = """Edge descriptor.

This class represents an edge in a :class:`~graph_tool.Graph`.

:class:`~graph_tool.Edge` instances are hashable, and are convertible to a
tuple, which contains the source and target vertices.
"""


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
                    e.__class__.__doc__ = _edge_doc

init_edge_classes()

# Add convenience function to vector classes


def _get_array_view(self):
    return self.get_array()[:]


def _set_array_view(self, v):
    self.get_array()[:] = v

vector_types = [Vector_bool, Vector_int32_t, Vector_int64_t, Vector_double,
                Vector_long_double]
for vt in vector_types:
    vt.a = property(_get_array_view, _set_array_view,
                    doc=r"""Shortcut to the `get_array` method as an attribute.""")
    vt.__repr__ = lambda self: self.a.__repr__()
Vector_string.a = None
Vector_string.get_array = lambda self: None


class GraphView(Graph):
    """
    A view of selected vertices or edges of another graph.

    This class uses shared data from another :class:`~graph_tool.Graph`
    instance, but allows for local filtering of vertices and/or edges, edge
    directionality or reversal. See :ref:`sec_graph_views` for more details and
    examples.

    The existence of a :class:`~graph_tool.GraphView` object does not affect the
    original graph, except if the graph view is modified (addition or removal of
    vertices or edges), in which case the modification is directly reflected in
    the original graph (and vice-versa), since they both point to the same
    underlying data. Because of this, instances of
    :class:`~graph_tool.PropertyMap` can be used interchangeably with a graph
    and its views.

    The argument ``g`` must be an instance of a :class:`~graph_tool.Graph`
    class. If specified, ``vfilt`` and ``efilt`` select which vertices and edges
    are filtered, respectively. These parameters can either be a
    boolean-valued :class:`~graph_tool.PropertyMap` or a
    :class:`~numpy.ndarray`, which specify which vertices/edges are selected, or
    an unary function which returns ``True`` if a given vertex/edge is to be
    selected, or ``False`` otherwise.

    The boolean parameter ``directed`` can be used to set the directionality of
    the graph view. If ``directed = None``, the directionality is inherited from
    ``g``.

    If ``reversed = True``, the direction of the edges is reversed.

    If ``vfilt`` or ``efilt`` is anything other than a
    :class:`~graph_tool.PropertyMap` instance, the instantiation running time is
    :math:`O(V)` and :math:`O(E)`, respectively. Otherwise, the running time is
    :math:`O(1)`.
    """

    def __init__(self, g, vfilt=None, efilt=None, directed=None,
                 reversed=False):
        self.__base = g if not isinstance(g, GraphView) else g.base
        Graph.__init__(self)
        # copy graph reference
        self._Graph__graph = libcore.GraphInterface(g._Graph__graph, True)
        self._Graph__properties = g._Graph__properties
        self._Graph__known_properties = g._Graph__known_properties

        # set already existent filters
        vf = g.get_vertex_filter()
        if vf[0] is not None:
            self.set_vertex_filter(vf[0], vf[1])
        ef = g.get_edge_filter()
        if ef[0] is not None:
            self.set_edge_filter(ef[0], ef[1])

        if vfilt is not None:
            if type(vfilt) is PropertyMap:
                self.set_vertex_filter(vfilt)
            else:
                vmap = self.new_vertex_property("bool")
                if issubclass(type(vfilt), numpy.ndarray):
                    vmap.a = vfilt
                else:
                    omap, inv = g.get_vertex_filter()
                    if omap is not None:
                        vmap.a = omap.a if not inv else omap.a ^ 1
                    for v in g.vertices():
                        vmap[v] = vfilt(v)
                self.set_vertex_filter(vmap)

        if efilt is not None:
            if type(efilt) is PropertyMap:
                self.set_edge_filter(efilt)
            else:
                emap = self.new_edge_property("bool")
                if issubclass(type(efilt), numpy.ndarray):
                    vmap.a = efilt
                else:
                    omap, inv = g.get_edge_filter()
                    if omap is not None:
                        emap.a = omap.a if not inv else omap.a ^ 1
                    for e in g.edges():
                        emap[e] = efilt(e)
                self.set_edge_filter(emap)

        if directed is not None:
            self.set_directed(directed)
        if reversed:
            self.set_reversed(not g.is_reversed())

    def __get_base(self):
        return self.__base
    base = property(__get_base, doc="Base graph.")
