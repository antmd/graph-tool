#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# graph_tool -- a general graph manipulation python module
#
# Copyright (C) 2006-2015 Tiago de Paula Peixoto <tiago@skewed.de>
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
   infect_vertex_property
   edge_difference
   perfect_prop_hash
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

from __future__ import division, absolute_import, print_function
import sys
if sys.version_info < (3,):
    range = xrange

__author__ = "Tiago de Paula Peixoto <tiago@skewed.de>"
__copyright__ = "Copyright 2006-2015 Tiago de Paula Peixoto"
__license__ = "GPL version 3 or above"
__URL__ = "http://graph-tool.skewed.de"

# import numpy and scipy before everything to avoid weird segmentation faults
# depending on the order things are imported.

import numpy
import numpy.ma
import scipy
import scipy.stats


from .dl_import import *
dl_import("from . import libgraph_tool_core as libcore")
__version__ = libcore.mod_info().version

from . import io  # sets up libcore io routines

import sys
import os
import re
import gzip
import weakref
import copy
import textwrap
import io

if sys.version_info < (3,):
    import StringIO

from .decorators import _wraps, _require, _attrs, _limit_args
from inspect import ismethod

__all__ = ["Graph", "GraphView", "Vertex", "Edge", "Vector_bool",
           "Vector_int16_t", "Vector_int32_t", "Vector_int64_t", "Vector_double",
           "Vector_long_double", "Vector_string", "value_types", "load_graph",
           "PropertyMap", "group_vector_property", "ungroup_vector_property",
           "infect_vertex_property", "edge_difference", "perfect_prop_hash",
           "seed_rng", "show_config", "PropertyArray", "openmp_enabled",
           "openmp_get_num_threads", "openmp_set_num_threads", "openmp_get_schedule",
           "openmp_set_schedule", "__author__", "__copyright__", "__URL__",
           "__version__"]

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
             "short": "int16_t",
             "int": "int32_t",
             "long": "int64_t",
             "long long": "int64_t",
             "unsigned long": "int64_t",
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
    if issubclass(t, numpy.int16):
        return "short"
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
    if vtype is object:
        return val
    return vtype(val)


def show_config():
    """Show ``graph_tool`` build configuration."""
    info = libcore.mod_info()
    print("version:", info.version)
    print("gcc version:", info.gcc_version)
    print("compilation flags:", info.cxxflags)
    print("install prefix:", info.install_prefix)
    print("python dir:", info.python_dir)
    print("graph filtering:", libcore.graph_filtering_enabled())
    print("openmp:", libcore.openmp_enabled())
    print("uname:", " ".join(os.uname()))

def terminal_size():
    import fcntl, termios, struct
    h, w, hp, wp = struct.unpack('HHHH',
        fcntl.ioctl(0, termios.TIOCGWINSZ,
        struct.pack('HHHH', 0, 0, 0, 0)))
    return w, h

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
        ``int16_t``                 ``short``
        ``int32_t``                 ``int``
        ``int64_t``                 ``long``, ``long long``
        ``double``                  ``float``
        ``long double``
        ``string``
        ``vector<bool>``            ``vector<uint8_t>``
        ``vector<int16_t>``         ``short``
        ``vector<int32_t>``         ``vector<int>``
        ``vector<int64_t>``         ``vector<long>``, ``vector<long long>``
        ``vector<double>``          ``vector<float>``
        ``vector<long double>``
        ``vector<string>``
        ``python::object``          ``object``
        =======================     ======================
    """
    def __init__(self, pmap, g, key_type):
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
        self.__register_map()

    def __key_trans(self, key):
        if self.key_type() == "g":
            return key._Graph__graph
        else:
            return key


    def __register_map(self):
        for g in [self.__g(), self.__base_g()]:
            if g is not None:
                g._Graph__known_properties[id(self)] = weakref.ref(self)

    def __unregister_map(self):
        for g in [self.__g(), self.__base_g()]:
            if g is not None and id(self) in g._Graph__known_properties:
                del g._Graph__known_properties[id(self)]

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
        try:
            vals = ", with values:\n%s" % str(self.fa)
        except ValueError:
            vals = ""
        return ("<PropertyMap object with key type '%s' and value type '%s',"
                + " for %s, at 0x%x%s>") % (k, self.value_type(), g, id(self),
                                            vals)

    def copy(self, value_type=None):
        """Return a copy of the property map. If ``value_type`` is specified,
        the value type is converted to the chosen type."""
        return self.get_graph().copy_property(self, value_type=value_type)

    def __copy__(self):
        return self.copy()

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
           instead. For vector and string objects, indirect array access is
           provided via the :func:`~graph_tool.PropertyMap.get_2d_array()` and
           :func:`~graph_tool.PropertyMap.set_2d_array()` member functions.

        .. warning::

           The returned array does not own the data, which belongs to the
           property map. Therefore, if the graph changes, the array may become
           *invalid* and any operation on it will fail with a
           :class:`ValueError` exception. Do **not** store the array if
           the graph is to be modified; store a **copy** instead.
        """
        a = self._get_data()
        if a is None:
            raise ValueError("Cannot get array for value type: " + self.value_type())
        return PropertyArray(a, prop_map=self)

    def _get_data(self):
        g = self.get_graph()
        if g is None:
            raise ValueError("Cannot get array for an orphaned property map")
        if self.__key_type == 'v':
            n = g._Graph__graph.GetNumberOfVertices(False)
        elif self.__key_type == 'e':
            n = g.max_edge_index
        else:
            n = 1
        a = self.__map.get_array(n)
        return a

    def __set_array(self, v):
        a = self.get_array()
        a[:] = v

    a = property(get_array, __set_array,
                 doc=r"""Shortcut to the :meth:`~PropertyMap.get_array` method
                 as an attribute. This makes assignments more convenient, e.g.:

                 >>> g = gt.Graph()
                 >>> g.add_vertex(10)
                 <...>
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
            if g.get_vertex_filter()[0] is not None:
                filt = (g.new_edge_property("bool"), filt[1])
                u = GraphView(g, directed=True, skip_properties=True)
                libcore.mark_edges(u._Graph__graph, _prop("e", u, filt[0]))
                if filt[1]:
                    filt[0].a = 1 - filt[0].a
            elif g._get_max_edge_index() != g.num_edges():
                filt = (g.new_edge_property("bool"), False)
                u = GraphView(g, directed=True, skip_properties=True)
                libcore.mark_edges(u._Graph__graph, _prop("e", u, filt[0]))
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
                try:
                    a[:] = v
                except ValueError:
                    a[:] = v[:len(a)]
            else:
                m = filt[0].a == (not filt[1])
                try:
                    a[m] = v
                except ValueError:
                    a[m] = v[:len(m)][m]

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
            if g.get_vertex_filter()[0] is not None:
                filt = (g.new_edge_property("bool"), filt[1])
                u = GraphView(g, directed=True, skip_properties=True)
                libcore.mark_edges(u._Graph__graph, _prop("e", g, filt[0]))
                if filt[1]:
                    filt[0].a = 1 - filt[0].a
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

    def get_2d_array(self, pos):
        r"""Return a two-dimensional array with a copy of the entries of the
        vector-valued property map. The parameter ``pos`` must be a sequence of
        integers which specifies the indexes of the property values which will
        be used. """

        if self.key_type() == "g":
            raise ValueError("Cannot create multidimensional array for graph property maps.")
        if "vector" not in self.value_type() and (len(pos) > 1 or pos[0] != 0):
            raise ValueError("Cannot create array of dimension %d (indexes %s) from non-vector property map of type '%s'." \
                             % (len(pos), str(pos), self.value_type()))
        if "string" in self.value_type():
            if "vector" in self.value_type():
                p = ungroup_vector_property(self, pos)
            else:
                p = [self]
            g = self.get_graph()
            if self.key_type() == "v":
                N = g.num_vertices()
            else:
                N = g.num_edges()
            a = [["" for j in range(N)] for i in range(len(p))]
            if self.key_type() == "v":
                iters = g.vertices()
            else:
                iters = g.edges()
            for v in iters:
                for i in range(len(p)):
                    a[i][int(v)] = p[i][v]
            if len(a) == 1:
                a = a[0]
            return numpy.array(a)

        try:
            return numpy.array(self.fa)
        except ValueError:
            p = ungroup_vector_property(self, pos)
            return numpy.array([x.a for x in p])

    def set_2d_array(self, a, pos=None):
        r"""Set the entries of the vector-valued property map from a
        two-dimensional array ``a``. If given, the parameter ``pos`` must be a
        sequence of integers which specifies the indexes of the property values
        which will be set."""

        if self.key_type() == "g":
            raise ValueError("Cannot set multidimensional array for graph property maps.")
        if "vector" not in self.value_type():
            if len(a.shape) != 1:
                raise ValueError("Cannot set array of shape %s to non-vector property map of type %s" % \
                                 (str(a.shape), self.value_type()))
            if self.value_type() != "string":
                self.fa = a
            else:
                g = self.get_graph()
                if self.key_type() == "v":
                    iters = g.vertices()
                else:
                    iters = g.edges()
                for i, v in enumerate(iters):
                    self[v] = a[i]
            return

        val = self.value_type()[7:-1]
        ps = []
        for i in range(a.shape[0]):
            ps.append(self.get_graph().new_property(self.key_type(), val))
            if self.value_type() != "string":
                ps[-1].fa = a[i]
            else:
                g = self.get_graph()
                if self.key_type() == "v":
                    iters = g.vertices()
                else:
                    iters = g.edges()
                for j, v in enumerate(iters):
                    ps[-1][v] = a[i, j]
        group_vector_property(ps, val, self, pos)

    def is_writable(self):
        """Return True if the property is writable."""
        return self.__map.is_writable()

    def __call__(self, a):
        p = self.copy()
        p.fa = a
        return p

    def __getstate__(self):
        g = self.get_graph()
        if g is None:
            raise ValueError("cannot pickle orphaned property map")
        value_type = self.value_type()
        key_type = self.key_type()
        if not self.is_writable():
            vals = None
        else:
            u = GraphView(g, skip_vfilt=True, skip_efilt=True)
            if key_type == "v":
                vals = [_convert(self, self[v]) for v in u.vertices()]
            elif key_type == "e":
                vals = [_convert(self, self[e]) for e in u.edges()]
            else:
                vals = _convert(self, self[g])

        state = dict(g=g, value_type=value_type,
                     key_type=key_type, vals=vals,
                     is_vindex=self is g.vertex_index,
                     is_eindex=self is g.edge_index)

        return state

    def __setstate__(self, state):
        g = state["g"]
        key_type = state["key_type"]
        value_type = state["value_type"]
        vals = state["vals"]

        if state["is_vindex"]:
            pmap = g.vertex_index
        elif state["is_eindex"]:
            pmap = g.edge_index
        else:
            u = GraphView(g, skip_vfilt=True, skip_efilt=True)
            if key_type == "v":
                pmap = g.new_vertex_property(value_type)
                for i, v in enumerate(u.vertices()):
                    pmap[v] = vals[i]
            elif key_type == "e":
                pmap = g.new_edge_property(value_type)
                for i, e in enumerate(u.edges()):
                    pmap[e] = vals[i]
            else:
                pmap = g.new_graph_property(value_type)
                pmap[g] = vals

        self.__map = pmap.__map
        self.__g = pmap.__g
        self.__base_g = pmap.__base_g
        self.__key_type = key_type
        self.__register_map()


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
    >>> gt.seed_rng(42)
    >>> g = gt.random_graph(100, lambda: (3, 3))
    >>> props = [g.new_vertex_property("int") for i in range(3)]
    >>> for i in range(3):
    ...    props[i].a = randint(0, 100, g.num_vertices())
    >>> gprop = gt.group_vector_property(props)
    >>> print(gprop[g.vertex(0)].a)
    [51 25  8]
    >>> print(array([p[g.vertex(0)] for p in props]))
    [51 25  8]
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
            u = GraphView(g, directed=True, reversed=g.is_reversed(),
                          skip_properties=True)
            libcore.group_vector_property(u._Graph__graph, _prop(k, g, vprop),
                                          _prop(k, g, p),
                                          i if pos == None else pos[i],
                                          k == 'e')
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
    >>> gt.seed_rng(42)
    >>> g = gt.random_graph(100, lambda: (3, 3))
    >>> prop = g.new_vertex_property("vector<int>")
    >>> for v in g.vertices():
    ...    prop[v] = randint(0, 100, 3)
    >>> uprops = gt.ungroup_vector_property(prop, [0, 1, 2])
    >>> print(prop[g.vertex(0)].a)
    [51 92 14]
    >>> print(array([p[g.vertex(0)] for p in uprops]))
    [51 92 14]
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
            u = GraphView(g, directed=True, reversed=g.is_reversed(),
                          skip_properties=True)
            libcore.ungroup_vector_property(u._Graph__graph,
                                            _prop(k, g, vprop),
                                            _prop(k, g, props[i]),
                                            p, k == 'e')
        else:
            if len(vprop[g]) <= pos[i]:
                vprop[g].resize(pos[i] + 1)
            props[i][g] = vprop[g][pos[i]]
    return props


def infect_vertex_property(g, prop, vals=None):
    """Propagate the `prop` values of vertices with value `val` to all their
    out-neighbours.

    Parameters
    ----------
    prop : :class:`~graph_tool.PropertyMap`
        Property map to be modified.
    vals : list (optional, default: `None`)
        List of values to be propagated. If not provided, all values
        will be propagated.

    Returns
    -------
    None : ``None``

    Examples
    --------
    >>> from numpy.random import seed
    >>> seed(42)
    >>> gt.seed_rng(42)
    >>> g = gt.random_graph(100, lambda: (3, 3))
    >>> prop = g.vertex_index.copy("int32_t")
    >>> gt.infect_vertex_property(g, prop, [10])
    >>> print(sum(prop.a == 10))
    4
    """
    libcore.infect_vertex_property(g._Graph__graph, _prop("v", g, prop),
                                   vals)


def edge_difference(g, prop, ediff=None):
    """Return an edge property map corresponding to the difference between the
    values of `prop` of target and source vertices of each edge.

    Parameters
    ----------
    prop : :class:`~graph_tool.PropertyMap`
        Vertex property map to be used to compute the difference..
    ediff : :class:`~graph_tool.PropertyMap` (optional, default: `None`)
        If provided, the difference values will be stored in this property map.

    Returns
    -------
    ediff : :class:`~graph_tool.PropertyMap`
        Edge differences.

    Examples
    --------
    >>> gt.seed_rng(42)
    >>> g = gt.random_graph(100, lambda: (3, 3))
    >>> ediff = gt.edge_difference(g, g.vertex_index)
    >>> print(ediff.a)
    [ 47  66  10  -1 -65 -56 -91 -62 -38  -1 -86   3 -46 -62 -23 -17 -75 -74
      23 -22   9   2 -35   8 -24 -16 -59 -60  32 -13 -43  30  26 -33  32  -8
      17 -29  -3 -38 -45  41  50  11 -37  58  13 -23  20  48 -17  64  38  22
      63  32  17   7 -10  34  71  20  43  -3  22  13  70  15  63  16  64  86
      37  43  32  -4 -34 -19  36  67  63  65  74  39  17  43  10  29  28  37
      67  93  61  81  27  -2  20  68  50  45  30  71  -1  39  48 -11  25  50
     -25  32  -2  22  -3  40 -29  -7   7  33   8  29 -40  47 -31  -6 -48  -1
     -13   4 -27 -57   9   2  22 -35  23  -2  21  25 -67 -39 -14 -82 -42  -7
     -74 -27 -43 -10 -22 -36  28  21   1 -15 -71 -63 -72 -14 -40 -78 -43 -23
     -83 -12 -43 -26 -48 -52 -58 -82 -14 -44 -62   1 -31 -51 -84 -42 -37 -59
       4 -21 -63 -60 -22 -77 -64 -11   4  18 -56  19   4 -12   3  23 -29 -14
     -15 -38 -20 -34 -62 -49 -47  29  32  -1 -53 -21  19  13  15 -21  12 -30
      -4 -33  66  39 -17   2  50  16  28 -59 -50 -58  69  53  68  -3  60  35
      68 -13  -6  48  71  58  14  58  -4  51  52  13  32  -8  -2  63  56  19
      47  -1  53  39  45 -18 -13  56  52 -16  44  49 -16  35 -11  49  14 -17
      33  11  -8 -31  38 -46  48 -34  47  33 -45 -29  16   3 -48 -44 -58  -9
     -12 -57 -21 -47  34  10 -23 -44 -37 -44   1  11]
    """
    val_t = prop.value_type()
    if val_t == "unsigned long":
        val_t = "int32_t"
    if ediff is None:
        ediff = g.new_edge_property(val_t)
    if ediff.value_type() != val_t:
        raise ValueError("'ediff' must be of the same value type as 'prop': " +
                         val_t)
    if not g.is_directed():
        g = GraphView(g, directed=True, skip_properties=True)
    libcore.edge_difference(g._Graph__graph, _prop("v", g, prop),
                            _prop("e", g, ediff))
    return ediff

@_limit_args({"htype": ["int8_t", "int32_t", "int64_t"]})
def perfect_prop_hash(props, htype="int32_t"):
    """Given a list of property maps `props` of the same type, a derived list of
    property maps with integral type `htype` is retured, where each value is
    replaced by a perfect (i.e. unique) hash value.

    .. note::
       The hash value is deterministic, but it will not be necessarily the same
       for different values of `props`.
    """

    val_types = set([p.value_type() for p in props])
    if len(val_types) > 1:
        raise ValueError("All properties must have the same value type")
    hprops = [p.get_graph().new_property(p.key_type(), htype) for p in props]

    eprops = [p for p in props if p.key_type() == "e"]
    heprops = [p for p in hprops if p.key_type() == "e"]

    vprops = [p for p in props if p.key_type() == "v"]
    hvprops = [p for p in hprops if p.key_type() == "v"]

    hdict = libcore.any()

    for eprop, heprop in zip(eprops, heprops):
        g = eprop.get_graph()
        g = GraphView(g, directed=True, skip_properties=True)
        libcore.perfect_ehash(g._Graph__graph, _prop('e', g, eprop),
                              _prop('e', g, heprop), hdict)

    for vprop, hvprop in zip(vprops, hvprops):
        g = vprop.get_graph()
        g = GraphView(g, directed=True, skip_properties=True)
        libcore.perfect_vhash(g._Graph__graph, _prop('v', g, vprop),
                              _prop('v', g, hvprop), hdict)

    return hprops

class PropertyDict(dict):
    """Wrapper for the dict of vertex, graph or edge properties, which sets the
    value on the property map when changed in the dict.

    .. note::

        The class is only an one-way proxy to the internally-kept properties. If
        you modify this object, the change will be propagated to the internal
        dictionary, but not vice-versa. Keep this in mind if you intend to keep
        a copy of the class instance.
    """
    def __init__(self, g, old, get_func, set_func, del_func):
        dict.__init__(self)
        dict.update(self, old)
        self.g = g
        self.get_func = get_func
        self.set_func = set_func
        self.del_func = del_func

    def __getitem__(self, key):
        if self.get_func != None:
            val = self.get_func(self.g, key)
            dict.__setitem__(self, key, val)
            return val
        else:
            raise KeyError("Property dict cannot be gotten")

    def __setitem__(self, key, val):
        if self.set_func != None:
            self.set_func(self.g, key, val)
        else:
            raise KeyError("Property dict cannot be set")
        dict.__setitem__(self, key, val)

    def __delitem__(self, key):
        self.del_func(self.g, key)
        dict.__delitem__(self, key)

    def clear(self):
        for k in self.keys():
            self.del_func(self.g, k)
        dict.clear(self)

################################################################################
# Graph class
# The main graph interface
################################################################################

from .libgraph_tool_core import Vertex, EdgeBase, Vector_bool, Vector_int16_t, \
    Vector_int32_t, Vector_int64_t, Vector_double, Vector_long_double, \
    Vector_string, new_vertex_property, new_edge_property, new_graph_property


class Graph(object):
    """Generic multigraph class.

    This class encapsulates either a directed multigraph (default or if
    ``directed=True``) or an undirected multigraph (if ``directed=False``),
    with optional internal edge, vertex or graph properties.

    If ``g`` is specified, the graph (and its internal properties) will be
    copied.

    If ``prune`` is set to ``True``, and ``g`` is specified, only the filtered
    graph will be copied, and the new graph object will not be
    filtered. Optionally, a tuple of three booleans can be passed as value to
    ``prune``, to specify a different behavior to vertex, edge, and reversal
    filters, respectively.

    If ``vorder`` is specified, it should correspond to a vertex
    :class:`~graph_tool.PropertyMap` specifying the ordering of the vertices in
    the copied graph.

    The graph is implemented as an `adjacency list`_, where both vertex and edge
    lists are C++ STL vectors.

    .. _adjacency list: http://en.wikipedia.org/wiki/Adjacency_list

    """

    def __init__(self, g=None, directed=True, prune=False, vorder=None):
        self.__properties = {}
        self.__known_properties = {}
        self.__filter_state = {"reversed": False,
                               "edge_filter": (None, False),
                               "vertex_filter": (None, False),
                               "directed": True}
        if g is None:
            self.__graph = libcore.GraphInterface()
            self.set_directed(directed)

            # internal index maps
            self.__vertex_index = \
                     PropertyMap(libcore.get_vertex_index(self.__graph), self, "v")
            self.__edge_index = \
                     PropertyMap(libcore.get_edge_index(self.__graph), self, "e")

        else:
            if isinstance(prune, bool):
                vprune = eprune = rprune = prune
            else:
                vprune, eprune, rprune = prune
            if not (vprune or eprune or rprune):
                gv = GraphView(g, skip_vfilt=True,
                               skip_efilt=True)
                if not rprune:
                    gv.set_reversed(False)
            else:
                gv = g

            # The filters may or may not not be in the internal property maps
            vfilt = g.get_vertex_filter()[0]
            efilt = g.get_edge_filter()[0]

            if (vorder is None and ((vfilt is None and efilt is None) or
                                    (not vprune and not eprune))):
                # Do a simpler, faster copy.
                self.__graph = libcore.GraphInterface(gv.__graph, False,
                                                      [], [], None)

                # internal index maps
                self.__vertex_index = \
                         PropertyMap(libcore.get_vertex_index(self.__graph), self, "v")
                self.__edge_index = \
                         PropertyMap(libcore.get_edge_index(self.__graph), self, "e")

                nvfilt = nefilt = None
                for k, m in g.properties.items():
                    nmap = self.copy_property(m, g=gv)
                    self.properties[k] = nmap
                    if m is vfilt:
                        nvfilt = nmap
                    if m is efilt:
                        nefilt = nmap
                if vfilt is not None:
                    if nvfilt is None:
                        nvfilt = self.copy_property(vfilt, g=gv)
                if efilt is not None:
                    if nefilt is None:
                        nefilt = self.copy_property(efilt, g=gv)
                self.set_filters(nefilt, nvfilt,
                                 inverted_edges=g.get_edge_filter()[1],
                                 inverted_vertices=g.get_vertex_filter()[1])
            else:

                # Copy all internal properties from original graph.
                vprops = []
                eprops = []
                ef_pos = vf_pos = None
                for k, m in gv.vertex_properties.items():
                    if not m.is_writable():
                        m = m.copy("int32_t")
                    if not vprune and m is vfilt:
                        vf_pos = len(vprops)
                    vprops.append([_prop("v", gv, m), libcore.any()])
                for k, m in gv.edge_properties.items():
                    if not m.is_writable():
                        m = m.copy("int32_t")
                    if not eprune and m is efilt:
                        ef_pos = len(eprops)
                    eprops.append([_prop("e", gv, m), libcore.any()])
                if not vprune and vf_pos is None and vfilt is not None:
                    vf_pos = len(vprops)
                    vprops.append([_prop("v", gv, vfilt), libcore.any()])
                if not eprune and ef_pos is None and efilt is not None:
                    ef_pos = len(eprops)
                    eprops.append([_prop("e", gv, efilt), libcore.any()])

                # The vertex ordering
                if vorder is None:
                    vorder = gv.new_vertex_property("int")
                    vorder.fa = numpy.arange(gv.num_vertices())

                # The actual copying of the graph and property maps
                self.__graph = libcore.GraphInterface(gv.__graph, False,
                                                      vprops,
                                                      eprops,
                                                      _prop("v", gv, vorder))
                # internal index maps
                self.__vertex_index = \
                         PropertyMap(libcore.get_vertex_index(self.__graph), self, "v")
                self.__edge_index = \
                         PropertyMap(libcore.get_edge_index(self.__graph), self, "e")

                # Put the copied properties in the internal dictionary
                for i, (k, m) in enumerate(gv.vertex_properties.items()):
                    pmap = new_vertex_property(m.value_type() if m.is_writable() else "int32_t",
                                               self.__graph.GetVertexIndex(),
                                               vprops[i][1])
                    self.vertex_properties[k] = PropertyMap(pmap, self, "v")

                for i, (k, m) in enumerate(gv.edge_properties.items()):
                    pmap = new_edge_property(m.value_type() if m.is_writable() else "int32_t",
                                             self.__graph.GetEdgeIndex(),
                                             eprops[i][1])
                    self.edge_properties[k] = PropertyMap(pmap, self, "e")

                for k, v in gv.graph_properties.items():
                    new_p = self.new_graph_property(v.value_type())
                    new_p[self] = v[gv]
                    self.graph_properties[k] = new_p

                epmap = vpmap = None
                if vf_pos is not None:
                    vpmap = new_vertex_property("bool",
                                                self.__graph.GetVertexIndex(),
                                                vprops[vf_pos][1])
                    vpmap = PropertyMap(vpmap, self, "v")
                if ef_pos is not None:
                    epmap = new_edge_property("bool",
                                              self.__graph.GetEdgeIndex(),
                                              eprops[ef_pos][1])
                    epmap = PropertyMap(epmap, self, "e")
                self.set_filters(epmap, vpmap,
                                 inverted_edges=g.get_edge_filter()[1],
                                 inverted_vertices=g.get_vertex_filter()[1])

            if not rprune:
                self.set_reversed(g.is_reversed())

            # directedness is always a filter
            self.set_directed(g.is_directed())

        # modification permissions
        self.__perms = {"add_edge": True, "del_edge": True,
                        "add_vertex": True, "del_vertex": True}

    def copy(self):
        """Return a deep copy of self. All :ref:`internal property maps <sec_internal_props>`
        are also copied."""
        return Graph(self)

    def __copy__(self):
        return self.copy()

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
        >>> vlist = list(g.add_vertex(5))
        >>> vlist2 = []
        >>> for v in g.vertices():
        ...     vlist2.append(v)
        ...
        >>> assert(vlist == vlist2)

        """
        return libcore.get_vertices(weakref.ref(self))

    def vertex(self, i, use_index=True):
        """Return the vertex with index ``i``. If ``use_index=False``, the
        ``i``-th vertex is returned (which can differ from the vertex with index
        ``i`` in case of filtered graphs). """
        vfilt = self.get_vertex_filter()
        if vfilt[0] is None or not use_index:
            return libcore.get_vertex(weakref.ref(self), int(i))
        try:
            self.set_vertex_filter(None)
            v = libcore.get_vertex(weakref.ref(self), int(i))
        finally:
            self.set_vertex_filter(vfilt[0], vfilt[1])
        if not v.is_valid():
            return v
        if vfilt[0] is not None and vfilt[0][v] == vfilt[1]:
            return None
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
        if s is None or t is None:
            return None
        efilt = self.get_edge_filter()
        edges = []
        for e in s.out_edges():
            if efilt[0] is not None and efilt[0][e] == efilt[1]:
                continue
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
        return libcore.get_edges(weakref.ref(self))

    def add_vertex(self, n=1):
        """Add a vertex to the graph, and return it. If ``n != 1``, ``n``
        vertices are inserted and an iterator over the new vertices is returned.
        This operation is :math:`O(n)`.
        """
        if n == 0:
            return (None for i in range(0, 0))

        self.__check_perms("add_vertex")
        v = libcore.add_vertex(weakref.ref(self), n)

        if n <= 1:
            return v
        else:
            pos = self.num_vertices() - n
            return (self.vertex(i) for i in range(pos, pos + n))

    def remove_vertex(self, vertex, fast=False):
        r"""Remove a vertex from the graph.

        .. note::

           If the option ``fast == False`` is given, this operation is
           :math:`O(N + E)` (this is the default). Otherwise it is
           :math:`O(k + k_{\text{last}})`, where :math:`k` is the (total)
           degree of the vertex being deleted, and :math:`k_{\text{last}}` is
           the (total) degree of the vertex with the largest index.

        .. warning::

           If ``fast == True``, the vertex being deleted is 'swapped' with the
           last vertex (i.e. with the largest index), which will in turn inherit
           the index of the vertex being deleted. All property maps associated
           with the graph will be properly updated, but the index ordering of
           the graph will no longer be the same.

        """
        self.__check_perms("del_vertex")
        vertex = self.vertex(int(vertex))
        index = self.vertex_index[vertex]
        back = self.__graph.GetNumberOfVertices(False) - 1

        # move / shift all known property maps
        if index != back:
            for pmap in self.__known_properties.values():
                if pmap() is not None and pmap().key_type() == "v" and pmap().is_writable():
                    if fast:
                        self.__graph.MoveVertexProperty(pmap()._PropertyMap__map.get_map(), index)
                    else:
                        self.__graph.ShiftVertexProperty(pmap()._PropertyMap__map.get_map(), index)

        libcore.remove_vertex(self.__graph, vertex, fast)

    def clear_vertex(self, vertex):
        """Remove all in and out-edges from the given vertex."""
        del_es = set()
        for e in vertex.all_edges():
            del_es.add(e)
        for e in del_es:
            self.remove_edge(e)

    def add_edge(self, source, target):
        """Add a new edge from ``source`` to ``target`` to the graph, and return
        it. This operation is :math:`O(1)`."""
        self.__check_perms("add_edge")
        e = libcore.add_edge(weakref.ref(self), self.vertex(int(source)),
                             self.vertex(int(target)))
        efilt = self.get_edge_filter()
        if efilt[0] is not None:
            efilt[0][e] = not efilt[1]
        return e

    def remove_edge(self, edge):
        r"""Remove an edge from the graph.

        .. note::

           This operation is normally :math:`O(k_s + k_t)`, where :math:`k_s`
           and :math:`k_s` are the total degrees of the source and target
           vertices, respectively. However, if :meth:`~Graph.set_fast_edge_removal`
           is set to `True`, this operation becomes :math:`O(1)`.

        .. warning::

           The relative ordering of the remaining edges in the graph is kept
           unchanged, unless :meth:`~Graph.set_fast_edge_removal` is set to
           `True`, in which case it can change.
        """
        self.__check_perms("del_edge")
        return libcore.remove_edge(self.__graph, edge)

    def add_edge_list(self, edge_list):
        """Add a list of edges to the graph, given by ``edge_list``, which can
        be a list of ``(source, target)`` pairs where both ``source`` and
        ``target`` are vertex indexes, or a :class:`~numpy.ndarray` of shape
        ``(E,2)``, where ``E`` is the number of edges, and each line specifies a 
        ``(source, target)`` pair. If the list references vertices which do not
        exist in the graph, they will be created."""
        self.__check_perms("add_edge")
        edges = numpy.asarray(edge_list)
        libcore.add_edge_list(self.__graph, edges)

    def set_fast_edge_removal(self, fast=True):
        r"""If ``fast == True`` the fast :math:`O(1)` removal of edges will be
        enabled. This requires an additional data structure of size :math:`O(E)`
        to be kept at all times.  If ``fast == False``, this data structure is
        destroyed."""
        self.__graph.SetKeepEpos(fast)

    def get_fast_edge_removal(self):
        r"""Return whether the fast :math:`O(1)` removal of edges is currently
        enabled."""
        return self.__graph.GetKeepEpos()

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
    def __set_property(self, t, k, v):
        if t == "g" and not isinstance(v, PropertyMap):
            self.__properties[(t, k)][self] = v
        else:
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
        props = dict([(k[1], v) for k, v in self.__properties.items() \
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
    vp = property(__get_vertex_properties,
                  doc="Alias to :attr:`~Graph.vertex_properties`.")

    # edge properties
    def __get_edge_properties(self):
        return PropertyDict(self, self.__get_specific_properties("e"),
                            lambda g, k: g.__properties[("e", k)],
                            lambda g, k, v: g.__set_property("e", k, v),
                            lambda g, k: g.__del_property("e", k))
    edge_properties = property(__get_edge_properties,
                                 doc="Dictionary of internal edge properties. The keys are the property names.")
    ep = property(__get_edge_properties,
                  doc="Alias to :attr:`~Graph.edge_properties`.")

    # graph properties
    def __get_graph_properties(self):
        return PropertyDict(self, self.__get_specific_properties("g"),
                            lambda g, k: g.__properties[("g", k)][g],
                            lambda g, k, v: g.__set_property("g", k, v),
                            lambda g, k: g.__del_property("g", k))
    graph_properties = property(__get_graph_properties,
                                 doc="Dictionary of internal graph properties. The keys are the property names.")
    gp = property(__get_graph_properties,
                  doc="Alias to :attr:`~Graph.graph_properties`.")

    def own_property(self, prop):
        """Return a version of the property map 'prop' (possibly belonging to
        another graph) which is owned by the current graph."""
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
        w = max([len(x[0]) for x in list(self.__properties.keys())]) + 4
        w = w if w > 14 else 14

        for k, v in self.__properties.items():
            if k[0] == "g":
                pref="%%-%ds (graph)   (type: %%s, val: " % w % \
                      (k[1], v.value_type())
                val = str(v[self])
                if len(val) > 1000:
                    val = val[:1000] + "..."
                tw = terminal_size()[0]
                val = textwrap.indent(textwrap.fill(val,
                                                    width=max(tw - len(pref), 1)),
                                      " " * len(pref))
                val = val[len(pref):]
                print("%s%s)" % (pref, val))
        for k, v in self.__properties.items():
            if k[0] == "v":
                print("%%-%ds (vertex)  (type: %%s)" % w % (k[1],
                                                            v.value_type()))
        for k, v in self.__properties.items():
            if k[0] == "e":
                print("%%-%ds (edge)    (type: %%s)" % w % (k[1],
                                                            v.value_type()))

    # index properties

    def _get_vertex_index(self):
        return self.__vertex_index
    vertex_index = property(_get_vertex_index,
                            doc="""Vertex index map.

                            It maps for each vertex in the graph an unique
                            integer in the range [0, :meth:`~graph_tool.Graph.num_vertices` - 1].

                            .. note::

                                Like :attr:`~graph_tool.Graph.edge_index`, this
                                is a special instance of a :class:`~graph_tool.PropertyMap`
                                class, which is **immutable**, and cannot be
                                accessed as an array.""")

    def _get_edge_index(self):
        return self.__edge_index
    edge_index = property(_get_edge_index, doc="""Edge index map.

                            It maps for each edge in the graph an unique
                            integer.

                            .. note::

                                Like :attr:`~graph_tool.Graph.vertex_index`, this
                                is a special instance of a :class:`~graph_tool.PropertyMap`
                                class, which is **immutable**, and cannot be
                                accessed as an array.

                                Additionally, the indexes may not necessarily
                                lie in the range [0, :meth:`~graph_tool.Graph.num_edges` - 1].
                                However this will always happen whenever no
                                edges are deleted from the graph.""")

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

    def new_property(self, key_type, value_type, vals=None):
        """Create a new (uninitialized) vertex property map of key type
        ``key_type`` (``v``, ``e`` or ``g``), value type ``value_type``, and
        return it. If provided, the values will be initialized by ``vals``,
        which should be a sequence.
        """
        if key_type == "v" or key_type == "vertex":
            return self.new_vertex_property(value_type, vals)
        if key_type == "e" or key_type == "edge":
            return self.new_edge_property(value_type, vals)
        if key_type == "g" or key_type == "graph":
            return self.new_graph_property(value_type, vals)
        raise ValueError("unknown key type: " + key_type)

    def new_vertex_property(self, value_type, vals=None):
        """Create a new (uninitialized) vertex property map of type ``value_type``,
        and return it. If provided, the values will be initialized by ``vals``,
        which should be either a sequence or a single value."""
        prop = PropertyMap(new_vertex_property(_type_alias(value_type),
                                               self.__graph.GetVertexIndex(),
                                               libcore.any()),
                           self, "v")
        if vals is not None:
            try:
                prop.a = vals
            except ValueError:
                for v, x in zip(self.vertices(), vals):
                    prop[v] = x
        return prop

    def new_edge_property(self, value_type, vals=None):
        """Create a new (uninitialized) edge property map of type
        ``value_type``, and return it. If provided, the values will be
        initialized by ``vals``, which should be a sequence or a single value."""
        prop = PropertyMap(new_edge_property(_type_alias(value_type),
                                             self.__graph.GetEdgeIndex(),
                                             libcore.any()),
                           self, "e")
        if vals is not None:
            try:
                prop.a = vals
            except ValueError:
                for e, x in zip(self.edges(), vals):
                    prop[e] = x
        return prop

    def new_graph_property(self, value_type, val=None):
        """Create a new graph property map of type ``value_type``, and return
        it. If ``val`` is not None, the property is initialized to its value."""
        prop = PropertyMap(new_graph_property(_type_alias(value_type),
                                              self.__graph.GetGraphIndex(),
                                              libcore.any()),
                           self, "g")
        if val is not None:
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

        u = self if g is None else g
        g = GraphView(u, directed=True, reversed=u.is_reversed(),
                      skip_properties=True)

        if src.key_type() == "v":
            self.__graph.CopyVertexProperty(g.__graph, _prop("v", g, src),
                                            _prop("v", self, tgt))
        elif src.key_type() == "e":
            self.__graph.CopyEdgeProperty(g.__graph, _prop("e", g, src),
                                            _prop("e", self, tgt))
        else:
            tgt[self] = src[g]
        return ret

    # degree property map
    @_limit_args({"deg": ["in", "out", "total"]})
    def degree_property_map(self, deg, weight=None):
        """Create and return a vertex property map containing the degree type
        given by ``deg``, which can be any of ``"in"``, ``"out"``, or ``"total"``.
        If provided, ``weight`` should be an edge :class:`~graph_tool.PropertyMap`
        containing the edge weights which should be summed."""
        pmap = self.__graph.DegreeMap(deg, _prop("e", self, weight))
        return PropertyMap(pmap, self, "v")

    # I/O operations
    # ==============
    def __get_file_format(self, file_name):
        fmt = None
        for f in ["gt", "graphml", "xml", "dot", "gml"]:
            names = ["." + f, ".%s.gz" % f, ".%s.bz2" % f, ".%s.xz" % f]
            for name in names:
                if file_name.endswith(name):
                    fmt = f
                    break
        if fmt is None:
            raise ValueError("cannot determine file format of: " + file_name)
        return fmt

    def load(self, file_name, fmt="auto", ignore_vp=None, ignore_ep=None,
             ignore_gp=None):
        """Load graph from ``file_name`` (which can be either a string or a file-like
        object). The format is guessed from ``file_name``, or can be specified
        by ``fmt``, which can be either "gt", "graphml", "xml", "dot" or "gml".
        (Note that "graphml" and "xml" are synonyms).

        If provided, the parameters ``ignore_vp``, ``ignore_ep`` and
        ``ignore_gp``, should contain a list of property names (vertex, edge or
        graph, respectively) which should be ignored when reading the file.

        .. warning::

           The only file formats which are capable of perfectly preserving the
           internal property maps are "gt" and "graphml". Because of this,
           they should be preferred over the other formats whenever possible.

        """

        if isinstance(file_name, str):
            file_name = os.path.expanduser(file_name)
            f = open(file_name) # throw the appropriate exception, if not found
        if fmt == 'auto' and isinstance(file_name, str):
            fmt = self.__get_file_format(file_name)
        elif fmt == "auto":
            fmt = "gt"
        if isinstance(file_name, str) and file_name.endswith(".xz"):
            try:
                import lzma
                file_name = lzma.open(file_name, mode="rb")
            except ImportError:
                raise ValueError("lzma compression is only available in Python >= 3.3")
        if fmt == "graphml":
            fmt = "xml"
        if ignore_vp is None:
            ignore_vp = []
        if ignore_ep is None:
            ignore_ep = []
        if ignore_gp is None:
            ignore_gp = []
        if isinstance(file_name, str):
            props = self.__graph.ReadFromFile(file_name, None, fmt, ignore_vp,
                                              ignore_ep, ignore_gp)
        else:
            props = self.__graph.ReadFromFile("", file_name, fmt, ignore_vp,
                                              ignore_ep, ignore_gp)
        for name, prop in props[0].items():
            self.vertex_properties[name] = PropertyMap(prop, self, "v")
        for name, prop in props[1].items():
            self.edge_properties[name] = PropertyMap(prop, self, "e")
        for name, prop in props[2].items():
            self.graph_properties[name] = PropertyMap(prop, self, "g")
        if "_Graph__save__vfilter" in self.graph_properties:
            self.set_vertex_filter(self.vertex_properties["_Graph__save__vfilter"],
                                   self.graph_properties["_Graph__save__vfilter"])
            del self.vertex_properties["_Graph__save__vfilter"]
            del self.graph_properties["_Graph__save__vfilter"]
        if "_Graph__save__efilter" in self.graph_properties:
            self.set_edge_filter(self.edge_properties["_Graph__save__efilter"],
                                 self.graph_properties["_Graph__save__efilter"])
            del self.edge_properties["_Graph__save__efilter"]
            del self.graph_properties["_Graph__save__efilter"]
        if "_Graph__reversed" in self.graph_properties:
            self.set_reversed(True)
            del self.graph_properties["_Graph__reversed"]



    def save(self, file_name, fmt="auto"):
        """Save graph to ``file_name`` (which can be either a string or a file-like
        object). The format is guessed from the ``file_name``, or can be
        specified by ``fmt``, which can be either "gt", "graphml", "xml", "dot"
        or "gml".  (Note that "graphml" and "xml" are synonyms).

        .. warning::

           The only file formats which are capable of perfectly preserving the
           internal property maps are "gt" and "graphml". Because of this,
           they should be preferred over the other formats whenever possible.

        """

        u = GraphView(self, reversed=self.is_reversed(), skip_vfilt=True,
                      skip_efilt=True)

        if self.get_vertex_filter()[0] is not None:
            u.graph_properties["_Graph__save__vfilter"] = self.new_graph_property("bool")
            u.vertex_properties["_Graph__save__vfilter"] =  self.get_vertex_filter()[0]
            u.graph_properties["_Graph__save__vfilter"] = self.get_vertex_filter()[1]
        if self.get_edge_filter()[0] is not None:
            u.graph_properties["_Graph__save__efilter"] = self.new_graph_property("bool")
            u.edge_properties["_Graph__save__efilter"] = self.get_edge_filter()[0]
            u.graph_properties["_Graph__save__efilter"] = self.get_edge_filter()[1]

        if self.is_reversed():
            u.graph_properties["_Graph__reversed"] = self.new_graph_property("bool")
            u.graph_properties["_Graph__reversed"] = True

        if type(file_name) == str:
            file_name = os.path.expanduser(file_name)
        if fmt == 'auto' and isinstance(file_name, str):
            fmt = self.__get_file_format(file_name)
        elif fmt == "auto":
            fmt = "gt"
        if fmt == "graphml":
            fmt = "xml"

        if isinstance(file_name, str) and file_name.endswith(".xz"):
            try:
                import lzma
                file_name = lzma.open(file_name, mode="wb")
            except ImportError:
                raise ValueError("lzma compression is only available in Python >= 3.3")

        props = [(name[1], prop._PropertyMap__map) for name, prop in \
                 self.__properties.items()]

        if isinstance(file_name, str):
            f = open(file_name, "w") # throw the appropriate exception, if
                                     # unable to open
            f.close()
            u.__graph.WriteToFile(file_name, None, fmt, props)
        else:
            u.__graph.WriteToFile("", file_name, fmt, props)


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

    def set_filters(self, eprop, vprop, inverted_edges=False, inverted_vertices=False):
        """Set the boolean properties for edge and vertex filters, respectively.
        Only the vertices and edges with value different than ``True`` are kept in
        the filtered graph. If either the ``inverted_edges`` or ``inverted_vertex``
        options are supplied with the value ``True``, only the edges or vertices
        with value ``False`` are kept. If any of the supplied property is ``None``,
        an empty filter is constructed which allows all edges or vertices."""

        if eprop is None and vprop is None:
            return

        if eprop is None:
            eprop = self.new_edge_property("bool")
            eprop.a = not inverted_edges

        if vprop is None:
            vprop = self.new_vertex_property("bool")
            vprop.a = not inverted_vertices

        self.__graph.SetVertexFilterProperty(_prop("v", self, vprop),
                                             inverted_vertices)
        self.__filter_state["vertex_filter"] = (vprop, inverted_vertices)

        self.__graph.SetEdgeFilterProperty(_prop("e", self, eprop),
                                           inverted_edges)
        self.__filter_state["edge_filter"] = (eprop, inverted_edges)

    def set_vertex_filter(self, prop, inverted=False):
        """Set the vertex boolean filter property. Only the vertices with value
        different than ``False`` are kept in the filtered graph. If the ``inverted``
        option is supplied with value ``True``, only the vertices with value
        ``False`` are kept. If the supplied property is ``None``, the filter is
        replaced by an uniform filter allowing all vertices."""

        if prop is not None and prop.value_type() != "bool":
            raise ValueError("filter property map must have 'bool' type")

        vfilt = prop
        efilt = None

        eprop = self.get_edge_filter()
        if eprop[0] is None and vfilt is not None:
            efilt = self.new_edge_property("bool")
            efilt.a = True
        if eprop[0] is not None and vfilt is None:
            vfilt = self.new_vertex_property("bool")
            vfilt.a = not inverted

        self.__graph.SetVertexFilterProperty(_prop("v", self, vfilt),
                                             inverted)
        self.__filter_state["vertex_filter"] = (vfilt, inverted)

        if efilt is not None:
            self.set_edge_filter(efilt)

    def get_vertex_filter(self):
        """Return a tuple with the vertex filter property and bool value
        indicating whether or not it is inverted."""
        return self.__filter_state["vertex_filter"]

    def set_edge_filter(self, prop, inverted=False):
        """Set the edge boolean filter property. Only the edges with value
        different than ``False`` are kept in the filtered graph. If the ``inverted``
        option is supplied with value ``True``, only the edges with value ``False``
        are kept. If the supplied property is ``None``, the filter is
        replaced by an uniform filter allowing all edges."""

        if prop is not None and prop.value_type() != "bool":
            raise ValueError("filter property map must have 'bool' type")

        efilt = prop
        vfilt = None

        vprop = self.get_vertex_filter()
        if vprop[0] is None and efilt is not None:
            vfilt = self.new_vertex_property("bool")
            vfilt.a = True
        if vprop[0] is not None and efilt is None:
            efilt = self.new_edge_property("bool")
            efilt.a = not inverted

        self.__graph.SetEdgeFilterProperty(_prop("e", self, efilt), inverted)
        self.__filter_state["edge_filter"] = (efilt, inverted)

        if vfilt is not None:
            self.set_vertex_filter(vfilt)

    def get_edge_filter(self):
        """Return a tuple with the edge filter property and bool value
        indicating whether or not it is inverted."""
        return self.__filter_state["edge_filter"]

    def clear_filters(self):
        """Remove vertex and edge filters, and set the graph to the unfiltered
        state."""
        self.__graph.SetVertexFilterProperty(_prop("v", self, None), False)
        self.__filter_state["vertex_filter"] = (None, False)
        self.__graph.SetEdgeFilterProperty(_prop("e", self, None), False)
        self.__filter_state["edge_filter"] = (None, False)

    def purge_vertices(self, in_place=False):
        """Remove all vertices of the graph which are currently being filtered
        out, and return it to the unfiltered state. This operation is not
        reversible.

        If the option ``in_place == True`` is given, the algorithm will remove
        the filtered vertices and re-index all property maps which are tied with
        the graph. This is a slow operation which has an :math:`O(N^2)`
        complexity.

        If ``in_place == False``, the graph and its vertex and edge property
        maps are temporarily copied to a new unfiltered graph, which will
        replace the contents of the original graph. This is a fast operation
        with an :math:`O(N + E)` complexity. This is the default behaviour if no
        option is given.

        """
        if in_place:
            old_indexes = self.vertex_index.copy("int64_t")
            self.__graph.PurgeVertices(_prop("v", self, old_indexes))
            self.set_vertex_filter(None)
            for pmap in self.__known_properties.values():
                if (pmap() is not None and pmap().key_type() == "v" and
                    pmap().is_writable() and
                    pmap() not in [self.vertex_index, self.edge_index]):
                    self.__graph.ReIndexVertexProperty(pmap()._PropertyMap__map.get_map(),
                                                       _prop("v", self, old_indexes))
        else:
            stamp = id(self)
            pmaps = []
            for pmap in self.__known_properties.values():
                if (pmap() is not None and pmap().key_type() in ["v", "e"] and
                    pmap() not in [self.vertex_index, self.edge_index]):
                    pmaps.append(pmap())
                    pname = "__tmp_purge_vertices_%d_%d" % (stamp, id(pmaps[-1]))
                    self.properties[(pmaps[-1].key_type(), pname)] = pmaps[-1]

            new_g = Graph(self, prune=(True, False, False))
            self.__graph = new_g.__graph
            self.set_vertex_filter(None)

            for pmap in pmaps:
                pname = "__tmp_purge_vertices_%d_%d" % (stamp, id(pmap))
                new_pmap = new_g.properties[(pmap.key_type(), pname)]
                pmap._PropertyMap__map = new_pmap._PropertyMap__map
                del self.properties[(pmap.key_type(), pname)]

            # update edge filter if set
            efilt = self.get_edge_filter()
            if efilt[0] is not None:
                self.set_edge_filter(efilt[0], efilt[1])

    def purge_edges(self):
        """Remove all edges of the graph which are currently being filtered out,
        and return it to the unfiltered state. This operation is not reversible."""
        self.__graph.PurgeEdges()
        self.set_edge_filter(None)

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
        return self.__graph.GetNumberOfVertices(True)

    def num_edges(self):
        """Get the number of edges.

        .. note::

            If the edges are being filtered, this operation is
            :math:`O(E)`. Otherwise it is :math:`O(1)`.
        """
        return self.__graph.GetNumberOfEdges(True)

    # Pickling support
    # ================

    def __getstate__(self):
        state = dict()
        if sys.version_info < (3,):
            sio = StringIO.StringIO()
        else:
            sio = io.BytesIO()
        stream = gzip.GzipFile(fileobj=sio, mode="wb")
        self.save(stream, "gt")
        stream.close()
        state["blob"] = sio.getvalue()
        return state

    def __setstate__(self, state):
        self.__init__()
        blob = state["blob"]
        if blob != "":
            try:
                if sys.version_info < (3,):
                    sio = StringIO.StringIO(blob)
                else:
                    sio = io.BytesIO(blob)
                stream = gzip.GzipFile(fileobj=sio, mode="rb")
                self.load(stream, "gt")
            except OSError:
                if sys.version_info < (3,):
                    sio = StringIO.StringIO(blob)
                else:
                    sio = io.BytesIO(blob)
                stream = gzip.GzipFile(fileobj=sio, mode="rb")
                self.load(stream, "xml")


def load_graph(file_name, fmt="auto", ignore_vp=None, ignore_ep=None,
               ignore_gp=None):
    """Load a graph from ``file_name`` (which can be either a string or a file-like object).

    The format is guessed from ``file_name``, or can be specified by ``fmt``,
    which can be either "gt", "graphml", "xml", "dot" or "gml".  (Note that
    "graphml" and "xml" are synonyms).

    If provided, the parameters ``ignore_vp``, ``ignore_ep`` and
    ``ignore_gp``, should contain a list of property names (vertex, edge or
    graph, respectively) which should be ignored when reading the file.

    .. warning::

       The only file formats which are capable of perfectly preserving the
       internal property maps are "gt" and "graphml". Because of this,
       they should be preferred over the other formats whenever possible.

    """
    g = Graph()
    g.load(file_name, fmt, ignore_vp, ignore_ep, ignore_gp)
    return g


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

    If either ``skip_properties``, ``skip_vfilt`` or ``skip_efilt`` is ``True``,
    then the internal properties, vertex filter or edge filter of the original
    graph are ignored, respectively.

    """

    def __init__(self, g, vfilt=None, efilt=None, directed=None,
                 reversed=False, skip_properties=False, skip_vfilt=False,
                 skip_efilt=False):
        self.__base = g if not isinstance(g, GraphView) else g.base
        Graph.__init__(self)
        # copy graph reference
        self._Graph__graph = libcore.GraphInterface(g._Graph__graph, True,
                                                    [], [],
                                                    _prop("v", g, g.vertex_index))

        if not skip_properties:
            for k, v in g.properties.items():
                self.properties[k] = self.own_property(v)

        # set already existent filters
        if not skip_vfilt:
            vf = g.get_vertex_filter()
            if vf[0] is not None:
                self.set_vertex_filter(vf[0], vf[1])
        if not skip_efilt:
            ef = g.get_edge_filter()
            if ef[0] is not None:
                self.set_edge_filter(ef[0], ef[1])

        if vfilt is not None:
            if type(vfilt) is PropertyMap:
                self.set_vertex_filter(vfilt)
            else:
                vmap = self.new_vertex_property("bool")
                if issubclass(type(vfilt), numpy.ndarray):
                    vmap.fa = vfilt
                else:
                    for v in g.vertices():
                        vmap[v] = vfilt(v)
                self.set_vertex_filter(vmap)

        if efilt is not None:
            if type(efilt) is PropertyMap:
                self.set_edge_filter(efilt)
            else:
                emap = self.new_edge_property("bool")
                if issubclass(type(efilt), numpy.ndarray):
                    emap.fa = efilt
                else:
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

    # pickling support
    def __getstate__(self):
        return Graph.__getstate__(self)

    def __setstate__(self, state):
        g = Graph()
        g.__setstate__(state)
        self.__init__(g)


def value_types():
    """Return a list of possible properties value types."""
    return libcore.get_property_types()

# Vertex and Edge Types
# =====================
from .libgraph_tool_core import Vertex, Edge, EdgeBase

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

def _in_degree(self, weight=None):
    """Return the in-degree of the vertex. If provided, ``weight`` should be a
    scalar edge property map, and the in-degree will correspond to the sum of
    the weights of the in-edges.
    """

    if weight is None:
        return self.__in_degree()
    else:
        return self.__weighted_in_degree(_prop("e", self.get_graph(), weight))

Vertex.in_degree = _in_degree

def _out_degree(self, weight=None):
    """Return the out-degree of the vertex. If provided, ``weight`` should be a
    scalar edge property map, and the out-degree will correspond to the sum of
    the weights of the out-edges.
    """

    if weight is None:
        return self.__out_degree()
    else:
        return self.__weighted_out_degree(_prop("e", self.get_graph(), weight))

Vertex.out_degree = _out_degree


def _vertex_repr(self):
    if not self.is_valid():
        return "<invalid Vertex object at 0x%x>" % (id(self))
    return "<Vertex object with index '%d' at 0x%x>" % (int(self), id(self))
Vertex.__repr__ = _vertex_repr

Vertex.__eq__ = lambda v1, v2 : int(v1) == (int(v2) if isinstance(v2, Vertex) else v2)
Vertex.__ne__ = lambda v1, v2 : int(v1) != (int(v2) if isinstance(v2, Vertex) else v2)
Vertex.__lt__ = lambda v1, v2 : int(v1) < (int(v2) if isinstance(v2, Vertex) else v2)
Vertex.__gt__ = lambda v1, v2 : int(v1) > (int(v2) if isinstance(v2, Vertex) else v2)
Vertex.__le__ = lambda v1, v2 : int(v1) <= (int(v2) if isinstance(v2, Vertex) else v2)
Vertex.__ge__ = lambda v1, v2 : int(v1) >= (int(v2) if isinstance(v2, Vertex) else v2)

_edge_doc = """Edge descriptor.

This class represents an edge in a :class:`~graph_tool.Graph`.

:class:`~graph_tool.Edge` instances are hashable, and are convertible to a
tuple, which contains the source and target vertices.
"""

def _edge_cmp(e1, e2):
    te1, te2 = tuple(e1), tuple(e2)
    g1 = e1.get_graph()
    g2 = e2.get_graph()
    if not g1.is_directed():
        te1 = sorted(te1)
    if not g2.is_directed():
        te2 = sorted(te2)
    te1 = (te1, g1.edge_index[e1])
    te2 = (te2, g2.edge_index[e2])
    if te1 < te2:
        return -1
    if te1 > te2:
        return 1
    return 0

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
                    e = next(g.edges())
                    e.__class__.__repr__ = _edge_repr
                    e.__class__.__iter__ = _edge_iter
                    e.__class__.__doc__ = _edge_doc

                    e.__class__.__eq__ = lambda e1, e2 : _edge_cmp(e1, e2) == 0
                    e.__class__.__ne__ = lambda e1, e2 : _edge_cmp(e1, e2) != 0
                    e.__class__.__lt__ = lambda e1, e2 : _edge_cmp(e1, e2) < 0
                    e.__class__.__gt__ = lambda e1, e2 : _edge_cmp(e1, e2) > 0
                    e.__class__.__le__ = lambda e1, e2 : _edge_cmp(e1, e2) <= 0
                    e.__class__.__ge__ = lambda e1, e2 : _edge_cmp(e1, e2) >= 0

init_edge_classes()

# some shenanigans to make it seem there is only a single edge class
EdgeBase.__doc__ = Edge.__doc__
EdgeBase.source = Edge.source
EdgeBase.target = Edge.target
EdgeBase.is_valid = Edge.is_valid
EdgeBase.get_graph = Edge.get_graph
Edge = EdgeBase
Edge.__name__ = "Edge"


# Add convenience function to vector classes
def _get_array_view(self):
    return self.get_array()[:]


def _set_array_view(self, v):
    self.get_array()[:] = v

vector_types = [Vector_bool, Vector_int16_t, Vector_int32_t, Vector_int64_t,
                Vector_double, Vector_long_double]
for vt in vector_types:
    vt.a = property(_get_array_view, _set_array_view,
                    doc=r"""Shortcut to the `get_array` method as an attribute.""")
    vt.__repr__ = lambda self: self.a.__repr__()
Vector_string.a = None
Vector_string.get_array = lambda self: None
Vector_string.__repr__ = lambda self: repr(list(self))


# Global RNG

_rng = libcore.get_rng((numpy.random.randint(0, sys.maxsize) + os.getpid()) % sys.maxsize)

def seed_rng(seed):
    "Seed the random number generator used by graph-tool's algorithms."
    import graph_tool
    graph_tool._rng = libcore.get_rng(int(seed))

def _get_rng():
    global _rng
    return _rng

# OpenMP Setup

def openmp_enabled():
    """Return `True` if OpenMP was enabled during compilation."""
    return libcore.openmp_enabled()

def openmp_get_num_threads():
    """Return the number of OpenMP threads."""
    return libcore.openmp_get_num_threads()

def openmp_set_num_threads(n):
    """Set the number of OpenMP threads."""
    return libcore.openmp_set_num_threads(n)

def openmp_get_schedule():
    """Return the runtime OpenMP schedule and chunk size. The schedule can by
    any of: `"static"`, `"dynamic"`, `"guided"`, `"auto"`."""
    return libcore.openmp_get_schedule()

def openmp_set_schedule(schedule, chunk=0):
    """Set the runtime OpenMP schedule and chunk size. The schedule can by
    any of: `"static"`, `"dynamic"`, `"guided"`, `"auto"`."""
    return libcore.openmp_set_schedule(schedule, chunk)

if openmp_enabled() and os.environ.get("OMP_SCHEDULE") is None:
    openmp_set_schedule("static", 0)
