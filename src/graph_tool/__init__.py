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

"""
``graph_tool`` - a general graph manipulation python module
===========================================================

Provides
   1. A Graph object for graph representation and manipulation
   2. Property maps for Vertex, Edge or Graph.
   3. Fast algorithms implemented in C++.

How to use the documentation
----------------------------

Documentation is available in two forms: docstrings provided
with the code, and the full documentation available in
`the graph-tool homepage <http://graph-tool.forked.de>`_.

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

Summary
-------

.. autosummary::
   :nosignatures:

   Graph
   Vertex
   Edge
   PropertyMap
   load_graph
   group_vector_property
   ungroup_vector_property
   value_types
   show_config

Classes
-------
"""

__author__ = "Tiago de Paula Peixoto <tiago@forked.de>"
__copyright__ = "Copyright 2007-2010 Tiago de Paula Peixoto"
__license__ = "GPL version 3 or above"
__URL__ = "http://graph-tool.forked.de"

# import numpy and scipy before everything to avoid weird segmentation faults
# depending on the order things are imported.

import numpy
import scipy
import scipy.stats

from . core import  __version__, Graph, Vector_bool, Vector_int32_t, \
     Vector_int64_t, Vector_double, Vector_long_double, Vector_string, \
     value_types, load_graph, PropertyMap, Vertex, Edge, \
     group_vector_property, ungroup_vector_property, show_config

__all__ = ["Graph", "Vertex", "Edge", "Vector_bool", "Vector_int32_t",
           "Vector_int64_t", "Vector_double", "Vector_long_double",
           "Vector_string", "value_types", "load_graph", "PropertyMap",
           "group_vector_property", "ungroup_vector_property", "show_config",
           "__author__", "__copyright__", "__URL__", "__version__"]
