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

"""
graph_tool - a general graph manipulation python module
----------

(include general description here)

Documentation is also available in the docstrings.

This module defines:
    Graph: A general graph class
    GraphError: Exceptions raised by the Graph class

Modules include:
    test: Unit testing routines.
"""

__author__="Tiago de Paula Peixoto <tiago@forked.de>"
__copyright__="Copyright 2008 Tiago de Paula Peixoto"
__license__="GPL version 3 or above"
__URL__="http://graph-tool.forked.de"

from . core import  __version__, Graph, GraphError, Vector_bool, \
     Vector_int32_t, Vector_int64_t, Vector_double, Vector_long_double

__all__ = ["Graph", "GraphError", "Vector_bool", "Vector_int32_t",
           "Vector_int64_t", "Vector_double", "Vector_long_double"]
