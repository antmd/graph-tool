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
Utility module which includes all the sub-modules in graph_tool
"""

from __future__ import division, absolute_import, print_function

from graph_tool import *
import graph_tool
from graph_tool.correlations import *
import graph_tool.correlations
from graph_tool.centrality import *
import graph_tool.centrality
try:
    from graph_tool.draw import *
    import graph_tool.draw
except ImportError:
    # Proceed despite errors with cairo, matplotlib, etc.
    pass
from graph_tool.stats import *
import graph_tool.stats
from graph_tool.generation import *
import graph_tool.generation
from graph_tool.stats import *
import graph_tool.stats
from graph_tool.clustering import *
import graph_tool.clustering
from graph_tool.community import *
import graph_tool.community
from graph_tool.run_action import *
import graph_tool.run_action
from graph_tool.topology import *
import graph_tool.topology
from graph_tool.flow import *
import graph_tool.flow
from graph_tool.spectral import *
import graph_tool.spectral
from graph_tool.search import *
import graph_tool.search
from graph_tool.util import *
import graph_tool.util
import graph_tool.collection
import graph_tool.collection as collection
