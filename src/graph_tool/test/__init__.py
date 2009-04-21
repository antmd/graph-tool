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
Unit testing routines for graph_tool.
"""

import unittest, basic, properties, io

def run(verbosity=2):
    """Run all tests with desired verbosity."""
    unittest.TextTestRunner(verbosity=verbosity).run(basic.suite())
    unittest.TextTestRunner(verbosity=verbosity).run(properties.suite())
    unittest.TextTestRunner(verbosity=verbosity).run(io.suite())
