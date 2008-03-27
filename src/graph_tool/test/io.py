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
Test graph IO operations
"""

from graph_tool import *
import unittest
from StringIO import StringIO
from basic import gen_suite
from properties import value_types, type_init

class TestGraphIO(unittest.TestCase):
    """Test graph I/O manipulations"""

    def setUp(self):
        self.g = Graph()

    def add_vertices(self):
        """add 10 vertices"""
        for i in range(0,10):
            v = self.g.add_vertex()
            self.assert_(str(v) == str(i))

    def clear_vertices(self):
        self.g.clear()

    def add_edges(self):
        """add 10x10 edges to form a complete graph with self-loops"""
        for i in self.g.vertices():
            for j in self.g.vertices():
                self.g.add_edge(i,j)

    def add_properties(self):
        """add a bunch of properties"""
        for t in value_types:
            # vertex properties
            prop = self.g.add_vertex_property("prop_" + t, t)
            for v in self.g.vertices():
                prop[v] = type_init[t]

            # edge properties
            prop = self.g.add_edge_property("prop_" + t, t)
            for e in self.g.edges():
                prop[e] = type_init[t]

            # graph properties
            self.g.add_graph_property("prop_" + t, t)
            self.g.graph_properties["prop_" + t] = type_init[t]

    def test_io(self):
        """Testing graph xml io"""
        self.add_vertices()
        self.add_edges()
        self.add_properties()

        stream = StringIO()
        self.g.save(stream, "xml")
        stream.seek(0)
        b = Graph()
        b.load(stream, "xml")
        self.graph_compare(self.g,b)

        self.clear_vertices()

    def graph_compare(self, g1, g2):
        """Compare two graphs"""
        self.assert_(g1.num_vertices() == g2.num_vertices())
        self.assert_(g1.num_edges() == g2.num_edges())
        for p in g1.graph_properties.keys():
            self.assert_(p in g2.graph_properties)
            self.assert_(g1.graph_properties[p] == g2.graph_properties[p])
        for p in g1.vertex_properties.keys():
            self.assert_(p in g2.vertex_properties)
            for v in g1.vertices():
                self.assert_(g1.vertex_properties[p][v] == \
                             g2.vertex_properties[p][v])
        for p in g1.edge_properties.keys():
            self.assert_(p in g2.edge_properties)
            for v in g1.edges():
                self.assert_(g1.edge_properties[p][v] ==\
                             g2.edge_properties[p][v])

def suite():
    return gen_suite(TestGraphIO)
