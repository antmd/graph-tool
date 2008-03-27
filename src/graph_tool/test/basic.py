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
Basic graph manipulation tests
"""

from graph_tool import *
import unittest
from numpy import *

class TestBasicManipulation(unittest.TestCase):
    """Test basic graph manipulation"""

    def setUp(self):
        self.g = Graph()

    def add_vertices(self):
        # add 10 vertices
        for i in range(0,10):
            v = self.g.add_vertex()
            self.assert_(str(v) == str(i))

    def clear_vertices(self):
        self.g.clear()

    def test_add_vertices(self):
        """Adding vertices"""
        self.add_vertices()
        self.assertEqual(self.g.num_vertices(), 10)
        self.clear_vertices()

    def test_iterate_vertices(self):
        """Iterating through vertices"""
        self.add_vertices()
        vset = set()
        for v in self.g.vertices():
            vset.add(v)
        self.assertEqual(len(vset), 10)
        vset.clear()
        for i in range(0,10):
            v = self.g.vertex(i)
            vset.add(v)
        self.assertEqual(len(vset), 10)
        self.clear_vertices()

    def add_edges(self):
        # add 10x10 edges to form a complete graph with self-loops
        for i in self.g.vertices():
            for j in self.g.vertices():
                self.g.add_edge(i,j)

    def test_add_edges(self):
        """Adding edges"""
        self.add_vertices()
        self.add_edges()
        N = self.g.num_vertices()
        self.assertEqual(self.g.num_edges(), N**2)
        self.clear_vertices()

    def test_iterate_edges(self):
        """Iterating through edges"""
        self.add_vertices()
        self.add_edges()
        eset = set()
        # iterate through edges and make sure the target and sources make sense
        for e in self.g.edges():
            eset.add(e)
            src = e.source()
            tgt = e.target()
            self.assert_(e in src.out_edges())
            self.assert_(tgt in [w.target() for w in src.out_edges()])
            if self.g.is_directed():
                self.assert_(e in tgt.in_edges())
                self.assert_(src in [w.source() for w in tgt.in_edges()])
            else:
                self.assert_(e in tgt.out_edges())
                self.assert_(src in [w.target() for w in tgt.out_edges()])
        self.assertEqual(len(eset), self.g.num_vertices()**2)
        # iterate through vertices and make sure the out- and in-edges make
        # sense
        for v in self.g.vertices():
            for e in v.out_edges():
                self.assertEqual(e.source(), v)
            for e in v.in_edges():
                self.assertEqual(e.target(), v)
        self.clear_vertices()

    def test_out_of_bounds(self):
        """Testing out of bound conditions"""
        # accessing non existing vertices and edges
        self.assertRaises(StopIteration, lambda: self.g.vertices().next())
        self.assertRaises(StopIteration, lambda: self.g.edges().next())
        self.assertRaises(GraphError, lambda: self.g.vertex(1000))
        self.add_vertices()
        vprop = self.g.add_vertex_property("proptest", "double")
        v = self.g.vertex(5)
        self.g.remove_vertex(v)

        # v is no more...
        self.assert_(v.is_valid() == False)
        self.assertRaises(GraphError, lambda: v.in_degree())
        self.assertRaises(GraphError, lambda: v.out_degree())
        self.assertRaises(GraphError, lambda: v.in_edges())
        self.assertRaises(GraphError, lambda: v.out_edges())
        self.assertRaises(GraphError, lambda: vprop[v])

        self.add_edges()
        eprop = self.g.add_edge_property("proptest", "double")
        e = self.g.edges().next()
        self.g.remove_edge(e)

        # e is no more...
        self.assert_(e.is_valid() == False)
        self.assertRaises(GraphError, lambda: e.target())
        self.assertRaises(GraphError, lambda: e.source())
        self.assertRaises(GraphError, lambda: eprop[e])

    def test_directed(self):
        """Sanity check for directed graph"""

        self.assertEqual(self.g.is_directed(), True)
        self.add_vertices()
        self.add_edges()
        v1 = self.g.vertex(4)
        v2 = self.g.vertex(5)
        e = None
        edel = set()
        for w in self.g.edges():
            if w.source() == v1 and w.target() == v2:
                e = w
            else:
                edel.add(w)
        for w in edel:
            self.g.remove_edge(w, reindex=False)
        self.g.reindex_edges()

        self.assert_(v1 not in [e.target() for e in v2.out_edges()])
        self.assert_(v1 in [e.source() for e in v2.in_edges()])
        self.assert_(v2 in [e.target() for e in v1.out_edges()])
        self.assert_(v2 not in [e.source() for e in v1.in_edges()])
        self.assert_(v1.in_degree() == 0 and v2.out_degree() == 0)
        self.assert_(v1.out_degree() == 1 and v2.in_degree() == 1)
        self.clear_vertices()

    def test_graph_copy(self):
        """Testing graph copying"""
        self.add_vertices()
        self.add_edges()
        orig_vprop = self.g.add_vertex_property("foo", "vector<string>")
        orig_vprop[self.g.vertex(0)] = ["apple", "orange"]
        orig_eprop = self.g.add_edge_property("foo", "bool")
        orig_eprop[self.g.edges().next()] = True
        self.g.add_graph_property("foo", "int32_t", 3)

        g_copy = Graph(self.g)

        vprop = g_copy.vertex_properties["foo"]
        eprop = g_copy.edge_properties["foo"]

        self.assertEqual(vprop.value_type(), orig_vprop.value_type())
        self.assertEqual(list(vprop[g_copy.vertex(0)]), list(orig_vprop[self.g.vertex(0)]))
        vprop[g_copy.vertex(0)] = []
        self.assertNotEqual(list(vprop[g_copy.vertex(0)]),
                            list(orig_vprop[self.g.vertex(0)]))

        self.assertEqual(eprop.value_type(), orig_eprop.value_type())
        self.assertEqual(eprop[g_copy.edges().next()],
                         orig_eprop[self.g.edges().next()])
        eprop[g_copy.edges().next()] = False
        self.assertNotEqual(eprop[g_copy.edges().next()],
                            orig_eprop[g_copy.edges().next()])

        self.assertEqual(self.g.get_graph_property("foo"),
                         g_copy.get_graph_property("foo"))
        g_copy.set_graph_property("foo", 42)
        self.assertNotEqual(self.g.get_graph_property("foo"),
                            g_copy.get_graph_property("foo"))

        self.assertEqual(self.g.is_directed(), g_copy.is_directed())
        self.assertEqual(g_copy.get_vertex_filter(), None)
        self.assertEqual(g_copy.is_reversed(), False)
        self.assertEqual(self.g.num_vertices(), g_copy.num_vertices())
        self.assertEqual(self.g.num_edges(), g_copy.num_edges())

        g_copy.remove_vertex(g_copy.vertex(0))
        self.assertEqual(self.g.num_vertices() - 1, g_copy.num_vertices())
        self.assertEqual(self.g.num_edges() - 19, g_copy.num_edges())

        self.clear_vertices()

def get_undirected_test(TestCase):
    class TestUndirectedManipulation(TestCase):
        """Test basic undirected graph manipulation"""
        def id(self):
            """Properly name our tests"""
            name = \
                 eval(TestCase.id(self).\
                      replace("graph_tool.test.basic.TestUndirectedManipulation",
                              "self").\
                      replace("graph_tool.test.basic.FilteredTestCase",
                              "self") + ".__doc__")
            return name + " (undirected)"
        def shortDescription(self):
            return self.id()
        def setUp(self):
            self.g = Graph()
            self.g.set_directed(False)

        def test_directed(self): # tests if it is _undirected_
            """Sanity check for undirected graph"""
            self.assertEqual(self.g.is_directed(), False)
            self.add_vertices()
            self.add_edges()
            v1 = self.g.vertex(4)
            v2 = self.g.vertex(5)
            e = None
            edel = set()
            for w in self.g.edges():
                if w.source() == v1 and w.target() == v2:
                    e = w
                else:
                    edel.add(w)
            for w in edel:
                self.g.remove_edge(w, reindex=False)
            self.g.reindex_edges()
            self.assert_(v1 in [e.target() for e in v2.out_edges()] and
                         v2 in [e.target() for e in v1.out_edges()])
            self.assert_(v1.in_degree() == 0 and v1.in_degree() == 0)
            self.clear_vertices()
    return TestUndirectedManipulation

def get_reversed_test(TestCase):
    class TestReversedManipulation(TestCase):
        """Test basic reversed graph manipulation"""
        def id(self):
            """Properly name our tests"""
            name = eval(TestCase.id(self).\
                        replace("graph_tool.test.basic.TestReversedManipulation",
                                "self").\
                      replace("graph_tool.test.basic.FilteredTestCase",
                              "self") + ".__doc__")
            return name + " (reversed)"

        def shortDescription(self):
            return self.id()

        def setUp(self):
            self.g = Graph()
            self.g.set_reversed(True)
            self.assertEqual(self.g.is_reversed(), True)

        def testReverse(self):
            """Sanity check for reversed graphs"""

            self.add_vertices()
            self.add_edges()
            v1 = self.g.vertex(4)
            v2 = self.g.vertex(5)
            e = None
            edel = set()
            for w in self.g.edges():
                if w.source() == v1 and w.target() == v2:
                    e = w
                else:
                    edel.add(w)
            for w in edel:
                self.g.remove_edge(w, reindex=False)
            self.g.reindex_edges()
            self.assert_(e.target() == v2 and e.source() == v1)
            self.g.set_reversed(not self.g.is_reversed())
            self.assert_(e.target() == v1 and e.source() == v2)
    return TestReversedManipulation

def get_filtered_test(BaseTestCase, edge_filter=True, vertex_filter=True,
                      invert_edge=False, invert_vertex=False):
    class FilteredTestCase(BaseTestCase):
        """Testing filtered graphs"""
        def id(self):
            """Properly name our tests"""
            if "FilteredTestCase" in BaseTestCase.id(self):
                name = eval(BaseTestCase.id(self).\
                            replace("graph_tool.test.basic.FilteredTestCase",
                                    "self") +\
                            ".__doc__")
            else:
                name = BaseTestCase.id(self)

            filter_name = "[ "
            filter_name += "edge_filter " if edge_filter else ""
            filter_name += "(inverted) " if invert_edge else ""
            filter_name += "vertex_filter " if vertex_filter else ""
            filter_name += "(inverted) " if invert_vertex else ""
            filter_name += "]"
            return name + " " + (60-len(name))*" " + filter_name
        def shortDescription(self):
             return self.id()
        def setUp(self):
            #print self.id()
            self.g = Graph()
            BaseTestCase.setUp(self)
            if vertex_filter:
                self.vprop = self.g.add_vertex_property("filter_prop","bool")
                self.g.set_vertex_filter("filter_prop", invert_vertex)
            if edge_filter:
                self.eprop = self.g.add_edge_property("filter_prop","bool")
                self.g.set_edge_filter("filter_prop", invert_edge)

        def add_vertices(self):
            BaseTestCase.add_vertices(self)
            if vertex_filter:
                self.g.reset_vertex_filter()
                for v in self.g.vertices():
                    self.vprop[v] = True if not invert_vertex else False
                self.g.set_vertex_filter("filter_prop", invert_vertex)

        def add_edges(self):
            BaseTestCase.add_edges(self)
            if edge_filter:
                self.g.reset_edge_filter()
                for e in self.g.edges():
                    self.eprop[e] = True if not invert_edge else False
                self.g.set_edge_filter("filter_prop", invert_edge)

        def test_filtering(self):
            "Test explicit vertex and/or edge filtering"
            self.add_vertices()
            self.add_edges()

            if edge_filter:
                edges = [e for e in self.g.edges()]
                self.assert_(len(edges) > 0)
                e = edges[5]
                v = e.source()
                ki = v.out_degree()
                EI = self.g.num_edges()
                self.eprop[e] = not self.eprop[e] # filter out e
                EF = self.g.num_edges()
                kf = v.out_degree()
                self.assertEqual(EF, EI - 1)
                self.assertEqual(kf, ki - 1)
                self.assertEqual(e.source(), v) # filtered edge is still
                                                # functional
            if vertex_filter:
                EI = self.g.num_edges()
                VI = self.g.num_vertices()
                v = self.g.vertex(5)
                ki = v.out_degree() + v.in_degree()
                self.vprop[v] = not self.vprop[v] # filter out v
                EF = self.g.num_edges()
                VF = self.g.num_vertices()
                kf = v.out_degree() + v.in_degree()
                self.assertEqual(EF, EI - ki + 1)
                self.assertEqual(VF, VI - 1)
                self.assertEqual(kf, ki - 2)
            self.clear_vertices()
    return FilteredTestCase

def gen_suite(TestCase):
    s = unittest.TestLoader().loadTestsFromTestCase(TestCase)
    s.addTest(unittest.TestLoader().\
              loadTestsFromTestCase(get_undirected_test(TestCase)))
    s.addTest(unittest.TestLoader().\
              loadTestsFromTestCase(get_reversed_test(TestCase)))
    for vertex_filtered in [True,False]:
        for edge_filtered in [True,False]:
            if vertex_filtered == edge_filtered == False:
                continue
            for invert_vertex in [True, False]:
                if not vertex_filtered and invert_vertex:
                    continue
                for invert_edge in [True, False]:
                    if not edge_filtered and invert_edge:
                        continue
                    testcase = get_filtered_test(TestCase,
                                                 edge_filtered, vertex_filtered,
                                                 invert_edge, invert_vertex)
                    s.addTest(unittest.TestLoader().\
                              loadTestsFromTestCase(testcase))
                    testcase = get_filtered_test\
                               (get_undirected_test(TestCase),
                                edge_filtered, vertex_filtered,
                                invert_edge, invert_vertex)
                    s.addTest(unittest.TestLoader().\
                              loadTestsFromTestCase(testcase))
                    testcase = get_filtered_test\
                               (get_reversed_test(TestCase),
                                edge_filtered, vertex_filtered,
                                invert_edge, invert_vertex)
                    s.addTest(unittest.TestLoader().\
                              loadTestsFromTestCase(testcase))
    return s

def suite():
    return gen_suite(TestBasicManipulation)
