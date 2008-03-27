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
Basic property map manipulation tests
"""

from graph_tool import *
import unittest
from numpy import *

from basic import gen_suite

value_types = ["bool", "int32_t", "int64_t", "double", "long double", "string",
               "vector<bool>", "vector<int32_t>", "vector<int64_t>",
               "vector<double>", "vector<long double>", "vector<string>",
               "python::object"]

type_equivalence = { "bool": int, # bool,
                     "int32_t": int,
                     "int64_t": long,
                     "double": float,
                     "long double": float,
                     "string": str,
                     "vector<bool>": Vector_bool,
                     "vector<int32_t>": Vector_int32_t,
                     "vector<int64_t>": Vector_int64_t,
                     "vector<double>": Vector_double,
                     "vector<long double>": Vector_long_double,
                     "vector<string>": Vector_string,
                     "python::object": dict #any object really
                     }

type_init = { "bool": True,
              "int32_t": int(42),
              "int64_t": long(42),
              "double": float(3.1416),
              "long double": float(3.1416),
              "string": "foo",
              "vector<bool>": [True, False, True],
              "vector<int32_t>": [1, 2, 3, 42],
              "vector<int64_t>": [long(1), long(2), long(99999999999999)],
              "vector<double>": [1.1, 1.2, 1.3, 1.4] ,
              "vector<long double>": [1.1, 1.2, 1.3, 1.4],
              "vector<string>": [ "lambs", "stoats", "orangutans",
                                  "breakfast cereals", "lima beans" ],
              "python::object": {3:[1,2,3], "foo":(3,2,42) }}

class TestPropertyManipulation(unittest.TestCase):
    """Test property map manipulation"""

    def setUp(self):
        self.g = Graph()

    def add_vertices(self):
        # add 10 vertices
        for i in range(0,10):
            v = self.g.add_vertex()
            self.assert_(str(v) == str(i))

    def clear_vertices(self):
        self.g.clear()

    def add_edges(self):
        # add 10x10 edges to form a complete graph with self-loops
        for i in self.g.vertices():
            for j in self.g.vertices():
                self.g.add_edge(i,j)

    def test_property_manipulation(self):
        """Testing property map creation and manipulation"""
        self.add_vertices()
        self.add_edges()

        # test properties of all types
        for t in value_types:
            # vertex properties
            self.g.add_vertex_property("prop_" + t, t)
            self.assertEqual(self.g.vertex_properties["prop_" + t].value_type(),
                             t)
            prop = self.g.vertex_properties["prop_" + t]
            for v in self.g.vertices():
                prop[v] = type_init[t]
                val = prop[v]
                self.assertEqual(type(val), type_equivalence[t])
                if t.startswith("vector"):
                    self.assertEqual(list(prop[v]), type_init[t])
                else:
                    self.assertEqual(prop[v], type_init[t])

            # edge properties
            prop = self.g.add_edge_property("prop_" + t, t)
            self.assertEqual(self.g.edge_properties["prop_" + t].value_type(),t)
            for e in self.g.edges():
                prop[e] = type_init[t]
                val = prop[e]
                self.assertEqual(type(val), type_equivalence[t])
                if t.startswith("vector"):
                    self.assertEqual(list(prop[e]), type_init[t])
                else:
                    self.assertEqual(prop[e], type_init[t])

            # graph properties
            self.g.add_graph_property("prop_" + t, t)
            self.g.graph_properties["prop_" + t] = type_init[t]
            val = self.g.graph_properties["prop_" + t]
            self.assertEqual(type(val), type_equivalence[t])
            val = self.g.graph_properties["prop_" + t]
            if t.startswith("vector"):
                self.assertEqual(list(val), type_init[t])
            else:
                self.assertEqual(val, type_init[t])

            # test property deletion
            del self.g.vertex_properties["prop_" + t]
            self.assertRaises(KeyError,
                              lambda: self.g.vertex_properties["prop_" + t])

            del self.g.edge_properties["prop_" + t]
            self.assertRaises(KeyError,
                              lambda: self.g.edge_properties["prop_" + t])

            del self.g.graph_properties["prop_" + t]
            self.assertRaises(KeyError,
                              lambda: self.g.graph_properties["prop_" + t])
        self.clear_vertices()

    def test_vector_property_manipulation(self):
        """Testing vector property manipulation"""
        self.add_vertices()
        self.add_edges()

        def test_prop(add_func, get_func, set_func, get_desc,
                      get_prop = lambda p,v: p[v]):
            add_func(t)
            prop = get_func(t)
            try:
                prop.__getattribute__("value_type")
                self.assertEqual(get_func(t).value_type(), t)
            except AttributeError:
                pass
            v = get_desc()
            vec = get_prop(prop,v)
            set_func(t,v,type_init[t])
            self.assertEqual(get_prop(get_func(t),v), vec)
            vec[0] = 42
            self.assertEqual(vec[0], 42)
            self.assertEqual(get_prop(get_func(t),v)[0], vec[0])

            # test numpy integration
            try:
                prop.__getattribute__("value_type")
                if "string" not in str(get_func(t).value_type()):
                    a = vec.ndarray()
                    self.assertEqual(len(vec), a.size)
                    self.assertEqual(vec[0], a[0])
                    self.assertEqual(list(vec), list(a))
                    a[0] = 21
                    self.assertEqual(vec[0], a[0])
                    self.assertEqual(list(vec), list(a))
            except AttributeError:
                pass

        def set_vertex_prop(t,v,val):
            self.g.vertex_properties["prop_" + t][v] = val
        def set_edge_prop(t,e,val):
            self.g.edge_properties["prop_" + t][e] = val
        def set_graph_prop(t,v,val):
            self.g.graph_properties["prop_" + t] = val

        # test properties of all types
        for t in [val for val in value_types if val.startswith("vector")\
                  and ("string" not in val)]:
            # vertex properties
            test_prop(lambda t: self.g.add_vertex_property("prop_" + t, t),
                      lambda t: self.g.vertex_properties["prop_" + t],
                      set_vertex_prop,
                      lambda: self.g.vertex(0))

            # edge properties
            test_prop(lambda t: self.g.add_edge_property("prop_" + t, t),
                      lambda t: self.g.edge_properties["prop_" + t],
                      set_edge_prop,
                      lambda: self.g.edges().next())

            # graph properties
            test_prop(lambda t: self.g.add_graph_property("prop_" + t, t),
                      lambda t: self.g.graph_properties["prop_" + t],
                      set_graph_prop,
                      lambda: self.g,
                      lambda p,v: p)
        self.clear_vertices()


def suite():
    return gen_suite(TestPropertyManipulation)
