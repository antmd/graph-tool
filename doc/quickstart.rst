Quick start using `graph-tool`
==============================

The :mod:`graph_tool` module provides a :class:`~graph_tool.Graph` class and
several algorithms that operate on it. The internals of this class, and of most
algorithms, are written in C++ for performance.

The module must be of course imported before it can be used. The package is
subdivided into several sub-modules. To import everything from all of them, one
can do:

.. doctest::

   >>> from graph_tool.all import *

In the following, it will always be assumed that the previous line was run.

Creating and manipulating graphs
--------------------------------

An empty graph can be created by instantiating a :class:`~graph_tool.Graph`
class:

.. doctest::

   >>> g = Graph()

By default, newly created graphs are always directed. To construct undirected
graphs, one must pass the ``directed`` parameter:

.. doctest::

   >>> ug = Graph(directed=False)

A graph can always be switched on-the-fly from directed to undirected (and vice
versa), with the :meth:`~graph_tool.Graph.set_directed` method. The "directedness" of the
graph can be queried with the :meth:`~graph_tool.Graph.is_directed` method,

.. doctest::

   >>> ug = Graph()
   >>> ug.set_directed(False)
   >>> assert(ug.is_directed() == False)

A graph can also be created from another graph, in which case the entire graph
(and its internal property maps, see :ref:`sec_property_maps`) is copied:

.. doctest::

   >>> g1 = Graph()
   >>> # ... populate g1 ...
   >>> g2 = Graph(g1)                 # g1 and g2 are copies

Once a graph is created, it can be populated with vertices and edges. A vertex
can be added with the :meth:`~graph_tool.Graph.add_vertex` method,

.. doctest::

   >>> v = g.add_vertex()

which returns an instance of a :class:`~graph_tool.Vertex` class, also called a
*vertex descriptor*. The :meth:`~graph_tool.Graph.add_vertex` method also
accepts an optional parameter which specifies the number of vertices to
create. If this value is greater than 1, it returns a list of vertices:

.. doctest::

   >>> vlist = g.add_vertex(10)
   >>> print len(vlist)
   10

Each vertex has an unique index, which is numbered from 0 to N-1, where N is the
number of vertices. This index can be obtained by using the
:attr:`~graph_tool.Graph.vertex_index` attribute of the graph (which is a
*property map*, see :ref:`sec_property_maps`), or by converting the vertex
descriptor to an int.

.. doctest::

   >>> v = g.add_vertex()
   >>> print g.vertex_index[v], int(v)
   11 11

There is no need to keep the vertex descriptor lying around to access them at a
later point: One can obtain the descriptor of a vertex with a given index using
the :meth:`~graph_tool.Graph.vertex` method,

.. doctest::

   >>> print g.vertex(8)
   8

Another option is to iterate through the vertices, as described in section
:ref:`sec_iteration`.

Once we have some vertices in the graph, we can create some edges between them
with the :meth:`~graph_tool.Graph.add_edge` method, which returns an edge
descriptor (an instance of the :class:`~graph_tool.Edge` class).

.. doctest::

   >>> v1 = g.add_vertex()
   >>> v2 = g.add_vertex()
   >>> e = g.add_edge(v1, v2)

Edges also have an unique index, which is given by the :attr:`~graph_tool.Graph.edge_index`
property:

.. doctest::

   >>> print g.edge_index[e]
   0

Unlike the vertices, edge indexes are not guaranteed to be continuous in any
range, but they are always unique.

Both vertex and edge descriptors have methods which query associate information,
such as :meth:`~graph_tool.Vertex.in_degree`,
:meth:`~graph_tool.Vertex.out_degree`, :meth:`~graph_tool.Edge.source` and
:meth:`~graph_tool.Edge.target`:

.. doctest::

   >>> v1 = g.add_vertex()
   >>> v2 = g.add_vertex()
   >>> e = g.add_edge(v1, v2)
   >>> print v1.out_degree(), v2.in_degree()
   1 1
   >>> assert(e.source() == v1 and e.target() == v2)

Edges and vertices can also be removed at any time with the
:meth:`~graph_tool.Graph.remove_vertex` and :meth:`~graph_tool.Graph.remove_edge` methods,

.. doctest::

   >>> e = g.add_edge(g.vertex(0), g.vertex(1))
   >>> g.remove_edge(e)                                      # e no longer exists
   >>> g.remove_vertex(g.vertex(1))              # the second vertex is also gone

.. _sec_iteration:

Iterating over vertices and edges
+++++++++++++++++++++++++++++++++

Algorithms must often iterate through the vertices, edges, out edge, etc. of the
graph. The :class:`~graph_tool.Graph` and :class:`~graph_tool.Edge` classes
provide the necessary iterators for doing so. The iterators always give back
edge or vertex descriptors.

In order to iterate through the vertices or edges of the graph, the
:meth:`~graph_tool.Graph.vertices` and :meth:`~graph_tool.Graph.edges` methods should be used, as such:

.. doctest::

   for v in g.vertices():
       print v
   for e in e.vertices():
       print e

The code above will print the vertices and edges of the graph in the order they
are found.

The out- and in-edges of a vertex, as well as the out- and in-neighbours can be
iterated through with the :meth:`~graph_tool.Vertex.out_edges`,
:meth:`~graph_tool.Vertex.in_edges`, :meth:`~graph_tool.Vertex.out_neighbours`
and :meth:`~graph_tool.Vertex.in_neighbours` respectively.

.. doctest::

   from itertools import izip
   for v in g.vertices():
      for e in v.out_edges():
          print e
      for e in v.out_neighbours():
          print e

      # the edge and neighbours order always match
      for e,w in izip(v.out_edges(), v.out_neighbours()):
          assert(e.target() == w)

.. warning:

   You should never remove vertex or edge descriptors when iterating over them,
   since this invalidates the iterators. If you plan to remove vertices or edges
   during iteration, you must first store them somewhere (such as in a list) and
   remove them only later. Removal during iteration will cause bad things to
   happen.

.. _sec_property_maps:

Property maps
-------------

Property maps are a way of associating additional to the vertices, edges or to
the graph itself. There are thus three types of property maps: vertex, edge and
graph. All of them are instances of the same class,
:class:`~graph_tool.PropertyMap`. Each property map has an associated *value
type*, which must be chosen from the predefined set:

.. tabularcolumns:: |l|l|

.. table::

    =======================     ================
     Type name                  Aliases
    =======================     ================
    ``bool``                    ``uint8_t``
    ``int32_t``                 ``int``   
    ``int64_t``                 ``long``  
    ``double``                  ``float``
    ``long double``             
    ``string``                  
    ``vector<bool>``            ``vector<uint8_t>``
    ``vector<int32_t>``         ``vector<int>``
    ``vector<int64_t>``         ``vector<long>``
    ``vector<double>``          ``vector<float>``
    ``vector<long double>``
    ``vector<string>``
    ``python::object``          ``object`` 
    =======================     ================

New property maps can be created for a given graph by calling the
:meth:`~graph_tool.Graph.new_vertex_property`, :meth:`~graph_tool.Graph.new_edge_property`, or
:meth:`~graph_tool.Graph.new_graph_property`, for each map type. The values are then
accessed by vertex or edge descriptors, or the graph itself, as such:

.. doctest::

    from itertools import izip
    from numpy.random import randint

    g = Graph()
    g.add_vertex(100)
    # insert some random links
    for s,t in izip(randint(0, 100, 100), randint(0, 100, 100)):
        g.add_edge(g.vertex(s), g.vertex(t))
    
    vprop_double = g.new_vertex_property("double")
    vprop_vint = g.new_vertex_property("vector<int>")

    eprop_dict = g.new_edge_property("object")

    gprop_bool = g.new_edge_property("bool")

    vprop_double[g.vertex(10)] = 3.1416

    vprop_vint[g.vertex(40)] = [1, 3, 42, 54]
    
    eprop_dict[g.edges().next()] = {"foo":"bar", "gnu":42}

    gprop_bool[g] = True

Property maps with scalar value types can also be accessed as a numpy
:class:`~numpy.ndarray`, with the :meth:`~graph_tool.PropertyMap.get_array`
method, i.e.,

.. doctest::

    from numpy.random import random

    # this assigns random values to the properties
    vprop_double.get_array()[:] = random(g.num_vertices())

Internal property maps
++++++++++++++++++++++

Any created property map can be made "internal" to the respective graph. This
means that it will be copied and saved to a file together with the
graph. Properties are internalized by including them in the graph's
dictionary-like attributes :attr:`~graph_tool.Graph.vertex_properties`,
:attr:`~graph_tool.Graph.edge_properties` or
:attr:`~graph_tool.Graph.graph_properties`. When inserted in the graph, the
property maps must have an unique name (between those of the same type):

.. doctest::

    >>> eprop = g.new_edge_property("string")
    >>> g.edge_properties["some name"] = eprop
    >>> g.list_properties()
    some name      (edge)    (type: string)


Graph I/O
---------

Graphs can be saved and loaded in two formats: `graphml
<http://graphml.graphdrawing.org/>`_ and `dot
<http://www.graphviz.org/doc/info/lang.html>`_. Graphml is the default and
preferred format. The dot format is also supported, but since it contains no
type information, all properties are read later as strings, and must be
converted per hand. Therefore you should always use graphml, except when
interfacing with another software which expects dot format.

A graph can be saved or loaded to a file with the :attr:`~graph_tool.Graph.save`
and :attr:`~graph_tool.Graph.load` methods, which take either a file name or a
file-like object. A graph can also be loaded from disk with the
:func:`~graph_tool.load_graph` function, as such:

.. doctest::

    g = Graph()
    #  ... fill the graph ...
    g.save("my_graph.xml.gz")    
    g2 = load_graph("my_graph.xml.gz")
    # g and g2 should be a copy of each other

Graph classes can also be pickled with the :mod:`pickle` module.


An Example: Building a Price Network
------------------------------------

.. literalinclude:: price.py
   :linenos:

.. testcode::
   :hide:

   from price import *
   clf()

.. testoutput::

    vertex: 36063 in-degree: 0 out-degree: 1 age: 36063
    vertex: 9075 in-degree: 4 out-degree: 1 age: 9075
    vertex: 5967 in-degree: 3 out-degree: 1 age: 5967
    vertex: 1113 in-degree: 7 out-degree: 1 age: 1113
    vertex: 25 in-degree: 84 out-degree: 1 age: 25
    vertex: 10 in-degree: 541 out-degree: 1 age: 10
    vertex: 5 in-degree: 140 out-degree: 1 age: 5
    vertex: 2 in-degree: 459 out-degree: 1 age: 2
    vertex: 1 in-degree: 520 out-degree: 1 age: 1
    vertex: 0 in-degree: 210 out-degree: 0 age: 0
    Nowhere else to go... We found the main hub!


.. image:: deg-hist.png

