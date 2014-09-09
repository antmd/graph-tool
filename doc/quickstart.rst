Quick start using `graph-tool`
==============================

The :mod:`graph_tool` module provides a :class:`~graph_tool.Graph` class
and several algorithms that operate on it. The internals of this class,
and of most algorithms, are written in C++ for performance, using the
`Boost Graph Library <http://www.boost.org>`_.

The module must be of course imported before it can be used. The package is
subdivided into several sub-modules. To import everything from all of them, one
can do:

.. testsetup::

   np.random.seed(42)
   gt.seed_rng(42)

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
graphs, one must pass a value to the ``directed`` parameter:

.. doctest::

   >>> ug = Graph(directed=False)

A graph can always be switched *on-the-fly* from directed to undirected
(and vice versa), with the :meth:`~graph_tool.Graph.set_directed`
method. The "directedness" of the graph can be queried with the
:meth:`~graph_tool.Graph.is_directed` method,

.. doctest::

   >>> ug = Graph()
   >>> ug.set_directed(False)
   >>> assert(ug.is_directed() == False)

A graph can also be created by providing another graph, in which case
the entire graph (and its internal property maps, see
:ref:`sec_property_maps`) is copied:

.. doctest::

   >>> g1 = Graph()
   >>> # ... construct g1 ...
   >>> g2 = Graph(g1)                 # g1 and g2 are copies

Above, ``g2`` is a "deep" copy of ``g1``, i.e. any modification of
``g2`` will not affect ``g1``.

Once a graph is created, it can be populated with vertices and edges. A
vertex can be added with the :meth:`~graph_tool.Graph.add_vertex`
method, which returns an instance of a :class:`~graph_tool.Vertex`
class, also called a *vertex descriptor*. For instance, the following
code creates two vertices, and returns vertex descriptors stored in the
variables ``v1`` and ``v2``.

.. doctest::

   >>> v1 = g.add_vertex()
   >>> v2 = g.add_vertex()

Edges can be added in an analogous manner, by calling the
:meth:`~graph_tool.Graph.add_edge` method, which returns an edge
descriptor (an instance of the :class:`~graph_tool.Edge` class):

.. doctest::

   >>> e = g.add_edge(v1, v2)

The above code creates a directed edge from ``v1`` to ``v2``. We can
visualize the graph we created so far with the
:func:`~graph_tool.draw.graph_draw` function.

.. doctest::

   >>> graph_draw(g, vertex_text=g.vertex_index, vertex_font_size=18,
   ...            output_size=(200, 200), output="two-nodes.png")
   <...>

.. doctest::
   :hide:

   graph_draw(g, vertex_text=g.vertex_index, vertex_font_size=18,
              output_size=(200, 200), output="two-nodes.pdf")


.. figure:: two-nodes.*
   :align: center

   A simple directed graph with two vertices and one edge, created by
   the commands above.

With vertex and edge descriptors, one can examine and manipulate the
graph in an arbitrary manner. For instance, in order to obtain the
out-degree of a vertex, we can simply call the
:meth:`~graph_tool.Vertex.out_degree` method:

.. doctest::

   >>> print(v1.out_degree())
   1

Analogously, we could have used the :meth:`~graph_tool.Vertex.in_degree`
method to query the in-degree.

.. note::

   For undirected graphs, the "out-degree" is synonym for degree, and
   in this case the in-degree of a vertex is always zero.

Edge descriptors have two useful methods, :meth:`~graph_tool.Edge.source`
and :meth:`~graph_tool.Edge.target`, which return the source and target
vertex of an edge, respectively.

.. doctest::

   >>> print(e.source(), e.target())
   0 1

The :meth:`~graph_tool.Graph.add_vertex` method also accepts an optional
parameter which specifies the number of vertices to create. If this
value is greater than 1, it returns an iterator on the added vertex
descriptors:

.. doctest::

   >>> vlist = g.add_vertex(10)
   >>> print(len(list(vlist)))
   10

Edges and vertices can also be removed at any time with the
:meth:`~graph_tool.Graph.remove_vertex` and :meth:`~graph_tool.Graph.remove_edge` methods,

.. doctest::

   >>> g.remove_edge(e)                               # e no longer exists
   >>> g.remove_vertex(v2)                # the second vertex is also gone

.. note::

   Removing a vertex is an :math:`O(N)` operation. The vertices are
   internally stored in a `STL vector <http://en.wikipedia.org/wiki/Vector_%28STL%29>`_,
   so removing an element somewhere in the middle of the list requires
   the shifting of the rest of the list. Thus, fast :math:`O(1)`
   removals are only possible either if one can guarantee that only
   vertices in the end of the list are removed (the ones last added to
   the graph), or if the relative vertex ordering is invalidated. This
   last behavior can be achieved by passing the option ``fast == True``,
   to :meth:`~graph_tool.Graph.remove_vertex`, which causes vertex
   being deleted to be 'swapped' with the last vertex (i.e. with the
   largest index), which will in turn inherit the index of the vertex
   being deleted.


   Removing an edge is an :math:`O(k_{s} + k_{t})` operation, where
   :math:`k_{s}` is the out-degree of the source vertex, and
   :math:`k_{t}` is the in-degree of the target vertex. This can be made
   faster by setting :meth:`~graph_tool.Graph.set_fast_edge_removal` to
   `True`, in which case it becomes :math:`O(1)`, at the expense of
   additional data of size :math:`O(E)`.

Each vertex in a graph has an unique index, which is numbered from 0 to
N-1, where N is the number of vertices. This index can be obtained by
using the :attr:`~graph_tool.Graph.vertex_index` attribute of the graph
(which is a *property map*, see :ref:`sec_property_maps`), or by
converting the vertex descriptor to an ``int``.

.. doctest::

   >>> v = g.add_vertex()
   >>> print(g.vertex_index[v])
   11
   >>> print(int(v))
   11

.. note::

   Removing a vertex will cause the index of any vertex with a larger
   index to be changed, so that the indexes *always* conform to the
   :math:`[0, N-1]` range. However, property map values (see
   :ref:`sec_property_maps`) are unaffected.

Since vertices are uniquely identifiable by their indexes, there is no
need to keep the vertex descriptor lying around to access them at a
later point. If we know its index, one can obtain the descriptor of a
vertex with a given index using the :meth:`~graph_tool.Graph.vertex`
method,

.. doctest::

   >>> v = g.vertex(8)

which takes an index, and returns a vertex descriptor. Edges cannot be
directly obtained by its index, but if the source and target vertices of
a given edge is known, it can be obtained with the
:meth:`~graph_tool.Graph.edge` method

.. doctest::

   >>> g.add_edge(g.vertex(2), g.vertex(3))
   <...>
   >>> e = g.edge(2, 3)


Another way to obtain edge or vertex descriptors is to *iterate* through
them, as described in section :ref:`sec_iteration`. This is in fact the
most useful way of obtaining vertex and edge descriptors.

Like vertices, edges also have unique indexes, which are given by the
:attr:`~graph_tool.Graph.edge_index` property:

.. doctest::

   >>> e = g.add_edge(g.vertex(0), g.vertex(1))
   >>> print(g.edge_index[e])
   1

Differently from vertices, edge indexes do not necessarily conform to
any specific range. If no edges are ever removed, the indexes will be in
the range :math:`[0, E-1]`, where :math:`E` is the number of edges, and
edges added earlier have lower indexes. However if an edge is removed,
its index will be "vacant", and the remaining indexes will be left
unmodified, and thus will not lie in the range :math:`[0, E-1]`.  If a
new edge is added, it will reuse old indexes, in an increasing order.

.. _sec_iteration:

Iterating over vertices and edges
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Algorithms must often iterate through vertices, edges, out-edges of a
vertex, etc. The :class:`~graph_tool.Graph` and
:class:`~graph_tool.Vertex` classes provide different types of iterators
for doing so. The iterators always point to edge or vertex descriptors.

Iterating over all vertices or edges
""""""""""""""""""""""""""""""""""""

In order to iterate through all the vertices or edges of a graph, the
:meth:`~graph_tool.Graph.vertices` and :meth:`~graph_tool.Graph.edges`
methods should be used:

.. doctest::

   for v in g.vertices():
       print(v)
   for e in g.edges():
       print(e)

The code above will print the vertices and edges of the graph in the order they
are found.

Iterating over the neighbourhood of a vertex
""""""""""""""""""""""""""""""""""""""""""""

The out- and in-edges of a vertex, as well as the out- and in-neighbours can be
iterated through with the :meth:`~graph_tool.Vertex.out_edges`,
:meth:`~graph_tool.Vertex.in_edges`, :meth:`~graph_tool.Vertex.out_neighbours`
and :meth:`~graph_tool.Vertex.in_neighbours` methods, respectively.

.. doctest::

   from itertools import izip
   for v in g.vertices():
      for e in v.out_edges():
          print(e)
      for w in v.out_neighbours():
          print(w)

      # the edge and neighbours order always match
      for e,w in izip(v.out_edges(), v.out_neighbours()):
          assert(e.target() == w)

The code above will print the out-edges and out-neighbours of all
vertices in the graph.

.. note::

   The ordering of the vertices and edges visited by the iterators is
   always the same as the order in which they were added to the graph
   (with the exception of the iterator returned by
   :meth:`~graph_tool.Graph.edges`). Usually, algorithms do not care
   about this order, but if it is ever needed, this inherent ordering
   can be relied upon.

.. warning::

   You should never remove vertex or edge descriptors when iterating
   over them, since this invalidates the iterators. If you plan to
   remove vertices or edges during iteration, you must first store them
   somewhere (such as in a list) and remove them only after no iterator
   is being used. Removal during iteration will cause bad things to
   happen.


.. _sec_property_maps:

Property maps
-------------

Property maps are a way of associating additional information to the
vertices, edges or to the graph itself. There are thus three types of
property maps: vertex, edge and graph. All of them are handled by the
same class, :class:`~graph_tool.PropertyMap`. Each created property map
has an associated *value type*, which must be chosen from the predefined
set:

.. tabularcolumns:: |l|l|

.. table::

    ========================     ======================
     Type name                   Alias
    ========================     ======================
    ``bool``                     ``uint8_t``
    ``int16_t``                  ``short``
    ``int32_t``                  ``int``
    ``int64_t``                  ``long``, ``long long``
    ``double``                   ``float``
    ``long double``
    ``string``
    ``vector<bool>``             ``vector<uint8_t>``
    ``vector<int16_t>``          ``vector<short>``
    ``vector<int32_t>``          ``vector<int>``
    ``vector<int64_t>``          ``vector<long>``, ``vector<long long>``
    ``vector<double>``           ``vector<float>``
    ``vector<long double>``
    ``vector<string>``
    ``python::object``           ``object``
    ========================     ======================

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

    vprop_double = g.new_vertex_property("double")            # Double-precision floating point
    vprop_double[g.vertex(10)] = 3.1416

    vprop_vint = g.new_vertex_property("vector<int>")         # Vector of ints
    vprop_vint[g.vertex(40)] = [1, 3, 42, 54]
    
    eprop_dict = g.new_edge_property("object")                # Arbitrary python object.
    eprop_dict[g.edges().next()] = {"foo": "bar", "gnu": 42}  # In this case, a dict.

    gprop_bool = g.new_edge_property("bool")                  # Boolean
    gprop_bool[g] = True

Property maps with scalar value types can also be accessed as a
:class:`numpy.ndarray`, with the
:meth:`~graph_tool.PropertyMap.get_array` method, or the
:attr:`~graph_tool.PropertyMap.a` attribute, i.e.,

.. doctest::

    from numpy.random import random

    # this assigns random values to the vertex properties
    vprop_double.get_array()[:] = random(g.num_vertices())

    # or more conveniently (this is equivalent to the above)
    vprop_double.a = random(g.num_vertices())

.. _sec_internal_props:

Internal property maps
^^^^^^^^^^^^^^^^^^^^^^

Any created property map can be made "internal" to the corresponding
graph. This means that it will be copied and saved to a file together
with the graph. Properties are internalized by including them in the
graph's dictionary-like attributes
:attr:`~graph_tool.Graph.vertex_properties`,
:attr:`~graph_tool.Graph.edge_properties` or
:attr:`~graph_tool.Graph.graph_properties` (or their aliases,
:attr:`~graph_tool.Graph.vp`, :attr:`~graph_tool.Graph.ep` or
:attr:`~graph_tool.Graph.gp`, respectively). When inserted in the graph,
the property maps must have an unique name (between those of the same
type):

.. doctest::

    >>> eprop = g.new_edge_property("string")
    >>> g.edge_properties["some name"] = eprop
    >>> g.list_properties()
    some name      (edge)    (type: string)

Internal graph property maps behave slightly differently. Instead of
returning the property map object, the value itself is returned from the
dictionaries:

.. doctest::

    >>> gprop = g.new_graph_property("int")
    >>> g.graph_properties["foo"] = gprop   # this sets the actual property map
    >>> g.graph_properties["foo"] = 42      # this sets its value
    >>> print(g.graph_properties["foo"])
    42
    >>> del g.graph_properties["foo"]       # the property map entry is deleted from the dictionary


.. _sec_graph_io:

Graph I/O
---------

Graphs can be saved and loaded in four formats: `graphml
<http://graphml.graphdrawing.org/>`_, `dot
<http://www.graphviz.org/doc/info/lang.html>`_, `gml
<http://www.fim.uni-passau.de/en/fim/faculty/chairs/theoretische-informatik/projects.html>`_
and a custom binary format ``gt`` (see :ref:`sec_gt_format`). 

.. warning::

    The binary format ``gt`` and ``graphml`` are the preferred formats,
    since they are by far the most complete. Both these formats are
    equally complete, but the ``gt`` format is faster and requires less
    storage.

    The ``dot`` and ``gml`` formats are fully supported, but since they
    contain no precise type information, all properties are read as
    strings (or also as double, in the case of ``gml``), and must be
    converted by hand to the desired type. Therefore you should always
    use either ``gt`` or ``graphml``, since they implement an exact
    bit-for-bit representation of all supported :ref:`sec_property_maps`
    types, except when interfacing with other software, or existing
    data, which uses ``dot`` or ``gml``.

A graph can be saved or loaded to a file with the :attr:`~graph_tool.Graph.save`
and :attr:`~graph_tool.Graph.load` methods, which take either a file name or a
file-like object. A graph can also be loaded from disc with the
:func:`~graph_tool.load_graph` function, as such:

.. doctest::

    g = Graph()
    #  ... fill the graph ...
    g.save("my_graph.xml.gz")    
    g2 = load_graph("my_graph.xml.gz")
    # g and g2 should be copies of each other

Graph classes can also be pickled with the :mod:`pickle` module.


An Example: Building a Price Network
------------------------------------

A Price network is the first known model of a "scale-free" graph,
invented in 1976 by `de Solla Price
<http://en.wikipedia.org/wiki/Derek_J._de_Solla_Price>`_. It is defined
dynamically, where at each time step a new vertex is added to the graph,
and connected to an old vertex, with probability proportional to its
in-degree. The following program implements this construction using
``graph-tool``.

.. note::

   Note that it would be much faster just to use the
   :func:`~graph_tool.generation.price_network` function, which is
   implemented in C++, as opposed to the script below which is in pure
   python. The code below is merely a demonstration on how to use the
   library.

.. literalinclude:: price.py
   :linenos:

The following is what should happen when the program is run.

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

Below is the degree distribution, with 100000 nodes. If you want to see
a broader power law, try to increase the number of vertices to something
like :math:`10^6` or :math:`10^7`.

.. figure:: price-deg-dist.*
   :align: center

   In-degree distribution of a price network with :math:`10^5` nodes.


We can draw the graph to see some other features of its topology. For that we
use the :func:`~graph_tool.draw.graph_draw` function.

.. testcode::

   g = load_graph("price.xml.gz")
   age = g.vertex_properties["age"]

   pos = sfdp_layout(g)
   graph_draw(g, pos, output_size=(1000, 1000), vertex_color=[1,1,1,0],
              vertex_fill_color=age, vertex_size=1, edge_pen_width=1.2,
              vcmap=matplotlib.cm.gist_heat_r, output="price.png")

.. figure:: price.*
   :align: center

   A Price network with :math:`10^5` nodes. The vertex colors represent
   the age of the vertex, from older (red) to newer (black).

.. _sec_graph_filtering:

Graph filtering
---------------

One of the very nice features of ``graph-tool`` is the "on-the-fly" filtering of
edges and/or vertices. Filtering means the temporary masking of vertices/edges,
which are in fact not really removed, and can be easily recovered. Vertices or
edges which are to be filtered should be marked with a
:class:`~graph_tool.PropertyMap` with value type ``bool``, and then set with
:meth:`~graph_tool.Graph.set_vertex_filter` or
:meth:`~graph_tool.Graph.set_edge_filter` methods. By default, vertex or edges
with value "1" are `kept` in the graphs, and those with value "0" are filtered
out. This behaviour can be modified with the ``inverted`` parameter of the
respective functions. All manipulation functions and algorithms will work as if
the marked edges or vertices were removed from the graph, with minimum overhead.

.. note::

    It is important to emphasize that the filtering functionality does not add
    any overhead when the graph is not being filtered. In this case, the
    algorithms run just as fast as if the filtering functionality didn't exist.

Here is an example which obtains the minimum spanning tree of a graph,
using the function :func:`~graph_tool.topology.min_spanning_tree` and
edge filtering.

.. testcode::
   :hide:

   from numpy.random import *
   seed(42)

.. testcode::

   g, pos = triangulation(random((500, 2)) * 4, type="delaunay")
   tree = min_spanning_tree(g)
   graph_draw(g, pos=pos, edge_color=tree, output="min_tree.pdf")

.. testcode::
   :hide:

   graph_draw(g, pos=pos, edge_color=tree, output_size=(400, 400),
              output="min_tree.png")


The ``tree`` property map has a bool type, with value "1" if the edge belongs to
the tree, and "0" otherwise. Below is an image of the original graph, with the
marked edges.

.. figure:: min_tree.*
   :align: center

We can now filter out the edges which don't belong to the minimum spanning tree.

.. testcode::

   g.set_edge_filter(tree)
   graph_draw(g, pos=pos, output="min_tree_filtered.pdf")

.. testcode::
   :hide:

   graph_draw(g, pos=pos, output_size=(400, 400), output="min_tree_filtered.png")

This is how the graph looks when filtered:

.. figure:: min_tree_filtered.*
   :align: center

Everything should work transparently on the filtered graph, simply as if the
masked edges were removed. For instance, the following code will calculate the
:func:`~graph_tool.centrality.betweenness` centrality of the edges and vertices,
and draws them as colors and line thickness in the graph.

.. testcode::

    bv, be = betweenness(g)
    be.a /= be.a.max() / 5
    graph_draw(g, pos=pos, vertex_fill_color=bv, edge_pen_width=be,
               output="filtered-bt.pdf")

.. testcode::
   :hide:

   graph_draw(g, pos=pos, vertex_fill_color=bv, edge_pen_width=be,
              output_size=(400, 400), output="filtered-bt.png")

.. figure:: filtered-bt.*
   :align: center


The original graph can be recovered by setting the edge filter to ``None``.

.. testcode::

    g.set_edge_filter(None)
    bv, be = betweenness(g)
    be.a /= be.a.max() / 5
    graph_draw(g, pos=pos, vertex_fill_color=bv, edge_pen_width=be,
               output="nonfiltered-bt.pdf")

.. testcode::
   :hide:

   graph_draw(g, pos=pos, vertex_fill_color=bv, edge_pen_width=be,
              output_size=(400, 400), output="nonfiltered-bt.png")

.. figure:: nonfiltered-bt.*
   :align: center

Everything works in analogous fashion with vertex filtering.

Additionally, the graph can also have its edges reversed with the
:meth:`~graph_tool.Graph.set_reversed` method. This is also an :math:`O(1)`
operation, which does not really modify the graph.

As mentioned previously, the directedness of the graph can also be changed
"on-the-fly" with the :meth:`~graph_tool.Graph.set_directed` method.

.. _sec_graph_views:

Graph views
^^^^^^^^^^^

It is often desired to work with filtered and unfiltered graphs
simultaneously, or to temporarily create a filtered version of graph for
some specific task. For these purposes, graph-tool provides a
:class:`~graph_tool.GraphView` class, which represents a filtered "view"
of a graph, and behaves as an independent graph object, which shares the
underlying data with the original graph. Graph views are constructed by
instantiating a :class:`~graph_tool.GraphView` class, and passing a
graph object which is supposed to be filtered, together with the desired
filter parameters. For example, to create a directed view of the graph
``g`` constructed above, one should do:

.. doctest::

    >>> ug = GraphView(g, directed=True)
    >>> ug.is_directed()
    True

Graph views also provide a much more direct and convenient approach to
vertex/edge filtering: To construct a filtered minimum spanning tree
like in the example above, one must only pass the filter property as the
"efilter" parameter:

.. doctest::

    >>> tv = GraphView(g, efilt=tree)

Note that this is an :math:`O(1)` operation, since it is equivalent (in
speed) to setting the filter in graph ``g`` directly, but in this case
the object ``g`` remains unmodified.

Like above, the result should be the isolated minimum spanning tree:

.. doctest::

    >>> bv, be = betweenness(tv)
    >>> be.a /= be.a.max() / 5
    >>> graph_draw(tv, pos=pos, vertex_fill_color=bv,
    ...            edge_pen_width=be, output="mst-view.pdf")
    <...>

.. testcode::
   :hide:

   graph_draw(tv, pos=pos, vertex_fill_color=bv,
              edge_pen_width=be, output_size=(400, 400),
              output="mst-view.png")

.. figure:: mst-view.*
   :align: center

   A view of the Delaunay graph, isolating only the minimum spanning tree.

.. note::

   :class:`~graph_tool.GraphView` objects behave *exactly* like regular
   :class:`~graph_tool.Graph` objects. In fact,
   :class:`~graph_tool.GraphView` is a subclass of
   :class:`~graph_tool.Graph`. The only difference is that a
   :class:`~graph_tool.GraphView` object shares its internal data with
   its parent :class:`~graph_tool.Graph` class. Therefore, if the
   original :class:`~graph_tool.Graph` object is modified, this
   modification will be reflected immediately in the
   :class:`~graph_tool.GraphView` object, and vice-versa.

For even more convenience, one can supply a function as filter
parameter, which takes a vertex or an edge as single parameter, and
returns ``True`` if the vertex/edge should be kept and ``False``
otherwise. For instance, if we want to keep only the most "central"
edges, we can do:

.. doctest::

    >>> bv, be = betweenness(g)
    >>> u = GraphView(g, efilt=lambda e: be[e] > be.a.max() / 2)

This creates a graph view ``u`` which contains only the edges of ``g``
which have a normalized betweenness centrality larger than half of the
maximum value. Note that, differently from the case above, this is an
:math:`O(E)` operation, where :math:`E` is the number of edges, since
the supplied function must be called :math:`E` times to construct a
filter property map. Thus, supplying a constructed filter map is always
faster, but supplying a function can be more convenient.

The graph view constructed above can be visualized as

.. doctest::

    >>> be.a /= be.a.max() / 5
    >>> graph_draw(u, pos=pos, vertex_fill_color=bv, output="central-edges-view.pdf")
    <...>

.. testcode::
   :hide:

   graph_draw(u, pos=pos, vertex_fill_color=bv, output_size=(400, 400),
              output="central-edges-view.png")


.. figure:: central-edges-view.*
   :align: center

   A view of the Delaunay graph, isolating only the edges with
   normalized betweenness centrality larger than 0.01.

Composing graph views
"""""""""""""""""""""

Since graph views are regular graphs, one can just as easily create
graph views `of graph views`. This provides a convenient way of
composing filters. For instance, in order to isolate the minimum
spanning tree of all vertices of the example above which have a degree
larger than four, one can do:


    >>> u = GraphView(g, vfilt=lambda v: v.out_degree() > 4)
    >>> tree = min_spanning_tree(u)
    >>> u = GraphView(u, efilt=tree)

The resulting graph view can be visualized as

.. doctest::

    >>> graph_draw(u, pos=pos, output="composed-filter.pdf")
    <...>

.. testcode::
   :hide:

   graph_draw(u, pos=pos, output_size=(400, 400), output="composed-filter.png")

.. figure:: composed-filter.*
   :align: center

   A composed view, obtained as the minimum spanning tree of all vertices
   in the graph which have a degree larger than four.