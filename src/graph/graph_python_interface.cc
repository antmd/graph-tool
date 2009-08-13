// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007  Tiago de Paula Peixoto <tiago@forked.de>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "graph_filtering.hh"
#include "graph.hh"
#include "graph_python_interface.hh"

#include <boost/python.hpp>
#include <boost/lambda/bind.hpp>
#include <set>

using namespace std;
using namespace boost;
using namespace graph_tool;

struct get_vertex_iterator
{
    template <class Graph>
    void operator()(Graph& g, GraphInterface& gi,
                    python::object& iter) const
    {
        typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
        iter = python::object(PythonIterator<PythonVertex,
                                             vertex_iterator>(gi, vertices(g)));
    }
};

python::object
GraphInterface::Vertices() const
{
    python::object iter;
    run_action<>()(const_cast<GraphInterface&>(*this),
                   lambda::bind<void>(get_vertex_iterator(), lambda::_1,
                                      lambda::var(const_cast<GraphInterface&>(*this)),
                                      lambda::var(iter)))();
    return iter;
}

struct get_vertex_soft
{
    template <class Graph>
    void operator()(Graph& g,
                    GraphInterface& gi,
                    size_t i, python::object& v) const
    {
        v = python::object(PythonVertex(gi, vertex(i, g)));
    }
};

struct get_vertex_hard
{
    template <class Graph>
    void operator()(Graph& g, GraphInterface& gi, size_t i,
                    python::object& v) const
    {
        size_t c = 0;
        typename graph_traits<Graph>::vertex_iterator vi, v_end;
        for (tie(vi, v_end) = vertices(g); vi != v_end; ++vi)
        {
            if (c == i)
            {
                v = python::object(PythonVertex(gi, *vi));
                return;
            }
            ++c;
        }
        v = python::object(PythonVertex(gi,
                                        graph_traits<Graph>::null_vertex()));
    }
};

python::object
GraphInterface::Vertex(size_t i) const
{
    python::object v;
    if (IsVertexFilterActive())
        run_action<>()(const_cast<GraphInterface&>(*this),
                       lambda::bind<void>(get_vertex_hard(), lambda::_1,
                                          lambda::var(const_cast<GraphInterface&>(*this)),
                                          i, lambda::var(v)))();
    else
        run_action<>()(const_cast<GraphInterface&>(*this),
                       lambda::bind<void>(get_vertex_soft(), lambda::_1,
                                          lambda::var(const_cast<GraphInterface&>(*this)), i,
                                          lambda::var(v)))();
    return v;
}

struct get_edge_iterator
{
    template <class Graph>
    void operator()(Graph& g, GraphInterface& gi,
                    python::object& iter) const
    {
        typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
        iter = python::object(PythonIterator<PythonEdge<Graph>,
                                             edge_iterator>(gi, edges(g)));
    }
};

python::object
GraphInterface::Edges() const
{
    python::object iter;
    run_action<>()(const_cast<GraphInterface&>(*this),
                   lambda::bind<void>(get_edge_iterator(), lambda::_1,
                                      lambda::var(const_cast<GraphInterface&>(*this)),
                                      lambda::var(iter)))();
    return iter;
}

python::object GraphInterface::AddVertex()
{
    return python::object(PythonVertex(*this, add_vertex(_mg)));
}

struct shift_vertex_property
{
    template <class Graph, class PropertyMap>
    void operator()(const Graph& g, size_t vi, PropertyMap pmap)
        const
    {
        size_t N = num_vertices(g);
        if (N <= 1 || vi >= N - 1)
            return;
        for (size_t i = vi; i < N-1; ++i)
        {
            pmap[vertex(i, g)] = pmap[vertex(i+1, g)];
        }
    }
};

void GraphInterface::RemoveVertex(const python::object& v)
{
    PythonVertex& pv = python::extract<PythonVertex&>(v);
    pv.CheckValid();
    vertex_t dv = pv.GetDescriptor();
    pv.SetValid(false);

    //remove vertex
    clear_vertex(dv, _mg);
    remove_vertex(dv, _mg);
}

struct add_new_edge
{
    template <class Graph, class EdgeIndexMap>
    void operator()(Graph& g, GraphInterface& gi, const PythonVertex& s,
                    const PythonVertex& t, EdgeIndexMap edge_index,
                    python::object& new_e) const
    {
        typename graph_traits<Graph>::edge_descriptor e =
            add_edge(s.GetDescriptor(), t.GetDescriptor(), g).first;
        new_e = python::object(PythonEdge<Graph>(gi, e));
        gi.AddEdgeIndex(e);
    }
};

void GraphInterface::AddEdgeIndex(const edge_t& e)
{
    if (!_free_indexes.empty())
    {
        _edge_index[e] = _free_indexes.front();
        _free_indexes.pop_front();
    }
    else
    {
        _edge_index[e] = _nedges;
        _max_edge_index = _nedges;
    }
    _nedges++;
}

python::object GraphInterface::AddEdge(const python::object& s,
                                       const python::object& t)
{
    PythonVertex& src = python::extract<PythonVertex&>(s);
    PythonVertex& tgt = python::extract<PythonVertex&>(t);
    src.CheckValid();
    tgt.CheckValid();
    python::object new_e;
    run_action<>()(*this, lambda::bind<void>(add_new_edge(), lambda::_1,
                                             lambda::var(*this), src, tgt,
                                             _edge_index,
                                             lambda::var(new_e)))();
    return new_e;
}

struct get_edge_descriptor
{
    template <class Graph>
    void operator()(const Graph& g, const python::object& e,
                    typename GraphInterface::edge_t& edge,
                    bool& found)
        const
    {
        PythonEdge<Graph>& pe = python::extract<PythonEdge<Graph>&>(e);
        pe.CheckValid();
        edge = pe.GetDescriptor();
        pe.SetValid(false);
        found = true;
    }
};

void GraphInterface::RemoveEdge(const python::object& e)
{
    edge_t de;
    bool found = false;
    run_action<>()(*this,
                   lambda::bind<void>(get_edge_descriptor(), lambda::_1,
                                      lambda::var(e), lambda::var(de),
                                      lambda::var(found)))();
    if (!found)
        throw ValueException("invalid edge descriptor");
    RemoveEdgeIndex(de);
}

void GraphInterface::RemoveEdgeIndex(const edge_t& e)
{
    size_t index = _edge_index[e];
    if (index == _max_edge_index)
    {
        if (_max_edge_index > 0)
            _max_edge_index--;

        while (!_free_indexes.empty() &&
               _max_edge_index == _free_indexes.back())
        {
            _free_indexes.pop_back();
            if (_max_edge_index > 0)
                _max_edge_index--;
        }
    }
    else
    {
        typeof(_free_indexes.begin()) pos =
                lower_bound(_free_indexes.begin(), _free_indexes.end(), index);
        _free_indexes.insert(pos, index);
    }
    _nedges--;
    remove_edge(e, _mg);
}

struct get_degree_map
{
    template <class Graph, class DegreeMap, class DegS>
    void operator()(const Graph& g, DegreeMap deg_map, DegS deg) const
    {
        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            deg_map[v] = deg(v, g);
        }
    }
};

python::object GraphInterface::DegreeMap(string deg) const
{
    typedef property_map_type::apply<double,
                                     GraphInterface::vertex_index_map_t>::type
        map_t;

    map_t deg_map(_vertex_index);
    deg_map.reserve(num_vertices(_mg));

    if (deg == "in")
        run_action<>()(const_cast<GraphInterface&>(*this),
                       lambda::bind<void>(get_degree_map(), lambda::_1,
                                          deg_map, in_degreeS()))();
    else if (deg == "out")
        run_action<>()(const_cast<GraphInterface&>(*this),
                       lambda::bind<void>(get_degree_map(), lambda::_1,
                                          deg_map, out_degreeS()))();
    else if (deg == "total")
        run_action<>()(const_cast<GraphInterface&>(*this),
                       lambda::bind<void>(get_degree_map(), lambda::_1,
                                          deg_map, total_degreeS()))();
    return python::object(PythonPropertyMap<map_t>(deg_map));
}

//
// Below are the functions with will properly register all the types to python,
// for every filter, type, etc.
//

// this will register all the Vertex/Edge classes to python
struct export_python_interface
{
    template <class Graph>
    void operator()(const Graph*, set<string>& v_iterators) const
    {
        using namespace boost::python;

        class_<PythonEdge<Graph> >
            ("Edge", "This class represents an edge in a graph", no_init)
            .def("source", &PythonEdge<Graph>::GetSource,
                 "Return the source vertex")
            .def("target", &PythonEdge<Graph>::GetTarget,
                 "Return the target vertex")
            .def("is_valid", &PythonEdge<Graph>::IsValid,
                 "Return whether the edge is valid")
            .def(python::self == python::self)
            .def(python::self != python::self)
            .def("__str__", &PythonEdge<Graph>::GetString)
            .def("__hash__", &PythonEdge<Graph>::GetHash);

        typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
        if (v_iterators.find(typeid(vertex_iterator).name()) ==
            v_iterators.end())
        {
            class_<PythonIterator<PythonVertex, vertex_iterator> >
                ("VertexIterator", no_init)
                .def("__iter__", objects::identity_function())
                .def("next", &PythonIterator<PythonVertex,
                                             vertex_iterator>::Next);
            v_iterators.insert(typeid(vertex_iterator).name());
        }

        typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
        class_<PythonIterator<PythonEdge<Graph>,
                              edge_iterator> >("EdgeIterator", no_init)
            .def("__iter__", objects::identity_function())
            .def("next", &PythonIterator<PythonEdge<Graph>,
                                         edge_iterator>::Next);

        typedef typename graph_traits<Graph>::out_edge_iterator
            out_edge_iterator;
        class_<PythonIterator<PythonEdge<Graph>,
                              out_edge_iterator> >("OutEdgeIterator", no_init)
            .def("__iter__", objects::identity_function())
            .def("next", &PythonIterator<PythonEdge<Graph>,
                                         out_edge_iterator>::Next);

        typedef typename graph_traits<Graph>::directed_category
            directed_category;
        typedef typename is_convertible<directed_category,
                                        boost::directed_tag>::type is_directed;
        if (is_directed::value)
        {
            typedef typename in_edge_iteratorS<Graph>::type in_edge_iterator;
            class_<PythonIterator<PythonEdge<Graph>,
                                  in_edge_iterator> >("InEdgeIterator", no_init)
                .def("__iter__", objects::identity_function())
                .def("next", &PythonIterator<PythonEdge<Graph>,
                                             in_edge_iterator>::Next);
        }
    }
};

void export_python_properties(GraphInterface&);

PythonPropertyMap<GraphInterface::vertex_index_map_t>
get_vertex_index(GraphInterface& g)
{
    return PythonPropertyMap<GraphInterface::vertex_index_map_t>
        (g.GetVertexIndex());
}

PythonPropertyMap<GraphInterface::edge_index_map_t>
get_edge_index(GraphInterface& g)
{
    return PythonPropertyMap<GraphInterface::edge_index_map_t>
        (g.GetEdgeIndex());
}

// register everything

void GraphInterface::ExportPythonInterface() const
{
    using namespace boost::python;

    class_<PythonVertex>
        ("Vertex", "This class represents a vertex in a graph", no_init)
        .def("in_degree", &PythonVertex::GetInDegree,
             "Return the in-degree")
        .def("out_degree", &PythonVertex::GetOutDegree,
             "Return the out-degree")
        .def("out_edges", &PythonVertex::OutEdges,
             "Return an iterator over the out-edges")
        .def("in_edges", &PythonVertex::InEdges,
             "Return an iterator over the in-edges")
        .def("is_valid", &PythonVertex::IsValid,
             "Return whether the vertex is valid")
        .def(python::self == python::self)
        .def(python::self != python::self)
        .def("__str__", &PythonVertex::GetString)
        .def("__int__", &PythonVertex::GetIndex)
        .def("__hash__", &PythonVertex::GetHash);

    set<string> v_iterators;
    typedef mpl::transform<graph_tool::detail::all_graph_views,
                           mpl::quote1<add_pointer> >::type graph_views;
    mpl::for_each<graph_views>(lambda::bind<void>(export_python_interface(),
                                                  lambda::_1,
                                                  lambda::var(v_iterators)));
    export_python_properties(const_cast<GraphInterface&>(*this));
    def("new_vertex_property",
        &new_property<GraphInterface::vertex_index_map_t>);
    def("new_edge_property", &new_property<GraphInterface::edge_index_map_t>);
    def("new_graph_property",
        &new_property<ConstantPropertyMap<size_t,graph_property_tag> >);

    def("get_vertex_index", get_vertex_index);
    def("get_edge_index", get_edge_index);
}
