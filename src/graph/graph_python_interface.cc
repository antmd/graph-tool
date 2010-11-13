// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007-2010 Tiago de Paula Peixoto <tiago@skewed.de>
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
#include <set>

using namespace std;
using namespace boost;
using namespace graph_tool;

namespace graph_tool
{

struct get_vertex_iterator
{
    template <class Graph>
    void operator()(Graph& g, python::object& pg,
                    python::object& iter) const
    {
        typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
        iter = python::object(PythonIterator<PythonVertex,
                                             vertex_iterator>(pg, vertices(g)));
    }
};

python::object get_vertices(python::object g)
{
    GraphInterface& gi = python::extract<GraphInterface&>(g());
    python::object iter;
    run_action<>()(gi, bind<void>(get_vertex_iterator(), _1,
                                  ref(g),
                                  ref(iter)))();
    return iter;
}

struct get_vertex_soft
{
    template <class Graph>
    void operator()(Graph& g,
                    python::object& pg,
                    size_t i, python::object& v) const
    {
        v = python::object(PythonVertex(pg, vertex(i, g)));
    }
};

struct get_vertex_hard
{
    template <class Graph>
    void operator()(Graph& g, python::object& pg, size_t i,
                    python::object& v) const
    {
        size_t c = 0;
        typename graph_traits<Graph>::vertex_iterator vi, v_end;
        for (tie(vi, v_end) = vertices(g); vi != v_end; ++vi)
        {
            if (c == i)
            {
                v = python::object(PythonVertex(pg, *vi));
                return;
            }
            ++c;
        }
        v = python::object(PythonVertex(pg,
                                        graph_traits<Graph>::null_vertex()));
    }
};

python::object get_vertex(python::object g, size_t i)
{
    GraphInterface& gi = python::extract<GraphInterface&>(g());
    python::object v;
    if (gi.IsVertexFilterActive())
        run_action<>()(gi,
                       bind<void>(get_vertex_hard(), _1, ref(g), i, ref(v)))();
    else
        run_action<>()(gi,
                       bind<void>(get_vertex_soft(), _1, ref(g), i, ref(v)))();
    return v;
}

struct get_edge_iterator
{
    template <class Graph>
    void operator()(Graph& g, const python::object& pg, python::object& iter)
        const
    {
        typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
        iter = python::object(PythonIterator<PythonEdge<Graph>,
                                             edge_iterator>(pg, edges(g)));
    }
};

python::object get_edges(python::object g)
{
    GraphInterface& gi = python::extract<GraphInterface&>(g());
    python::object iter;
    run_action<>()(gi,
                   bind<void>(get_edge_iterator(), _1, ref(g), ref(iter)))();
    return iter;
}

python::object add_vertex(python::object g)
{
    GraphInterface& gi = python::extract<GraphInterface&>(g());
    return python::object(PythonVertex(g, add_vertex(gi.GetGraph())));
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

void remove_vertex(GraphInterface& gi, const python::object& v)
{
    PythonVertex& pv = python::extract<PythonVertex&>(v);
    pv.CheckValid();
    GraphInterface::vertex_t dv = pv.GetDescriptor();
    pv.SetValid(false);

    //remove vertex
    clear_vertex(dv, gi.GetGraph());
    remove_vertex(dv, gi.GetGraph());
}

struct add_new_edge
{
    template <class Graph, class EdgeIndexMap>
    void operator()(Graph& g, python::object& pg, GraphInterface& gi,
                    const PythonVertex& s, const PythonVertex& t,
                    EdgeIndexMap edge_index, python::object& new_e) const
    {
        typename graph_traits<Graph>::edge_descriptor e =
            add_edge(s.GetDescriptor(), t.GetDescriptor(), g).first;
        new_e = python::object(PythonEdge<Graph>(pg, e));
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

python::object add_edge(python::object g, const python::object& s,
                        const python::object& t)
{
    PythonVertex& src = python::extract<PythonVertex&>(s);
    PythonVertex& tgt = python::extract<PythonVertex&>(t);
    src.CheckValid();
    tgt.CheckValid();
    GraphInterface& gi = python::extract<GraphInterface&>(g());
    python::object new_e;
    run_action<>()(gi, bind<void>(add_new_edge(), _1, ref(g), ref(gi), src, tgt,
                                  gi.GetEdgeIndex(), ref(new_e)))();
    return new_e;
}

struct get_edge_descriptor
{
    template <class Graph>
    void operator()(const Graph& g, const python::object& e,
                    typename GraphInterface::edge_t& edge,
                    bool& found)  const
    {
        PythonEdge<Graph>& pe = python::extract<PythonEdge<Graph>&>(e);
        pe.CheckValid();
        pe.SetValid(false);
        typename graph_traits<Graph>::out_edge_iterator e_begin, e_end;
        tie(e_begin, e_end) = out_edges(source(pe.GetDescriptor(),g),g);
        while(e_begin != e_end && *e_begin != pe.GetDescriptor())
            ++e_begin;
        if (e_begin == e_end)
            return; // invalid edge descriptor
        edge = pe.GetDescriptor();
        found = true;
    }
};

void remove_edge(GraphInterface& gi, const python::object& e)
{
    GraphInterface::edge_t de;
    bool found = false;
    run_action<>()(gi, bind<void>(get_edge_descriptor(), _1, ref(e), ref(de),
                                  ref(found)))();
    if (!found)
        throw ValueException("invalid edge descriptor");
    gi.RemoveEdgeIndex(de);
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
                       bind<void>(get_degree_map(), _1,
                                  deg_map, in_degreeS()))();
    else if (deg == "out")
        run_action<>()(const_cast<GraphInterface&>(*this),
                       bind<void>(get_degree_map(), _1,
                                  deg_map, out_degreeS()))();
    else if (deg == "total")
        run_action<>()(const_cast<GraphInterface&>(*this),
                       bind<void>(get_degree_map(), _1,
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
            ("Edge", no_init)
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

} // namespace graph_tool

// register everything

void export_python_properties();

void export_python_interface()
{
    using namespace boost::python;

    class_<PythonVertex>
        ("Vertex", no_init)
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
    mpl::for_each<graph_views>(bind<void>(graph_tool::export_python_interface(),
                                          _1,ref(v_iterators)));
    export_python_properties();
    def("new_vertex_property",
        &new_property<GraphInterface::vertex_index_map_t>);
    def("new_edge_property", &new_property<GraphInterface::edge_index_map_t>);
    def("new_graph_property",
        &new_property<ConstantPropertyMap<size_t,graph_property_tag> >);

    def("get_vertex", get_vertex);
    def("get_vertices", get_vertices);
    def("get_edges", get_edges);
    def("add_vertex", graph_tool::add_vertex);
    def("add_edge", graph_tool::add_edge);
    def("remove_vertex", graph_tool::remove_vertex);
    def("remove_edge", graph_tool::remove_edge);

    def("get_vertex_index", get_vertex_index);
    def("get_edge_index", get_edge_index);
}
