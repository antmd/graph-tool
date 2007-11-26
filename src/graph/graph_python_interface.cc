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

#include <boost/python.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include "graph_python_interface.hh"
#include "graph.hh"
#include "graph_filtering.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;


struct get_vertex_iterator
{
    template <class Graph>
    void operator()(const Graph& g,
                    python::object& iter) const
    {
        typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
        iter = python::object(PythonIterator<Graph,
                                             PythonVertex<Graph>,
                                             vertex_iterator>(g, vertices(g)));
    }
};

python::object
GraphInterface::Vertices() const
{
    python::object iter;
    check_filter(*this, lambda::bind<void>(get_vertex_iterator(),
                                           lambda::_1,
                                           lambda::var(iter)),
                 reverse_check(), directed_check());
    return iter;
}

struct get_edge_iterator
{
    template <class Graph>
    void operator()(const Graph& g,
                    python::object& iter) const
    {
        typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
        iter = python::object(PythonIterator<Graph,
                                             PythonEdge<Graph>,
                                             edge_iterator>(g, edges(g)));
    }
};

python::object
GraphInterface::Edges() const
{
    python::object iter;
    check_filter(*this, lambda::bind<void>(get_edge_iterator(),
                                           lambda::_1,
                                           lambda::var(iter)),
                 reverse_check(), directed_check());
    return iter;
}

// this will register all the Vertex/Edge classes to python
struct export_python_interface
{
    template <class Graph>
    void operator()(const Graph&) const
    {
        using namespace boost::python;

        class_<PythonVertex<Graph> >("Vertex", no_init)
            .def("in_degree", &PythonVertex<Graph>::GetInDegree)
            .def("out_degree", &PythonVertex<Graph>::GetOutDegree)
            .def("out_edges", &PythonVertex<Graph>::OutEdges)
            .def("in_edges", &PythonVertex<Graph>::InEdges)
            .def(python::self == python::self)
            .def(python::self != python::self)
            .def("__hash__", &PythonVertex<Graph>::GetHash);

        class_<PythonEdge<Graph> > ("Edge", no_init)
            .def("source", &PythonEdge<Graph>::GetSource)
            .def("target", &PythonEdge<Graph>::GetTarget)
            .def(python::self == python::self)
            .def(python::self != python::self)
            .def("__hash__", &PythonEdge<Graph>::GetHash);

        typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
        class_<PythonIterator<Graph,PythonVertex<Graph>,
                              vertex_iterator> >("VertexIterator",no_init)
            .def("__iter__", objects::identity_function())
            .def("next", &PythonIterator<Graph,PythonVertex<Graph>,
                                         vertex_iterator>::Next);

        typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
        class_<PythonIterator<Graph,PythonEdge<Graph>,
                              edge_iterator> >("EdgeIterator", no_init)
            .def("__iter__", objects::identity_function())
            .def("next", &PythonIterator<Graph,PythonEdge<Graph>,
                                         edge_iterator>::Next);

        typedef typename graph_traits<Graph>::out_edge_iterator
            out_edge_iterator;
        class_<PythonIterator<Graph,PythonEdge<Graph>,
                              out_edge_iterator> >("OutEdgeIterator", no_init)
            .def("__iter__", objects::identity_function())
            .def("next", &PythonIterator<Graph,PythonEdge<Graph>,
                                         out_edge_iterator>::Next);


        typedef typename graph_traits<Graph>::directed_category
            directed_category;
        typedef typename is_convertible<directed_category,
                                        boost::directed_tag>::type is_directed;
        if (is_directed::value)
        {
            typedef typename in_edge_iteratorS<Graph>::type in_edge_iterator;
            class_<PythonIterator<Graph,PythonEdge<Graph>,
                                  in_edge_iterator> >("InEdgeIterator", no_init)
                .def("__iter__", objects::identity_function())
                .def("next", &PythonIterator<Graph,PythonEdge<Graph>,
                                             in_edge_iterator>::Next);
        }
    }
};

void GraphInterface::ExportPythonInterface() const
{
    check_filter(*this, lambda::bind<void>(export_python_interface(),
                                           lambda::_1),
                 reverse_check(), directed_check(), true);
}
