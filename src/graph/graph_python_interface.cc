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

struct add_new_vertex
{
    template <class Graph>
        void operator()(Graph& g, python::object& new_v) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        vertex_t v = add_vertex(g);
        new_v = python::object(PythonVertex<Graph>(g, v));
    }
};

python::object GraphInterface::AddVertex()
{
    python::object new_v;
    check_filter(*this, lambda::bind<void>(add_new_vertex(), lambda::_1,
                                           lambda::var(new_v)),
                 reverse_check(), directed_check());
    return new_v;
}

struct get_vertex_descriptor
{
    template <class Graph>
    void operator()(Graph& g, python::object& v,
                    typename graph_traits<Graph>::vertex_descriptor& vertex)
        const
    {
        PythonVertex<Graph>& pv = python::extract<PythonVertex<Graph>&>(v);
        vertex = pv.GetDescriptor();
    }
};

void GraphInterface::RemoveVertex(python::object v)
{
    graph_traits<multigraph_t>::vertex_descriptor dv;
    check_filter(*this, lambda::bind<void>(get_vertex_descriptor(), lambda::_1,
                                           lambda::var(v), lambda::var(dv)),
                 reverse_check(), directed_check());

    //shift properties
    size_t N = num_vertices(_mg);
    for (size_t i = _vertex_index[dv]; i < N-1; ++i)
    {
        graph_traits<multigraph_t>::vertex_descriptor v = vertex(i, _mg);
        for (typeof(_properties.begin()) p = _properties.begin();
             p != _properties.end(); ++p)
            if (p->second->key() == typeid(vertex_t))
                try
                {
                    p->second->put(v, p->second->get(vertex(i+1, _mg)));
                }
                catch (dynamic_const_put_error) {} // index prop. is const

    }

    //remove vertex
    if (in_degree(dv, _mg) + out_degree(dv, _mg) > 0)
    {
        clear_vertex(dv, _mg);
        ReIndexEdges();
    }
    remove_vertex(dv, _mg);
}

struct add_new_edge
{
    template <class Graph>
        void operator()(Graph& g, const python::object& s,
                        const python::object& t, python::object& new_e) const
    {
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;

        PythonVertex<Graph>& src = python::extract<PythonVertex<Graph>&>(s);
        PythonVertex<Graph>& tgt = python::extract<PythonVertex<Graph>&>(t);
        edge_t e = add_edge(src.GetDescriptor(), tgt.GetDescriptor(), g).first;
        new_e = python::object(PythonEdge<Graph>(g, e));
    }
};

python::object GraphInterface::AddEdge(python::object s, python::object t)
{
    python::object new_e;
    check_filter(*this, lambda::bind<void>(add_new_edge(), lambda::_1,
                                           lambda::var(s), lambda::var(t),
                                           lambda::var(new_e)),
                 reverse_check(), directed_check());
    return new_e;
}

struct get_edge_descriptor
{
    template <class Graph>
    void operator()(Graph& g, python::object& v,
                    typename graph_traits
                    <GraphInterface::multigraph_t>::edge_descriptor& edge)
        const
    {
        PythonEdge<Graph>& pv = python::extract<PythonEdge<Graph>&>(v);
        edge = pv.GetDescriptor();
    }
};

void GraphInterface::RemoveEdge(python::object e)
{
    graph_traits<multigraph_t>::edge_descriptor de;
    check_filter(*this, lambda::bind<void>(get_edge_descriptor(), lambda::_1,
                                           lambda::var(e), lambda::var(de)),
                 reverse_check(), directed_check());
    remove_edge(de, _mg);
    ReIndexEdges();
}

struct get_property_map
{
    get_property_map(const string& name, dynamic_property_map& dp,
                     python::object& pmap)
        : _name(name), _dp(dp), _pmap(pmap) {}

    template <class ValueType>
    void operator()(ValueType) const
    {
        if (typeid(ValueType) == _dp.value())
            _pmap = python::object(PythonPropertyMap<ValueType>(_name, _dp));
    }

    const string& _name;
    dynamic_property_map& _dp;
    python::object& _pmap;
};


python::dict
GraphInterface::GetVertexProperties() const
{
    typedef graph_traits<multigraph_t>::vertex_descriptor vertex_t;;

    python::dict props;
    for(typeof(_properties.begin()) iter = _properties.begin();
        iter != _properties.end(); ++iter)
        if (iter->second->key() == typeid(vertex_t))
        {
            python::object pmap;
            mpl::for_each<value_types>(get_property_map
                                       (iter->first, *iter->second, pmap));
            props[iter->first] = pmap;
        }
    return props;
}

python::dict
GraphInterface::GetEdgeProperties() const
{
    typedef graph_traits<multigraph_t>::edge_descriptor edge_t;;

    python::dict props;
    for(typeof(_properties.begin()) iter = _properties.begin();
        iter != _properties.end(); ++iter)
        if (iter->second->key() == typeid(edge_t))
        {
            python::object pmap;
            mpl::for_each<value_types>(get_property_map
                                       (iter->first, *iter->second, pmap));
            props[iter->first] = pmap;
        }
    return props;
}

python::dict
GraphInterface::GetGraphProperties() const
{
    python::dict props;
    for(typeof(_properties.begin()) iter = _properties.begin();
        iter != _properties.end(); ++iter)
        if (iter->second->key() == typeid(graph_property_tag))
        {
            python::object pmap;
            mpl::for_each<value_types>(get_property_map
                                       (iter->first, *iter->second, pmap));
            props[iter->first] = pmap;
        }
    return props;
}


//
// Below are the functions with will properly register all the types to python,
// for every filter, type, etc.
//

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
            .def("__str__", &PythonVertex<Graph>::GetString)
            .def("__hash__", &PythonVertex<Graph>::GetHash);

        class_<PythonEdge<Graph> > ("Edge", no_init)
            .def("source", &PythonEdge<Graph>::GetSource)
            .def("target", &PythonEdge<Graph>::GetTarget)
            .def(python::self == python::self)
            .def(python::self != python::self)
            .def("__str__", &PythonEdge<Graph>::GetString)
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

// this will register the property maps types and all its possible access
// functions to python
struct export_property_map
{
    export_property_map(const string& name, const GraphInterface &gi)
        : _name(name), _gi(gi) {}

    template <class ValueType>
    struct export_access
    {
        typedef PythonPropertyMap<ValueType> pmap_t;

        export_access(python::class_<pmap_t>& pclass)
            : _pclass(pclass) {}

        template <class Graph>
        void operator()(const Graph&) const
        {
            _pclass
                .def("__getitem__",
                     &pmap_t::template GetValue<PythonVertex<Graph> >)
                .def("__setitem__",
                     &pmap_t::template SetValue<PythonVertex<Graph> >)
                .def("__getitem__",
                     &pmap_t::template GetValue<PythonEdge<Graph> >)
                .def("__setitem__",
                     &pmap_t::template SetValue<PythonEdge<Graph> >)
                .def("__getitem__",
                     &pmap_t::template GetValue<GraphInterface>)
                .def("__setitem__",
                     &pmap_t::template SetValue<GraphInterface>);
        }

        python::class_<PythonPropertyMap<ValueType> >& _pclass;
    };

    template <class ValueType>
    void operator()(ValueType) const
    {
        typedef PythonPropertyMap<ValueType> pmap_t;
        python::class_<pmap_t> pclass(_name.c_str(), python::no_init);
        pclass.def("__hash__", &pmap_t::GetHash);
        pclass.def("get_type", &pmap_t::GetType);
        check_filter(_gi, lambda::bind<void>(export_access<ValueType>(pclass),
                                            lambda::_1),
                     reverse_check(), directed_check(), true);
    }

    string _name;
    const GraphInterface& _gi;
};

// register everything

void GraphInterface::ExportPythonInterface() const
{
    check_filter(*this, lambda::bind<void>(export_python_interface(),
                                           lambda::_1),
                 reverse_check(), directed_check(), true);
    mpl::for_each<value_types>(export_property_map("PropertyMap", *this));
}
