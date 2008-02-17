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

#include "graph.hh"
#include "graph_python_interface.hh"

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

using namespace std;
using namespace graph_tool;
using namespace boost;
using namespace boost::python;

struct LibInfo
{
    string GetName()      const {return PACKAGE_NAME;}
    string GetAuthor()    const {return AUTHOR;}
    string GetCopyright() const {return COPYRIGHT;}
    string GetVersion()   const {return VERSION " (commit " GIT_COMMIT
                                        ", " GIT_COMMIT_DATE ")";}
    string GetLicense()   const {return "GPL version 3 or above";}
    string GetCXXFLAGS()  const {return CXXFLAGS;}
    string GetInstallPrefix() const {return INSTALL_PREFIX;}
    string GetPythonDir() const {return PYTHON_DIR;}
};

template <class ValueType>
struct vector_from_list
{
    vector_from_list()
    {
        converter::registry::push_back
            (&convertible, &construct,
             boost::python::type_id<vector<ValueType> >());
    }

    static void* convertible(PyObject* obj_ptr)
    {
        handle<> x(borrowed(obj_ptr));
        object o(x);
        size_t N = len(o);
        for (size_t i = 0; i < N; ++i)
        {
            extract<ValueType> elem(o[i]);
            if (!elem.check())
                return 0;
        }
        return obj_ptr;
    }

    static void construct(PyObject* obj_ptr,
                          converter::rvalue_from_python_stage1_data* data)
    {
        handle<> x(borrowed(obj_ptr));
        object o(x);
        vector<ValueType> value;
        size_t N = len(o);
        for (size_t i = 0; i < N; ++i)
            value.push_back(extract<ValueType>(o[i]));
        void* storage =
            ( (boost::python::converter::rvalue_from_python_storage
               <vector<ValueType> >*) data)->storage.bytes;
        new (storage) vector<ValueType>(value);
        data->convertible = storage;
    }
};

struct export_vector_types
{
    template <class ValueType>
    void operator()(ValueType) const
    {
        string type_name = get_type_name<>()(typeid(ValueType));
        if (type_name == "long double")
            type_name = "long_double";
        string name = "Vector_" + type_name;
        class_<vector<ValueType> >(name.c_str())
            .def(vector_indexing_suite<vector<ValueType> >());
        vector_from_list<ValueType>();
    }
};

// exception translation
static PyObject* pyex =
    PyErr_NewException((char *) "libgraph_tool_core.GraphError",
                       PyExc_Exception, NULL);

void graph_exception_translator(const GraphException& e)
{
    PyObject* message = PyString_FromString(e.what());
    PyObject_SetAttrString(pyex, "message", message);
    PyErr_SetString(pyex, e.what());
}

template <class Exception>
void translate(const Exception& e)
{
    PyErr_SetString(PyExc_RuntimeError, e.what());
}

void raise_error(const string& msg)
{
    throw GraphException(msg);
}

template <class T1, class T2>
struct pair_to_tuple
{
    static PyObject* convert(const pair<T1,T2>& p)
    {
        boost::python::tuple t = boost::python::make_tuple(p.first,p.second);
        return incref(t.ptr());
    }
};

BOOST_PYTHON_MODULE(libgraph_tool_core)
{
    GraphInterface().ExportPythonInterface();

    PyModule_AddObject(python::detail::current_scope, "GraphError", pyex);
    register_exception_translator<GraphException>(graph_exception_translator);

    def("raise_error", &raise_error);

    mpl::for_each<mpl::push_back<scalar_types,string>::type>(export_vector_types());

    class_<GraphInterface>("GraphInterface", init<>())
        .def(init<GraphInterface>())
        .def("GetNumberOfVertices", &GraphInterface::GetNumberOfVertices)
        .def("GetNumberOfEdges", &GraphInterface::GetNumberOfEdges)
        .def("SetDirected", &GraphInterface::SetDirected)
        .def("GetDirected", &GraphInterface::GetDirected)
        .def("SetReversed", &GraphInterface::SetReversed)
        .def("GetReversed", &GraphInterface::GetReversed)
        .def("SetVertexFilterProperty",
             &GraphInterface::SetVertexFilterProperty)
        .def("GetVertexFilterProperty",
             &GraphInterface::GetVertexFilterProperty)
        .def("IsVertexFilterActive", &GraphInterface::IsVertexFilterActive)
        .def("SetEdgeFilterProperty",
             &GraphInterface::SetEdgeFilterProperty)
        .def("GetEdgeFilterProperty",
             &GraphInterface::GetEdgeFilterProperty)
        .def("IsEdgeFilterActive", &GraphInterface::IsEdgeFilterActive)
        .def("AddEdgeProperty",  &GraphInterface::AddEdgeProperty)
        .def("AddVertexProperty",  &GraphInterface::AddVertexProperty)
        .def("AddGraphProperty",  &GraphInterface::AddGraphProperty)
        .def("RemoveEdgeProperty",  &GraphInterface::RemoveEdgeProperty)
        .def("RemoveVertexProperty",  &GraphInterface::RemoveVertexProperty)
        .def("RemoveGraphProperty",  &GraphInterface::RemoveGraphProperty)
        .def("PurgeVertices",  &GraphInterface::PurgeVertices)
        .def("PurgeEdges",  &GraphInterface::PurgeEdges)
        .def("ReIndexEdges",  &GraphInterface::ReIndexEdges)
        .def("InsertEdgeIndexProperty",
             &GraphInterface::InsertEdgeIndexProperty)
        .def("InsertVertexIndexProperty",
             &GraphInterface::InsertVertexIndexProperty)
        .def("WriteToFile", &GraphInterface::WriteToFile)
        .def("ReadFromFile",&GraphInterface::ReadFromFile)
        .def("Vertices", &GraphInterface::Vertices)
        .def("Vertex", &GraphInterface::Vertex)
        .def("Edges", &GraphInterface::Edges)
        .def("AddVertex", &GraphInterface::AddVertex)
        .def("AddEdge", &GraphInterface::AddEdge)
        .def("RemoveVertex", &GraphInterface::RemoveVertex)
        .def("RemoveEdge", &GraphInterface::RemoveEdge)
        .def("Clear", &GraphInterface::Clear)
        .def("GetVertexProperties", &GraphInterface::GetVertexProperties)
        .def("GetEdgeProperties", &GraphInterface::GetEdgeProperties)
        .def("GetGraphProperties", &GraphInterface::GetGraphProperties)
        .def("InitSignalHandling", &GraphInterface::InitSignalHandling);

    to_python_converter<pair<string,bool>, pair_to_tuple<string,bool> >();

    class_<LibInfo>("mod_info")
        .add_property("name", &LibInfo::GetName)
        .add_property("author", &LibInfo::GetAuthor)
        .add_property("copyright", &LibInfo::GetCopyright)
        .add_property("version", &LibInfo::GetVersion)
        .add_property("license", &LibInfo::GetLicense)
        .add_property("cxxflags", &LibInfo::GetCXXFLAGS)
        .add_property("install_prefix", &LibInfo::GetInstallPrefix)
        .add_property("python_dir", &LibInfo::GetPythonDir);
}
