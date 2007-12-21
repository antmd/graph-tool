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
#include <boost/tuple/tuple.hpp>
#include "graph.hh"
#include "graph_python_interface.hh"

using namespace std;
using namespace graph_tool;
using namespace boost;
using namespace boost::python;

// some conversions...

template <class T1, class T2>
struct pair_to_tuple
{
    static PyObject* convert(const pair<T1,T2>& p)
    {
        boost::python::tuple t = boost::python::make_tuple(p.first,p.second);
        return incref(t.ptr());
    }
};

template <class T1, class T2, class T3>
struct tuple_to_tuple
{
    static PyObject* convert(const boost::tuple<T1,T2,T3>& p)
    {
        boost::python::tuple t =
            boost::python::make_tuple(get<0>(p),get<1>(p),get<2>(p));
        return incref(t.ptr());
    }
};


template <class T1, class T2, class T3>
struct tuple_from_tuple
{
    tuple_from_tuple()
    {
        converter::registry::push_back
            (&convertible, &construct,
             boost::python::type_id<boost::tuple<T1,T2,T3> >());
    }

    static void* convertible(PyObject* obj_ptr)
    {
        handle<> x(borrowed(obj_ptr));
        object o(x);
        extract<T1> first(o[0]);
        extract<T2> second(o[1]);
        extract<T2> third(o[3]);
        if (!first.check() || !second.check() || !third.check())
            return 0;
        return obj_ptr;
    }

    static void construct(PyObject* obj_ptr,
                          converter::rvalue_from_python_stage1_data* data)
    {
        handle<> x(borrowed(obj_ptr));
        object o(x);
        boost::tuple<T1,T2,T3> value;
        get<0>(value) = extract<T1>(o[0]);
        get<1>(value) = extract<T2>(o[1]);
        get<2>(value) = extract<T2>(o[2]);
        void* storage =
            ((boost::python::converter::rvalue_from_python_storage
              <boost::tuple<T1,T2,T3> >*) data)->storage.bytes;
        new (storage) boost::tuple<T1,T2,T3>(value);
        data->convertible = storage;
    }
};

template <class T1, class T2>
struct pair_from_tuple
{
    pair_from_tuple()
    {
        converter::registry::push_back(&convertible, &construct,
                                       boost::python::type_id<pair<T1,T2> >());
    }

    static void* convertible(PyObject* obj_ptr)
    {
        handle<> x(borrowed(obj_ptr));
        object o(x);
        extract<T1> first(o[0]);
        extract<T2> second(o[1]);
        if (!first.check() || !second.check())
            return 0;
        return obj_ptr;
    }

    static void construct(PyObject* obj_ptr,
                          converter::rvalue_from_python_stage1_data* data)
    {
        handle<> x(borrowed(obj_ptr));
        object o(x);
        pair<T1,T2> value;
        value.first = extract<T1>(o[0]);
        value.second = extract<T2>(o[1]);
        void* storage =
            ( (boost::python::converter::rvalue_from_python_storage
               <pair<T1,T2> >*) data)->storage.bytes;
        new (storage) pair<T1,T2>(value);
        data->convertible = storage;
    }
};

template <class ValueType>
struct variant_from_python
{
    variant_from_python()
    {
        converter::registry::push_back
            (&convertible, &construct,
             boost::python::type_id<GraphInterface::deg_t>());
    }

    static void* convertible(PyObject* obj_ptr)
    {
        handle<> x(borrowed(obj_ptr));
        object o(x);
        extract<ValueType> str(o);
        if (!str.check())
            return 0;
        return obj_ptr;
    }

    static void construct(PyObject* obj_ptr,
                          converter::rvalue_from_python_stage1_data* data)
    {
        handle<> x(borrowed(obj_ptr));
        object o(x);
        ValueType value = extract<ValueType>(o);
        GraphInterface::deg_t deg = value;
        void* storage =
            ( (boost::python::converter::rvalue_from_python_storage
               <GraphInterface::deg_t>*) data)->storage.bytes;
        new (storage) GraphInterface::deg_t(deg);
        data->convertible = storage;
    }
};



template <class Hist>
struct hist_to_dict
{
    static PyObject* convert(const Hist& h)
    {
        dict hist;
        for(typeof(h.begin()) iter = h.begin(); iter != h.end(); ++iter)
            hist[iter->first] = iter->second;
        return incref(hist.ptr());
    }
};

struct pos_t_to_tuple
{
    static PyObject* convert(const pos_t& p)
    {
        boost::python::tuple t = boost::python::make_tuple(p.x,p.y);
        return incref(t.ptr());
    }
};

struct pos_t_from_tuple
{
    pos_t_from_tuple()
    {
        converter::registry::push_back(&convertible, &construct,
                                       boost::python::type_id<pos_t>());
    }

    static void* convertible(PyObject* obj_ptr)
    {
        handle<> x(borrowed(obj_ptr));
        object o(x);
        extract<double> first(o[0]);
        extract<double> second(o[1]);
        if (!first.check() || !second.check())
            return 0;
        return obj_ptr;
    }

    static void construct(PyObject* obj_ptr,
                          converter::rvalue_from_python_stage1_data* data)
    {
        handle<> x(borrowed(obj_ptr));
        object o(x);
        pos_t value;
        value.x = extract<double>(o[0]);
        value.y = extract<double>(o[1]);
        void* storage =
            ( (boost::python::converter::rvalue_from_python_storage
               <pos_t>*) data)->storage.bytes;
        new (storage) pos_t(value);
        data->convertible = storage;
    }
};

struct LibInfo
{
    string GetName()      const {return PACKAGE_NAME;}
    string GetAuthor()    const {return AUTHOR;}
    string GetCopyright() const {return COPYRIGHT;}
    string GetVersion()   const {return VERSION " (commit " GIT_COMMIT ", " GIT_COMMIT_DATE ")";}
    string GetLicense()   const {return "GPL version 3 or above";}
    string GetCXXFLAGS()  const {return CXXFLAGS;}
    string GetInstallPrefix() const {return INSTALL_PREFIX;}
    string GetPythonDir() const {return PYTHON_DIR;}
};

// overloads
void  (GraphInterface::*ReadFromFile1) (string) =
    &GraphInterface::ReadFromFile;
void  (GraphInterface::*ReadFromFile2) (string, string) =
    &GraphInterface::ReadFromFile;
void  (GraphInterface::*WriteToFile1)  (string) =
    &GraphInterface::WriteToFile;
void  (GraphInterface::*WriteToFile2)  (string, string)  =
    &GraphInterface::WriteToFile;

BOOST_PYTHON_MODULE(libgraph_tool)
{
    GraphInterface().ExportPythonInterface();

    class_<GraphInterface>("GraphInterface")
        .def("GenerateCorrelatedConfigurationalModel",
             &GraphInterface::GenerateCorrelatedConfigurationalModel)
        .def("GetNumberOfVertices", &GraphInterface::GetNumberOfVertices)
        .def("GetNumberOfEdges", &GraphInterface::GetNumberOfEdges)
        .def("GetVertexHistogram", &GraphInterface::GetVertexHistogram)
        .def("GetEdgeHistogram", &GraphInterface::GetEdgeHistogram)
        .def("LabelComponents", &GraphInterface::LabelComponents)
        .def("LabelParallelEdges", &GraphInterface::LabelParallelEdges)
        .def("GetCombinedVertexHistogram",
             &GraphInterface::GetCombinedVertexHistogram)
        .def("GetAverageCombinedVertexCorrelation",
             &GraphInterface::GetAverageCombinedVertexCorrelation)
        .def("GetVertexCorrelationHistogram",
             &GraphInterface::GetVertexCorrelationHistogram)
        .def("GetEdgeVertexCorrelationHistogram",
             &GraphInterface::GetEdgeVertexCorrelationHistogram)
        .def("GetAverageNearestNeighboursCorrelation",
             &GraphInterface::GetAverageNearestNeighboursCorrelation)
        .def("GetAssortativityCoefficient",
             &GraphInterface::GetAssortativityCoefficient)
        .def("GetScalarAssortativityCoefficient",
             &GraphInterface::GetScalarAssortativityCoefficient)
        .def("GetGlobalClustering", &GraphInterface::GetGlobalClustering)
        .def("SetLocalClusteringToProperty",
             &GraphInterface::SetLocalClusteringToProperty)
        .def("SetExtendedClusteringToProperty",
             &GraphInterface::SetExtendedClusteringToProperty)
        .def("GetDistanceHistogram", &GraphInterface::GetDistanceHistogram)
        .def("GetSampledDistanceHistogram",
             &GraphInterface::GetSampledDistanceHistogram)
        .def("GetReciprocity", &GraphInterface::GetReciprocity)
        .def("GetMinimumSpanningTree",
             &GraphInterface::GetMinimumSpanningTree)
        .def("GetLineGraph", &GraphInterface::GetLineGraph)
        .def("GetBetweenness", &GraphInterface::GetBetweenness)
        .def("GetCentralPointDominance",
             &GraphInterface::GetCentralPointDominance)
        .def("GetCommunityStructure",
             &GraphInterface::GetCommunityStructure)
        .def("GetCommunityNetwork", &GraphInterface::GetCommunityNetwork)
        .def("GetModularity", &GraphInterface::GetModularity)
        .def("RandomRewire", &GraphInterface::RandomRewire)
        .def("SetDirected", &GraphInterface::SetDirected)
        .def("GetDirected", &GraphInterface::GetDirected)
        .def("SetReversed", &GraphInterface::SetReversed)
        .def("GetReversed", &GraphInterface::GetReversed)
        .def("SetVertexFilterProperty",
             &GraphInterface::SetVertexFilterProperty)
        .def("SetVertexFilterRange", &GraphInterface::SetVertexFilterRange)
        .def("IsVertexFilterActive", &GraphInterface::IsVertexFilterActive)
        .def("SetEdgeFilterProperty",
             &GraphInterface::SetEdgeFilterProperty)
        .def("SetEdgeFilterRange", &GraphInterface::SetEdgeFilterRange)
        .def("IsEdgeFilterActive", &GraphInterface::IsEdgeFilterActive)
        .def("EditEdgeProperty",  &GraphInterface::EditEdgeProperty)
        .def("EditVertexProperty",  &GraphInterface::EditVertexProperty)
        .def("EditGraphProperty",  &GraphInterface::EditGraphProperty)
        .def("RemoveEdgeProperty",  &GraphInterface::RemoveEdgeProperty)
        .def("RemoveVertexProperty",  &GraphInterface::RemoveVertexProperty)
        .def("RemoveGraphProperty",  &GraphInterface::RemoveGraphProperty)
        .def("PurgeVertices",  &GraphInterface::PurgeVertices)
        .def("PurgeEdges",  &GraphInterface::PurgeEdges)
        .def("InsertEdgeIndexProperty",
             &GraphInterface::InsertEdgeIndexProperty)
        .def("InsertVertexIndexProperty",
             &GraphInterface::InsertVertexIndexProperty)
        .def("ComputeGraphLayoutGursoy",
             &GraphInterface::ComputeGraphLayoutGursoy)
        .def("ComputeGraphLayoutSpringBlock",
             &GraphInterface::ComputeGraphLayoutSpringBlock)
        .def("WriteToFile", WriteToFile1)
        .def("WriteToFile", WriteToFile2)
        .def("ReadFromFile", ReadFromFile1)
        .def("ReadFromFile", ReadFromFile2)
        .def("Vertices", &GraphInterface::Vertices)
        .def("Edges", &GraphInterface::Edges)
        .def("AddVertex", &GraphInterface::AddVertex)
        .def("AddEdge", &GraphInterface::AddEdge)
        .def("RemoveVertex", &GraphInterface::RemoveVertex)
        .def("RemoveEdge", &GraphInterface::RemoveEdge)
        .def("GetVertexProperties", &GraphInterface::GetVertexProperties)
        .def("GetEdgeProperties", &GraphInterface::GetEdgeProperties)
        .def("GetGraphProperties", &GraphInterface::GetGraphProperties)
        .def("InitSignalHandling", &GraphInterface::InitSignalHandling);

    enum_<GraphInterface::degree_t>("Degree")
        .value("In", GraphInterface::IN_DEGREE)
        .value("Out", GraphInterface::OUT_DEGREE)
        .value("Total", GraphInterface::TOTAL_DEGREE);

    enum_<GraphInterface::comm_corr_t>("CommCorr")
        .value("ErdosReyni", GraphInterface::ERDOS_REYNI)
        .value("Uncorrelated", GraphInterface::UNCORRELATED)
        .value("Correlated", GraphInterface::CORRELATED);

    variant_from_python<string>();
    variant_from_python<GraphInterface::degree_t>();
    to_python_converter<pair<double,double>, pair_to_tuple<double,double> >();
    pair_from_tuple<double,double>();
    to_python_converter<boost::tuple<double,double,double>,
                        tuple_to_tuple<double,double,double> >();
    tuple_from_tuple<double,double,double>();
    to_python_converter<pos_t, pos_t_to_tuple>();
    pos_t_from_tuple();
    pair_from_tuple<bool,bool>();
    to_python_converter<hist_t,hist_to_dict<hist_t> >();
    to_python_converter<hist2d_t,hist_to_dict<hist2d_t> >();
    to_python_converter<hist3d_t,hist_to_dict<hist3d_t> >();
    to_python_converter<avg_corr_t,hist_to_dict<avg_corr_t> >();

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
