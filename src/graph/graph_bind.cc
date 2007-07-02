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
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#include <boost/python.hpp>
#include <boost/tuple/tuple.hpp>
#include "graph.hh"

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
        boost::python::tuple t = boost::python::make_tuple(get<0>(p),get<1>(p),get<2>(p));
        return incref(t.ptr());
    }
};


template <class T1, class T2, class T3>
struct tuple_from_tuple
{
    tuple_from_tuple()
    {
        converter::registry::push_back(&convertible, &construct,  boost::python::type_id<boost::tuple<T1,T2,T3> >());
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

    static void construct(PyObject* obj_ptr, converter::rvalue_from_python_stage1_data* data)
    {          
        handle<> x(borrowed(obj_ptr));
        object o(x);
        boost::tuple<T1,T2,T3> value;
        get<0>(value) = extract<T1>(o[0]);
        get<1>(value) = extract<T2>(o[1]);
        get<2>(value) = extract<T2>(o[2]);
        void* storage = ( (boost::python::converter::rvalue_from_python_storage<boost::tuple<T1,T2,T3> >*) data)->storage.bytes;
        new (storage) boost::tuple<T1,T2,T3>(value);
        data->convertible = storage;
    }
};

template <class T1, class T2>
struct pair_from_tuple
{
    pair_from_tuple()
    {
        converter::registry::push_back(&convertible, &construct,  boost::python::type_id<pair<T1,T2> >());
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

    static void construct(PyObject* obj_ptr, converter::rvalue_from_python_stage1_data* data)
    {          
        handle<> x(borrowed(obj_ptr));
        object o(x);
        pair<T1,T2> value;
        value.first = extract<T1>(o[0]);
        value.second = extract<T2>(o[1]);
        void* storage = ( (boost::python::converter::rvalue_from_python_storage<pair<T1,T2> >*) data)->storage.bytes;
        new (storage) pair<T1,T2>(value);
        data->convertible = storage;
    }
};

template <class ValueType> 
struct variant_from_python
{
    variant_from_python()
    {
        converter::registry::push_back(&convertible, &construct, boost::python::type_id<GraphInterface::deg_t>());
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

    static void construct(PyObject* obj_ptr, converter::rvalue_from_python_stage1_data* data)
    {          
        handle<> x(borrowed(obj_ptr));
        object o(x);
        ValueType value = extract<ValueType>(o);
        GraphInterface::deg_t deg = value;
        void* storage = ( (boost::python::converter::rvalue_from_python_storage<GraphInterface::deg_t>*) data)->storage.bytes;
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


class GraphInterfaceWrap: public GraphInterface
{
public:
    void GenerateCorrelatedConfigurationalModel(size_t N, object pjk, object ceil_pjk, object inv_ceil_pjk, double ceil_pjk_bound,
                                                object corr, object ceil_corr, object inv_ceil_corr, double ceil_corr_bound, bool undirected_corr, 
                                                size_t seed, bool verbose) 
    {
        GraphInterface& base = *this;
        base.GenerateCorrelatedConfigurationalModel(N, pjk_t(python_function(pjk)), pjk_t(python_function(ceil_pjk)), 
                                                    inv_ceil_t(python_function(inv_ceil_pjk)),
                                                    ceil_pjk_bound, corr_t(python_function(corr)), corr_t(python_function(ceil_corr)), 
                                                    inv_corr_t(python_function(inv_ceil_corr)), ceil_corr_bound, undirected_corr, seed, verbose);
    }

    struct python_function
    {
        python_function(object o): _o(o) {}
        double operator()(size_t j, size_t k)
        {
            return extract<double>(_o(j,k));
        }
        double operator()(size_t jl, size_t kl, size_t j, size_t k)
        {
            return extract<double>(_o(jl,kl,j,k));
        }
        pair<size_t,size_t> operator()(double r1, double r2)
        {
            object retval = _o(r1,r2);
            return make_pair(size_t(max(int(extract<int>(retval[0])),0)), size_t(max(int(extract<int>(retval[1])),0)));
        }
        pair<size_t,size_t> operator()(double r1, double r2, size_t j, size_t k)
        {
            object retval = _o(r1,r2,j,k);
            return make_pair(size_t(max(int(extract<int>(retval[0])),0)), size_t(max(int(extract<int>(retval[1])),0)));
        }
        object _o;
    };

};

struct LibInfo
{
    string GetName()      const {return PACKAGE_NAME;}
    string GetAuthor()    const {return AUTHOR;}
    string GetCopyright() const {return COPYRIGHT;}
    string GetVersion()   const {return VERSION " (r" SVN_REVISION ")";}
};

// overloads
void  (GraphInterfaceWrap::*ReadFromFile1) (string)          = &GraphInterfaceWrap::ReadFromFile;
void  (GraphInterfaceWrap::*ReadFromFile2) (string, string)  = &GraphInterfaceWrap::ReadFromFile;
void  (GraphInterfaceWrap::*WriteToFile1)  (string)          = &GraphInterfaceWrap::WriteToFile;
void  (GraphInterfaceWrap::*WriteToFile2)  (string, string)  = &GraphInterfaceWrap::WriteToFile;


BOOST_PYTHON_MODULE(libgraph_tool)
{
    class_<GraphInterfaceWrap>("GraphInterface")
        .def("GenerateCorrelatedConfigurationalModel", &GraphInterfaceWrap::GenerateCorrelatedConfigurationalModel)
        .def("GetNumberOfVertices", &GraphInterfaceWrap::GetNumberOfVertices)
        .def("GetNumberOfEdges", &GraphInterfaceWrap::GetNumberOfEdges)
        .def("GetVertexHistogram", &GraphInterfaceWrap::GetVertexHistogram)
        .def("GetEdgeHistogram", &GraphInterfaceWrap::GetEdgeHistogram)
        .def("LabelComponents", &GraphInterfaceWrap::LabelComponents)
        .def("LabelParallelEdges", &GraphInterfaceWrap::LabelParallelEdges)
        .def("GetCombinedVertexHistogram", &GraphInterfaceWrap::GetCombinedVertexHistogram)
        .def("GetAverageCombinedVertexCorrelation", &GraphInterfaceWrap::GetAverageCombinedVertexCorrelation)
        .def("GetVertexCorrelationHistogram", &GraphInterfaceWrap::GetVertexCorrelationHistogram)
        .def("GetEdgeVertexCorrelationHistogram", &GraphInterfaceWrap::GetEdgeVertexCorrelationHistogram)
        .def("GetAverageNearestNeighboursCorrelation", &GraphInterfaceWrap::GetAverageNearestNeighboursCorrelation)
        .def("GetAssortativityCoefficient", &GraphInterfaceWrap::GetAssortativityCoefficient)
        .def("GetScalarAssortativityCoefficient", &GraphInterfaceWrap::GetScalarAssortativityCoefficient)
        .def("GetGlobalClustering", &GraphInterfaceWrap::GetGlobalClustering)
        .def("SetLocalClusteringToProperty", &GraphInterfaceWrap::SetLocalClusteringToProperty)
        .def("SetExtendedClusteringToProperty", &GraphInterfaceWrap::SetExtendedClusteringToProperty)
        .def("GetDistanceHistogram", &GraphInterfaceWrap::GetDistanceHistogram)
        .def("GetSampledDistanceHistogram", &GraphInterfaceWrap::GetSampledDistanceHistogram)
        .def("GetReciprocity", &GraphInterfaceWrap::GetReciprocity)
        .def("GetMinimumSpanningTree", &GraphInterfaceWrap::GetMinimumSpanningTree)
        .def("GetLineGraph", &GraphInterfaceWrap::GetLineGraph)
        .def("GetBetweenness", &GraphInterfaceWrap::GetBetweenness)
        .def("GetCentralPointDominance", &GraphInterfaceWrap::GetCentralPointDominance)
        .def("GetCommunityStructure", &GraphInterfaceWrap::GetCommunityStructure)
        .def("GetCommunityNetwork", &GraphInterfaceWrap::GetCommunityNetwork)
        .def("GetModularity", &GraphInterfaceWrap::GetModularity)
        .def("SetDirected", &GraphInterfaceWrap::SetDirected)
        .def("GetDirected", &GraphInterfaceWrap::GetDirected)
        .def("SetReversed", &GraphInterfaceWrap::SetReversed)
        .def("GetReversed", &GraphInterfaceWrap::GetReversed)
        .def("SetVertexFilterProperty", &GraphInterfaceWrap::SetVertexFilterProperty)
        .def("GetVertexFilterProperty", &GraphInterfaceWrap::GetVertexFilterProperty)
        .def("SetVertexFilterRange", &GraphInterfaceWrap::SetVertexFilterRange)
        .def("GetVertexFilterRange", &GraphInterfaceWrap::GetVertexFilterRange)
        .def("IsVertexFilterActive", &GraphInterfaceWrap::IsVertexFilterActive)
        .def("SetGenericVertexFilter",  &GraphInterfaceWrap::SetGenericVertexFilter)
        .def("SetEdgeFilterProperty", &GraphInterfaceWrap::SetEdgeFilterProperty)
        .def("GetEdgeFilterProperty", &GraphInterfaceWrap::GetEdgeFilterProperty)
        .def("SetEdgeFilterRange", &GraphInterfaceWrap::SetEdgeFilterRange)
        .def("GetEdgeFilterRange", &GraphInterfaceWrap::GetEdgeFilterRange)
        .def("IsEdgeFilterActive", &GraphInterfaceWrap::IsEdgeFilterActive)
        .def("SetGenericEdgeFilter",  &GraphInterfaceWrap::SetGenericEdgeFilter)
        .def("EditEdgeProperty",  &GraphInterfaceWrap::EditEdgeProperty)
        .def("EditVertexProperty",  &GraphInterfaceWrap::EditVertexProperty)
        .def("EditGraphProperty",  &GraphInterfaceWrap::EditGraphProperty)
        .def("RemoveEdgeProperty",  &GraphInterfaceWrap::RemoveEdgeProperty)
        .def("RemoveVertexProperty",  &GraphInterfaceWrap::RemoveVertexProperty)
        .def("RemoveGraphProperty",  &GraphInterfaceWrap::RemoveGraphProperty)
        .def("ListProperties",  &GraphInterfaceWrap::ListProperties)
        .def("InsertEdgeIndexProperty",  &GraphInterfaceWrap::InsertEdgeIndexProperty)
        .def("InsertVertexIndexProperty",  &GraphInterfaceWrap::InsertVertexIndexProperty)
        .def("ComputeGraphLayoutGursoy", &GraphInterfaceWrap::ComputeGraphLayoutGursoy)
        .def("ComputeGraphLayoutSpringBlock", &GraphInterfaceWrap::ComputeGraphLayoutSpringBlock)
        .def("WriteToFile", WriteToFile1)
        .def("WriteToFile", WriteToFile2)
        .def("ReadFromFile", ReadFromFile1)
        .def("ReadFromFile", ReadFromFile2)
        .def("InitSignalHandling", &GraphInterfaceWrap::InitSignalHandling);
        
    enum_<GraphInterfaceWrap::degree_t>("Degree")
        .value("In", GraphInterfaceWrap::IN_DEGREE)
        .value("Out", GraphInterfaceWrap::OUT_DEGREE)
        .value("Total", GraphInterfaceWrap::TOTAL_DEGREE);

    enum_<GraphInterfaceWrap::comm_corr_t>("CommCorr")
        .value("ErdosReyni", GraphInterfaceWrap::ERDOS_REYNI)
        .value("Uncorrelated", GraphInterfaceWrap::UNCORRELATED)
        .value("Correlated", GraphInterfaceWrap::CORRELATED);

    variant_from_python<string>();
    variant_from_python<GraphInterfaceWrap::degree_t>();
    to_python_converter<pair<double,double>, pair_to_tuple<double,double> >();
    pair_from_tuple<double,double>();
    to_python_converter<boost::tuple<double,double,double>, tuple_to_tuple<double,double,double> >();
    tuple_from_tuple<double,double,double>();
    to_python_converter<GraphInterfaceWrap::hist_t, hist_to_dict<GraphInterfaceWrap::hist_t> >();
    to_python_converter<GraphInterfaceWrap::hist2d_t, hist_to_dict<GraphInterfaceWrap::hist2d_t> >();
    to_python_converter<GraphInterfaceWrap::hist3d_t, hist_to_dict<GraphInterfaceWrap::hist3d_t> >();
    to_python_converter<GraphInterfaceWrap::avg_corr_t, hist_to_dict<GraphInterfaceWrap::avg_corr_t> >();

    class_<LibInfo>("mod_info")
        .add_property("name", &LibInfo::GetName)
        .add_property("author", &LibInfo::GetAuthor)
        .add_property("copyright", &LibInfo::GetCopyright)
        .add_property("version", &LibInfo::GetVersion);
}
