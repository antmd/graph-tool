// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2015 Tiago de Paula Peixoto <tiago@skewed.de>
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

#ifndef NUMPY_BIND_OLD_HH
#define NUMPY_BIND_OLD_HH

#include <vector>
#include <boost/python.hpp>

#include <boost/array.hpp>
#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>

#include <boost/type_traits.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/map.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/for_each.hpp>

using namespace std;

typedef boost::mpl::map<
    boost::mpl::pair<bool, boost::mpl::int_<NPY_BOOL> >,
    boost::mpl::pair<int8_t, boost::mpl::int_<NPY_INT8> >,
    boost::mpl::pair<uint8_t, boost::mpl::int_<NPY_UINT8> >,
    boost::mpl::pair<int16_t, boost::mpl::int_<NPY_INT16> >,
    boost::mpl::pair<uint16_t, boost::mpl::int_<NPY_UINT16> >,
    boost::mpl::pair<int32_t, boost::mpl::int_<NPY_INT32> >,
    boost::mpl::pair<uint32_t, boost::mpl::int_<NPY_UINT32> >,
    boost::mpl::pair<int64_t, boost::mpl::int_<NPY_INT64> >,
    boost::mpl::pair<uint64_t, boost::mpl::int_<NPY_UINT64> >,
    boost::mpl::pair<float, boost::mpl::int_<NPY_FLOAT> >,
    boost::mpl::pair<double, boost::mpl::int_<NPY_DOUBLE> >,
    boost::mpl::pair<long double, boost::mpl::int_<NPY_LONGDOUBLE> >,
    boost::mpl::pair<std::complex<float>, boost::mpl::int_<NPY_CFLOAT> >,
    boost::mpl::pair<std::complex<double>, boost::mpl::int_<NPY_CDOUBLE> >,
    boost::mpl::pair<std::complex<long double>, boost::mpl::int_<NPY_CLONGDOUBLE> > 
    > numpy_types;

template <class ValueType>
boost::python::object wrap_vector_owned(vector<ValueType>& vec)
{
    int val_type = boost::mpl::at<numpy_types,ValueType>::type::value;
    npy_intp size[1];
    size[0] = vec.size();
    PyArrayObject* ndarray;
    if (vec.empty())
    {
        ndarray = (PyArrayObject*) PyArray_SimpleNew(1, size, val_type);
    }
    else
    {
        ValueType* new_data = new ValueType[vec.size()];
        memcpy(new_data, &vec[0], vec.size()*sizeof(ValueType));
        ndarray = (PyArrayObject*) PyArray_SimpleNewFromData(1, size, val_type,
                                                             new_data);
    }
    ndarray->flags = NPY_ALIGNED | NPY_C_CONTIGUOUS | NPY_OWNDATA |
        NPY_WRITEABLE;
    boost::python::handle<> x((PyObject*) ndarray);
    boost::python::object o(x);
    return o;
}

template <class ValueType>
boost::python::object wrap_vector_not_owned(vector<ValueType>& vec)
{
    PyArrayObject* ndarray;
    int val_type = boost::mpl::at<numpy_types,ValueType>::type::value;
    npy_intp size = vec.size();
    if (vec.empty())
        return wrap_vector_owned(vec); // return an _owned_ array of size one.
    else
        ndarray = (PyArrayObject*) PyArray_SimpleNewFromData(1, &size, val_type,
                                                             &vec[0]);
    ndarray->flags = NPY_ALIGNED | NPY_C_CONTIGUOUS | NPY_WRITEABLE;
    boost::python::handle<> x((PyObject*) ndarray);
    boost::python::object o(x);
    return o;
}


template <class ValueType, int Dim>
boost::python::object wrap_multi_array_owned(boost::multi_array<ValueType,Dim>& array)
{
    ValueType* new_data = new ValueType[array.num_elements()];
    memcpy(new_data, array.data(), array.num_elements()*sizeof(ValueType));
    int val_type = boost::mpl::at<numpy_types,ValueType>::type::value;
    npy_intp shape[Dim];
    for (int i = 0; i < Dim; ++i)
        shape[i] = array.shape()[i];
    PyArrayObject* ndarray =
        (PyArrayObject*) PyArray_SimpleNewFromData(Dim, shape, val_type,
                                                   new_data);
    ndarray->flags = NPY_ALIGNED | NPY_C_CONTIGUOUS | NPY_OWNDATA |
        NPY_WRITEABLE;
    boost::python::handle<> x((PyObject*) ndarray);
    boost::python::object o(x);
    return o;
}

template <class ValueType, int Dim>
boost::python::object wrap_multi_array_not_owned(boost::multi_array<ValueType,Dim>& array)
{
    int val_type = boost::mpl::at<numpy_types,ValueType>::type::value;
    PyArrayObject* ndarray =
        (PyArrayObject*) PyArray_SimpleNewFromData(Dim, array.shape(), val_type,
                                                   array.origin());
    ndarray->flags = NPY_ALIGNED | NPY_C_CONTIGUOUS | NPY_WRITEABLE;
    boost::python::handle<> x((PyObject*) ndarray);
    boost::python::object o(x);
    return o;
}

// get multi_array_ref from numpy ndarrays

template <class ValueType, size_t dim>
class numpy_multi_array: public boost::multi_array_ref<ValueType,dim>
{
    typedef boost::multi_array_ref<ValueType,dim> base_t;
public:
    template <class ExtentList, class StrideList>
    explicit numpy_multi_array(typename base_t::element* data,
                               const ExtentList& sizes,
                               const StrideList& strides)
        :base_t(data, sizes)
    {
        for (int i = 0; i < dim; ++i)
            base_t::stride_list_[i] = strides[i];
    }
};

struct invalid_numpy_conversion:
    public std::exception
{
    string _error;
public:
    invalid_numpy_conversion(const string& error) {_error = error;}
    ~invalid_numpy_conversion() throw () {}
    const char * what () const throw () {return _error.c_str();}
};

template <class ValueType, size_t dim>
boost::multi_array_ref<ValueType,dim> get_array(boost::python::object points)
{
    PyArrayObject* pa = (PyArrayObject*) points.ptr();

    if (pa->nd != dim)
        throw invalid_numpy_conversion("invalid array dimension!");

    if (boost::mpl::at<numpy_types,ValueType>::type::value != pa->descr->type_num)
    {
        using boost::python::detail::gcc_demangle;
        boost::python::handle<> x(boost::python::borrowed((PyObject*)  pa->descr->typeobj));
        boost::python::object dtype(x);
        string type_name = boost::python::extract<string>(boost::python::str(dtype));
        string error = "invalid array value type: " + type_name;
        error += " (id: " + boost::lexical_cast<string>(pa->descr->type_num) + ")";
        error += ", wanted: " + string(gcc_demangle(typeid(ValueType).name()));
        error += " (id: " + boost::lexical_cast<string>(boost::mpl::at<numpy_types,ValueType>::type::value) + ")";
        throw invalid_numpy_conversion(error);
    }

    vector<size_t> shape(pa->nd);
    for (int i = 0; i < pa->nd; ++i)
        shape[i] = pa->dimensions[i];

    vector<size_t> stride(dim);
    for (size_t i = 0; i < dim; ++i)
        stride[i] = pa->strides[i] / sizeof(ValueType);

    return numpy_multi_array<ValueType,dim>((ValueType *) pa->data,
                                            shape, stride);
}

#endif // NUMPY_BIND_OLD_HH
