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

#ifndef NUMPY_BIND_HH
#define NUMPY_BIND_HH

#include <vector>
#include <boost/python.hpp>

// numpy unique symbol weirdness
#define PY_ARRAY_UNIQUE_SYMBOL graph_tool_numpy
#ifndef NUMPY_EXPORT
#define NO_IMPORT_ARRAY
#endif
#include "arrayobject.h"

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
using namespace boost;
using namespace boost::python;

typedef mpl::map<
    mpl::pair<uint8_t, mpl::int_<NPY_BYTE> >,
    mpl::pair<uint32_t, mpl::int_<NPY_UINT32> >,
    mpl::pair<int32_t, mpl::int_<NPY_INT32> >,
    mpl::pair<int64_t, mpl::int_<NPY_INT64> >,
    mpl::pair<uint64_t, mpl::int_<NPY_UINT64> >,
    mpl::pair<double,  mpl::int_<NPY_DOUBLE> >,
    mpl::pair<long double, mpl::int_<NPY_LONGDOUBLE> >
    > numpy_types;

template <class ValueType>
python::object wrap_vector_owned(vector<ValueType>& vec)
{
    int val_type = mpl::at<numpy_types,ValueType>::type::value;
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
    handle<> x((PyObject*) ndarray);
    object o(x);
    return o;
}

template <class ValueType>
python::object wrap_vector_not_owned(vector<ValueType>& vec)
{
    PyArrayObject* ndarray;
    int val_type = mpl::at<numpy_types,ValueType>::type::value;
    npy_intp size = vec.size();
    if (vec.empty())
        return wrap_vector_owned(vec); // return an _owned_ array of size one.
    else
        ndarray = (PyArrayObject*) PyArray_SimpleNewFromData(1, &size, val_type,
                                                             &vec[0]);
    ndarray->flags = NPY_ALIGNED | NPY_C_CONTIGUOUS | NPY_WRITEABLE;
    handle<> x((PyObject*) ndarray);
    object o(x);
    return o;
}


template <class ValueType, int Dim>
python::object wrap_multi_array_owned(multi_array<ValueType,Dim>& array)
{
    ValueType* new_data = new ValueType[array.num_elements()];
    memcpy(new_data, array.data(), array.num_elements()*sizeof(ValueType));
    int val_type = mpl::at<numpy_types,ValueType>::type::value;
    npy_intp shape[Dim];
    for (int i = 0; i < Dim; ++i)
        shape[i] = array.shape()[i];
    PyArrayObject* ndarray =
        (PyArrayObject*) PyArray_SimpleNewFromData(Dim, shape, val_type,
                                                   new_data);
    ndarray->flags = NPY_ALIGNED | NPY_C_CONTIGUOUS | NPY_OWNDATA |
        NPY_WRITEABLE;
    handle<> x((PyObject*) ndarray);
    object o(x);
    return o;
}

template <class ValueType, int Dim>
python::object wrap_multi_array_not_owned(multi_array<ValueType,Dim>& array)
{
    int val_type = mpl::at<numpy_types,ValueType>::type::value;
    PyArrayObject* ndarray =
        (PyArrayObject*) PyArray_SimpleNewFromData(Dim, array.shape(), val_type,
                                                   array.origin());
    ndarray->flags = NPY_ALIGNED | NPY_C_CONTIGUOUS | NPY_WRITEABLE;
    handle<> x((PyObject*) ndarray);
    object o(x);
    return o;
}

#endif // NUMPY_BIND_HH
