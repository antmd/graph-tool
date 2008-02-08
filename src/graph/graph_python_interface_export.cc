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

using namespace std;
using namespace boost;
using namespace graph_tool;

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
        void operator()(Graph*) const
        {
            _pclass
                .def("__getitem__",
                     &pmap_t::template GetValue<PythonVertex>)
                .def("__setitem__",
                     &pmap_t::template SetValue<PythonVertex>)
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
        
        typedef mpl::transform<graph_tool::detail::all_graph_views,
                               mpl::quote1<add_pointer> >::type graph_views;
                      
        mpl::for_each<graph_views>(export_access<ValueType>(pclass));
    }

    string _name;
    const GraphInterface& _gi;
};

void export_python_properties(const GraphInterface& gi)
{
    mpl::for_each<value_types>(export_property_map("PropertyMap", gi));
}
