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

#include <boost/python.hpp>
#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_util.hh"
#include "numpy_bind.hh"

#include "graph_selectors.hh"
#include "graph_properties.hh"

#include "graph_incidence.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

void incidence(GraphInterface& g, boost::any vindex, boost::any eindex,
               python::object odata, python::object oi,
               python::object oj)
{
    if (!belongs<vertex_scalar_properties>()(vindex))
        throw ValueException("index vertex property must have a scalar value type");
    if (!belongs<edge_scalar_properties>()(eindex))
        throw ValueException("index edge property must have a scalar value type");

    multi_array_ref<double,1> data = get_array<double,1>(odata);
    multi_array_ref<int32_t,1> i = get_array<int32_t,1>(oi);
    multi_array_ref<int32_t,1> j = get_array<int32_t,1>(oj);
    run_action<>()
        (g, std::bind(get_incidence(),
                      placeholders::_1,  placeholders::_2,  placeholders::_3,
                      std::ref(data), std::ref(i), std::ref(j)),
         vertex_scalar_properties(),
         edge_scalar_properties())(vindex, eindex);

}
