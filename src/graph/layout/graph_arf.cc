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

#include "graph_filtering.hh"
#include "graph.hh"
#include "graph_properties.hh"

#include <boost/lambda/bind.hpp>

#include "graph_arf.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

void arf_layout(GraphInterface& g, boost::any pos, boost::any weight, double d,
                double a, double dt, size_t max_iter, double epsilon,
                size_t dim)
{
    typedef ConstantPropertyMap<int32_t,GraphInterface::edge_t> weight_map_t;
    typedef boost::mpl::push_back<edge_scalar_properties, weight_map_t>::type
        edge_props_t;

    if(weight.empty())
        weight = weight_map_t(1);
    run_action<graph_tool::detail::never_directed>()
        (g, std::bind(get_arf_layout(), placeholders::_1, placeholders::_2,
                      placeholders::_3, a, d, dt, epsilon, max_iter, dim),
         vertex_floating_vector_properties(), edge_props_t())(pos, weight);
}

#include <boost/python.hpp>

void export_arf()
{
    boost::python::def("arf_layout", &arf_layout);
}
