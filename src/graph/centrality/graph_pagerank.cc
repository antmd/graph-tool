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

#include <boost/python.hpp>
#include <boost/bind.hpp>

#include "graph.hh"
#include "graph_selectors.hh"
#include "graph_pagerank.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

size_t pagerank(GraphInterface& g, boost::any rank, boost::any pers,
                boost::any weight, double d, double epsilon, size_t max_iter)
{
    if (!belongs<vertex_floating_properties>()(rank))
        throw ValueException("rank vertex property must have a floating-point value type");

    if (!pers.empty() && !belongs<vertex_scalar_properties>()(pers))
        throw ValueException("personalization vertex property must have a scalar value type");

    typedef ConstantPropertyMap<double, GraphInterface::vertex_t> pers_map_t;
    typedef boost::mpl::push_back<vertex_scalar_properties, pers_map_t>::type
        pers_props_t;

    if(pers.empty())
        pers = pers_map_t(1.0 / g.GetNumberOfVertices());

    typedef ConstantPropertyMap<double, GraphInterface::edge_t> weight_map_t;
    typedef boost::mpl::push_back<edge_scalar_properties, weight_map_t>::type
        weight_props_t;

    if (!weight.empty() && !belongs<edge_scalar_properties>()(weight))
        throw ValueException("weight edge property must have a scalar value type");

    if(weight.empty())
        weight = weight_map_t(1.0);

    size_t iter;
    run_action<>()
        (g, std::bind(get_pagerank(),
                      placeholders::_1, g.GetVertexIndex(), placeholders::_2,
                      placeholders::_3, placeholders::_4, d,
                      epsilon, max_iter, std::ref(iter)),
         vertex_floating_properties(),
         pers_props_t(), weight_props_t())(rank, pers, weight);
    return iter;
}


void export_pagerank()
{
    using namespace boost::python;
    def("get_pagerank", &pagerank);
}
