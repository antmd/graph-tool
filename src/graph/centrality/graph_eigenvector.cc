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
#include "graph_selectors.hh"
#include "graph_eigenvector.hh"

using namespace std;
using namespace graph_tool;

long double eigenvector(GraphInterface& g, boost::any w, boost::any c,
                        double epsilon, size_t max_iter)
{
    if (!w.empty() && !belongs<writable_edge_scalar_properties>()(w))
        throw ValueException("edge property must be writable");
    if (!belongs<vertex_floating_properties>()(c))
        throw ValueException("vertex property must be of floating point"
                             " value type");

    typedef ConstantPropertyMap<int, GraphInterface::edge_t> weight_map_t;
    typedef boost::mpl::push_back<writable_edge_scalar_properties, weight_map_t>::type
        weight_props_t;

    if(w.empty())
        w = weight_map_t(1);

    long double eig = 0;
    run_action<>()
        (g, std::bind(get_eigenvector(), placeholders::_1, g.GetVertexIndex(),
                      placeholders::_2, placeholders::_3, epsilon, max_iter,
                      std::ref(eig)),
         weight_props_t(),
         vertex_floating_properties())(w, c);
    return eig;
}

#include <boost/python.hpp>

void export_eigenvector()
{
    using namespace boost::python;
    def("get_eigenvector", &eigenvector);
}
