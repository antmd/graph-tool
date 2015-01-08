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

#include "graph.hh"
#include "graph_selectors.hh"
#include "graph_hits.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

struct get_hits_dispatch
{
    template <class Graph, class VertexIndex, class WeightMap,
              class CentralityMap>
    void operator()(Graph& g, VertexIndex vertex_index, WeightMap w,
                    CentralityMap x, boost::any ay, double epsilon,
                    size_t max_iter, long double& eig) const
    {
        try
        {
            typename CentralityMap::checked_t y = any_cast<typename CentralityMap::checked_t>(ay);
            get_hits()(g, vertex_index, w, x,
                       y.get_unchecked(num_vertices(g)), epsilon, max_iter,
                       eig);
        }
        catch (bad_any_cast&)
        {
            throw GraphException("x and y vertex properties must be of the same type.");
        }
    }
};


long double hits(GraphInterface& g, boost::any w, boost::any x, boost::any y,
                 double epsilon, size_t max_iter)
{
    if (!w.empty() && !belongs<writable_edge_scalar_properties>()(w))
        throw ValueException("edge property must be writable");
    if (!belongs<vertex_floating_properties>()(x))
        throw ValueException("vertex property must be of floating point"
                             " value type");

    typedef ConstantPropertyMap<int, GraphInterface::edge_t> weight_map_t;
    typedef boost::mpl::push_back<writable_edge_scalar_properties, weight_map_t>::type
        weight_props_t;

    if(w.empty())
        w = weight_map_t(1);

    long double eig = 0;
    run_action<>()
        (g, std::bind(get_hits_dispatch(), placeholders::_1, g.GetVertexIndex(),
                      placeholders::_2,  placeholders::_3, y, epsilon, max_iter,
                      std::ref(eig)),
         weight_props_t(),
         vertex_floating_properties())(w, x);
    return eig;
}

void export_hits()
{
    using namespace boost::python;
    def("get_hits", &hits);
}
