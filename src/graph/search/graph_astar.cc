// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007-2010 Tiago de Paula Peixoto <tiago@skewed.de>
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
#include "graph_python_interface.hh"

#include <boost/python.hpp>
#include <boost/graph/astar_search.hpp>

#include "graph.hh"
#include "graph_selectors.hh"
#include "graph_util.hh"

#include "graph_astar.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

struct do_astar_search
{
    template <class Graph, class DistanceMap, class WeightMap>
    void operator()(const Graph& g, size_t s, DistanceMap dist,
                    boost::any pred_map, WeightMap weight,
                    AStarVisitorWrapper vis, pair<AStarCmp, AStarCmb> cmp,
                    pair<python::object, python::object> range,
                    pair<python::object, python::object> h) const
    {
        typedef typename property_traits<DistanceMap>::value_type dtype_t;
        dtype_t z = python::extract<dtype_t>(range.first);
        dtype_t i = python::extract<dtype_t>(range.second);
        typedef typename property_map_type::
            apply<int32_t, typeof(get(vertex_index, g))>::type pred_t;
        pred_t pred = any_cast<pred_t>(pred_map);
        astar_search(g, vertex(s, g), AStarH<dtype_t>(h.first, h.second),
                     visitor(vis).weight_map(weight).
                     predecessor_map(pred).
                     distance_map(dist).distance_compare(cmp.first).
                     distance_combine(cmp.second).distance_inf(i).
                     distance_zero(z));
    }
};


void a_star_search(GraphInterface& g, python::object gi, size_t source,
                   boost::any dist_map, boost::any pred_map, boost::any weight,
                   python::object vis, python::object cmp, python::object cmb,
                   python::object zero, python::object inf, python::object h)
{
    run_action<graph_tool::detail::all_graph_views,mpl::true_>()
        (g, bind<void>(do_astar_search(), _1, source, _2, pred_map, _3,
                        AStarVisitorWrapper(gi, vis), make_pair(AStarCmp(cmp),
                                                                AStarCmb(cmb)),
                        make_pair(zero, inf), make_pair(gi, h)),
         writable_vertex_scalar_properties(),
         edge_scalar_properties())
        (dist_map, weight);
}

void export_astar()
{
    using namespace boost::python;
    def("astar_search", &a_star_search);
}
