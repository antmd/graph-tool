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

template <class T>
python::object operator |(const python::object& a, const T& b)
{
    return a / b;
}

struct do_astar_search
{
    template <class Graph, class DistanceMap>
    void operator()(const Graph& g, size_t s, DistanceMap dist,
                    boost::any pred_map, boost::any aweight,
                    AStarVisitorWrapper vis, pair<AStarCmp, AStarCmb> cmp,
                    pair<python::object, python::object> range,
                    pair<python::object, python::object> h) const
    {
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;
        typedef typename property_traits<DistanceMap>::value_type dtype_t;
        dtype_t z = python::extract<dtype_t>(range.first);
        dtype_t i = python::extract<dtype_t>(range.second);
        typedef typename property_map_type::
            apply<int32_t, typeof(get(vertex_index, g))>::type pred_t;
        pred_t pred = any_cast<pred_t>(pred_map);
        checked_vector_property_map<default_color_type,
                                    typeof(get(vertex_index, g))>
            color(get(vertex_index, g));
        checked_vector_property_map<dtype_t,
                                    typeof(get(vertex_index, g))>
            cost(get(vertex_index, g));
        DynamicPropertyMapWrap<dtype_t, edge_t> weight(aweight,
                                                       edge_properties());
        astar_search(g, vertex(s, g), AStarH<dtype_t>(h.first, h.second),
                     vis, pred, cost, dist, weight, get(vertex_index, g), color,
                     cmp.first, cmp.second, i, z);
   }
};


void a_star_search(GraphInterface& g, python::object gi, size_t source,
                   boost::any dist_map, boost::any pred_map, boost::any weight,
                   python::object vis, python::object cmp, python::object cmb,
                   python::object zero, python::object inf, python::object h)
{
    run_action<graph_tool::detail::all_graph_views,mpl::true_>()
        (g, std::bind(do_astar_search(),  placeholders::_1, source,
                      placeholders::_2, pred_map, weight,
                      AStarVisitorWrapper(gi, vis), make_pair(AStarCmp(cmp),
                                                              AStarCmb(cmb)),
                      make_pair(zero, inf), make_pair(gi, h)),
         writable_vertex_properties())(dist_map);
}

void export_astar()
{
    using namespace boost::python;
    def("astar_search", &a_star_search);
}
