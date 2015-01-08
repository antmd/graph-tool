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


struct do_astar_search
{
    template <class Graph, class DistanceMap>
    void operator()(const Graph& g, size_t s, DistanceMap dist,
                    pair<boost::any, boost::any> pc, boost::any aweight,
                    AStarVisitorWrapper vis, pair<AStarCmp, AStarCmb> cmp,
                    pair<python::object, python::object> range,
                    pair<python::object, python::object> h) const
    {

        typedef typename property_traits<DistanceMap>::value_type dtype_t;
        dtype_t z = python::extract<dtype_t>(range.first);
        dtype_t i = python::extract<dtype_t>(range.second);

        checked_vector_property_map<default_color_type,
                                    typeof(get(vertex_index, g))>
            color(get(vertex_index, g));
        typedef typename property_map_type::
            apply<int32_t, typeof(get(vertex_index, g))>::type pred_t;
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;
        DynamicPropertyMapWrap<dtype_t, edge_t> weight(aweight,
                                                       edge_properties());
        astar_search_no_init(g, vertex(s, g),
                             AStarH<dtype_t>(h.first, h.second), vis,
                             any_cast<pred_t>(pc.first),
                             any_cast<DistanceMap>(pc.second), dist, weight,
                             color, get(vertex_index, g), cmp.first,
                             cmp.second, i, z);
    }
};


void a_star_search_implicit(GraphInterface& g, python::object gi, size_t source,
                            boost::any dist_map, boost::any pred,
                            boost::any cost, boost::any weight,
                            python::object vis, python::object cmp,
                            python::object cmb, python::object zero,
                            python::object inf, python::object h)
{
    run_action<graph_tool::detail::all_graph_views,mpl::true_>()
        (g, std::bind(do_astar_search(), placeholders::_1, source,
                      placeholders::_2, make_pair(pred, cost),
                      weight,
                      AStarVisitorWrapper(gi, vis), make_pair(AStarCmp(cmp),
                                                              AStarCmb(cmb)),
                      make_pair(zero, inf), make_pair(gi, h)),
         writable_vertex_properties())(dist_map);
}

void export_astar_implicit()
{
    using namespace boost::python;
    def("astar_search_implicit", &a_star_search_implicit);
}
