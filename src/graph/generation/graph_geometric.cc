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

#include "graph.hh"
#include "graph_filtering.hh"

#include "graph_geometric.hh"

#include <boost/python.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;


typedef graph_tool::detail::get_all_graph_views
    ::apply<graph_tool::detail::filt_scalar_type, mpl::bool_<false>,
            mpl::bool_<true>,mpl::bool_<false>,
            mpl::bool_<true>,mpl::bool_<true> >::type graph_views;

typedef property_map_types::apply<mpl::vector<vector<double> >,
                                  GraphInterface::vertex_index_map_t,
                                  mpl::bool_<false> >::type prop_types;


void geometric(GraphInterface& gi, python::object opoints, double r,
               python::object orange, bool periodic, boost::any pos)
{
    python::object shape = opoints.attr("shape");
    size_t size = python::extract<size_t>(shape[0]);
    vector<vector<double> > points(size);
    vector<pair<double, double> > range(python::len(orange));

    size = python::extract<size_t>(shape[1]);
    for(size_t i = 0; i < points.size(); ++i)
    {
        points[i].resize(size);
        for (size_t j = 0; j < points[i].size(); ++j)
            points[i][j] = python::extract<double>(opoints[i][j]);
    }

    for(size_t i = 0; i < range.size(); ++i)
    {
        range[i].first = python::extract<double>(orange[i][0]);
        range[i].second = python::extract<double>(orange[i][1]);
    }

    run_action<graph_views>()(gi, std::bind(get_geometric(), placeholders::_1,
                                            placeholders::_2, std::ref(points),
                                            std::ref(range), r,
                                            periodic),
                              prop_types())(pos);
}
