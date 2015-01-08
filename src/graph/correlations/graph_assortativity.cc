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
#include "graph_selectors.hh"
#include "graph_properties.hh"
#include "graph.hh"

#include "graph_assortativity.hh"

using namespace std;
using namespace graph_tool;

pair<double,double>
assortativity_coefficient(GraphInterface& gi,
                          GraphInterface::deg_t deg)
{
    double a, a_err;
    run_action<>()(gi,std::bind(get_assortativity_coefficient(),
                                placeholders::_1, placeholders::_2,
                                std::ref(a), std::ref(a_err)),
                   scalar_selectors())
        (degree_selector(deg));
    return make_pair(a, a_err);
}

pair<double,double>
scalar_assortativity_coefficient(GraphInterface& gi,
                                 GraphInterface::deg_t deg)
{
    double a, a_err;
    run_action<>()(gi, std::bind(get_scalar_assortativity_coefficient(),
                                 placeholders::_1, placeholders::_2,
                                 std::ref(a), std::ref(a_err)),
                   scalar_selectors())
        (degree_selector(deg));
    return make_pair(a, a_err);
}

#include <boost/python.hpp>

using namespace boost::python;

void export_assortativity()
{
    def("assortativity_coefficient", &assortativity_coefficient);
    def("scalar_assortativity_coefficient", &scalar_assortativity_coefficient);
}
