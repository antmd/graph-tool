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
#include <boost/variant/get.hpp>

using namespace graph_tool;
using namespace boost;

// retrieves the appropriate degree selector
boost::any graph_tool::degree_selector(GraphInterface::deg_t deg)
{
    boost::any sel;

    GraphInterface::degree_t* d = boost::get<GraphInterface::degree_t>(&deg);

    if (d != 0)
    {
        mpl::for_each<selectors>
            (std::bind(get_degree_selector(), std::placeholders::_1, *d,
                       std::ref(sel)));
    }
    else
    {
        boost::any* d = boost::get<boost::any>(&deg);
        bool found = false;
        mpl::for_each<vertex_properties>
            (std::bind(get_scalar_selector(), std::placeholders::_1, *d,
                       std::ref(sel), std::ref(found)));
        if (!found)
            throw ValueException("invalid degree selector");
    }
    return sel;
}
