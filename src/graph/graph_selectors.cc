// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007  Tiago de Paula Peixoto <tiago@forked.de>
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

#include <boost/variant/get.hpp>
#include "graph.hh"
#include "graph_selectors.hh"

using namespace graph_tool;
using namespace boost;

// retrieves the appropriate degree selector
boost::any graph_tool::degree_selector(GraphInterface::deg_t deg)
{
    try
    {
        boost::any degS;
        mpl::for_each<selectors>
            (get_degree_selector(boost::get<GraphInterface::degree_t>(deg),
                                 degS));
        return degS;
    }
    catch (bad_get)
    {
        return boost::get<boost::any>(deg);
    }
}
