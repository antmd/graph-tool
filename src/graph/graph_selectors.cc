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

#include <boost/lambda/bind.hpp>
#include <boost/variant/get.hpp>
#include "graph.hh"
#include "graph_selectors.hh"

using namespace graph_tool;
using namespace boost;
using namespace boost::lambda;

// retrieves the appropriate degree selector
boost::any graph_tool::degree_selector(GraphInterface::deg_t deg)
{
    boost::any sel;
    try
    {
        mpl::for_each<selectors>
            (bind<void>(get_degree_selector(), _1,
                        boost::get<GraphInterface::degree_t>(deg),
                        var(sel)));
    }
    catch (bad_get)
    {
        bool found = false;
        mpl::for_each<vertex_properties>
            (bind<void>(get_scalar_selector(), _1, boost::get<boost::any>(deg),
                        var(sel), var(found)));
        if (!found)
            throw GraphException("invalid degree selector");
    }
    return sel;
}
