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


#include "graph_filtering.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"
#include "graph.hh"

#include "graph_assortativity.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;


pair<double,double>
GraphInterface::GetAssortativityCoefficient(GraphInterface::deg_t deg) const
{
    double a, a_err;
    try
    {
        run_action<>()(*this, 
                       bind<void>(get_assortativity_coefficient(), _1, _2,
                                  var(a), var(a_err)), all_selectors())
            (degree_selector(deg, _properties));
    }
    catch (dynamic_get_failure &e)
    {
        throw GraphException("error getting scalar property: " +
                             string(e.what()));
    }

    return make_pair(a, a_err);
}

pair<double,double>
GraphInterface::GetScalarAssortativityCoefficient(GraphInterface::deg_t deg)
    const
{
    double a, a_err;
    try
    {
        run_action<>()(*this, bind<void>(get_scalar_assortativity_coefficient(),
                                         _1, _2, var(a), var(a_err)),
                   all_selectors())
            (degree_selector(deg, _properties));
    }
    catch (dynamic_get_failure &e)
    {
        throw GraphException("error getting scalar property: " +
                             string(e.what()));
    }

    return make_pair(a, a_err);
}
