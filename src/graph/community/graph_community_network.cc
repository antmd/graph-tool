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
#include "graph.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"

#include <boost/bind.hpp>
#include <boost/bind/placeholders.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/python.hpp>

#include "graph_community_network.hh"

using namespace std;
using namespace boost;

using namespace graph_tool;


void community_network(GraphInterface& gi, GraphInterface& cgi,
                       boost::any community_property, boost::any vertex_count,
                       boost::any edge_count, boost::any weight)
{
    typedef DynamicPropertyMapWrap<double,GraphInterface::edge_t> weight_map_t;
    typedef ConstantPropertyMap<double,GraphInterface::edge_t> no_weight_map_t;
    typedef mpl::vector<weight_map_t,no_weight_map_t> weight_properties;

    if (weight.empty())
        weight = no_weight_map_t(1.0);
    else
        weight = weight_map_t(weight, edge_scalar_properties());

    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        vcount_t;
    vcount_t vcount(gi.GetVertexIndex());
    try
    {
        vcount = any_cast<vcount_t>(vertex_count);
    }
    catch (bad_any_cast&)
    {
        throw ValueException("invalid vertex count property");
    }

    typedef property_map_types::apply<mpl::vector<int32_t,double>,
                                      GraphInterface::edge_index_map_t,
                                      mpl::bool_<false> >::type
        ecount_properties;

    if (!belongs<ecount_properties>()(edge_count))
        throw ValueException("invalid edge count property");

     run_action<>()(gi, bind<void>(get_community_network(), _1,
                                   ref(cgi.GetGraph()), cgi.GetVertexIndex(),
                                   cgi.GetEdgeIndex(), _2,
                                   _3, vcount, _4),
                   vertex_properties(), weight_properties(),
                   ecount_properties())
        (community_property, weight, edge_count);
}
