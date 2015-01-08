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
#include "graph_properties.hh"
#include "graph_filtering.hh"
#include "graph_selectors.hh"

#include "graph_properties_group.hh"

#include <boost/python/extract.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;


typedef boost::mpl::transform<value_types,make_vector>::type vector_types;
typedef property_map_types::apply<vector_types,
                                  GraphInterface::vertex_index_map_t,
                                  boost::mpl::bool_<false> >::type
    vertex_vector_properties;

typedef property_map_types::apply<vector_types,
                                  GraphInterface::edge_index_map_t,
                                  boost::mpl::bool_<false> >::type
    edge_vector_properties;

void group_vector_property(GraphInterface& g, boost::any vector_prop,
                           boost::any prop, size_t pos, bool edge)
{
    if (edge)
        run_action<graph_tool::detail::always_directed_never_reversed>()
            (g, bind(do_group_vector_property<boost::mpl::true_,boost::mpl::true_>(),
                     placeholders::_1, placeholders::_2, placeholders::_3, pos),
             edge_vector_properties(), edge_properties())
            (vector_prop, prop);
    else
        run_action<graph_tool::detail::always_directed_never_reversed>()
            (g, bind(do_group_vector_property<boost::mpl::true_,boost::mpl::false_>(),
                     placeholders::_1, placeholders::_2, placeholders::_3, pos),
             vertex_vector_properties(), vertex_properties())
            (vector_prop, prop);
}
