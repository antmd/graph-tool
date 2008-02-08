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
#include "histogram.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"


#include <boost/lambda/bind.hpp>
#include <tr1/unordered_set>
#include <iostream>
#include <iomanip>

#include "graph_community_network.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

void GraphInterface::GetCommunityNetwork(string property, string size_property,
                                         string out_file, string format) const
{
    try
    {
        boost::any vertex_prop = prop(property, _vertex_index, _properties);

        dynamic_properties_copy edge_properties;
        for (typeof(_properties.begin()) iter = _properties.begin();
             iter != _properties.end(); ++iter)
            if (iter->second->key() == typeid(edge_t))
                edge_properties.insert(iter->first,
                                       auto_ptr<dynamic_property_map>
                                           (iter->second));

        run_action<>()(*this, 
                       bind<void>(get_community_network(), _1, _2, property,
                                  size_property, var(edge_properties),
                                  out_file, format), vertex_scalar_properties())
            (vertex_prop);
    }
    catch (property_not_found)
    {
        throw GraphException("edge property " + property + " not found");
    }
}
