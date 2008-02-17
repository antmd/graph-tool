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

#include "graph_clustering.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

pair<double,double>
GraphInterface::GetGlobalClustering()
{
    double c, c_err;
    bool directed = _directed;
    _directed = false;
    run_action<detail::never_directed>()
        (*this, bind<void>(get_global_clustering(), _1, var(c), var(c_err)))();
    _directed = directed;
    return make_pair(c,c_err);
}

void GraphInterface::SetLocalClusteringToProperty(string property)
{
    boost::any vertex_prop;

    try
    {
        find_property_map(_properties, property, typeid(vertex_t));
        vertex_prop = prop(property, _vertex_index, _properties);
        if (!belongs<vertex_floating_properties>()(vertex_prop))
            throw GraphException("vertex property " + property +
                                 " is not of floating type");
    }
    catch (property_not_found)
    {
        typedef vector_property_map<double,vertex_index_map_t> clust_map_t;
        clust_map_t clust_map(num_vertices(_mg), _vertex_index);
        _properties.property(property, clust_map);
        vertex_prop = clust_map;
    }

    bool directed = _directed;
    _directed = false;
    run_action<detail::never_directed>()
        (*this, bind<void>(set_clustering_to_property(), _1, _2))(vertex_prop);
    _directed = directed;
}
