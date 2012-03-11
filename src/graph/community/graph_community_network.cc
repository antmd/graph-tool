// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007-2011 Tiago de Paula Peixoto <tiago@skewed.de>
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

typedef ConstantPropertyMap<int32_t,GraphInterface::edge_t> no_eweight_map_t;
typedef ConstantPropertyMap<int32_t,GraphInterface::vertex_t> no_vweight_map_t;
typedef DynamicPropertyMapWrap<int32_t,GraphInterface::vertex_t> viweight_map_t;
typedef DynamicPropertyMapWrap<double,GraphInterface::vertex_t> vweight_map_t;
typedef DynamicPropertyMapWrap<int32_t,GraphInterface::edge_t> eiweight_map_t;
typedef DynamicPropertyMapWrap<double,GraphInterface::edge_t> eweight_map_t;
typedef DynamicPropertyMapWrap<python::object,GraphInterface::vertex_t> voweight_map_t;
typedef DynamicPropertyMapWrap<python::object,GraphInterface::edge_t> eoweight_map_t;

struct get_community_network_dispatch
{
    template <class Graph, class CommunityGraph, class CommunityMap,
              class VertexWeightMap, class EdgeWeightMap, class EdgeIndex,
              class VertexIndex>
    void operator()(const Graph& g, CommunityGraph& cg,
                    VertexIndex cvertex_index, EdgeIndex cedge_index,
                    CommunityMap s_map, boost::any acs_map,
                    VertexWeightMap vweight, EdgeWeightMap eweight,
                    pair<boost::any,boost::any> count) const
    {
        typedef typename get_prop_type<CommunityMap, VertexIndex>::type
            comm_t;
        comm_t cs_map = boost::any_cast<comm_t>(acs_map);

        typedef typename mpl::if_<is_same<no_vweight_map_t, VertexWeightMap>,
                                  viweight_map_t, VertexWeightMap>::type vweight_t;
        typedef typename mpl::if_<is_same<no_eweight_map_t, EdgeWeightMap>,
                                  eiweight_map_t, EdgeWeightMap>::type eweight_t;

        vweight_t vertex_count = boost::any_cast<vweight_t>(count.first);
        eweight_t edge_count = boost::any_cast<eweight_t>(count.second);

        get_community_network()(g, cg, cvertex_index, cedge_index, s_map,
                                cs_map, vweight, eweight, vertex_count,
                                edge_count);
    }

    struct get_checked_t
    {
        template <class PropertyMap>
        struct apply
        {
            typedef typename PropertyMap::checked_t type;
        };
    };

    struct get_identity
    {
        template <class PropertyMap>
        struct apply
        {
            typedef PropertyMap type;
        };
    };

    template <class PropertyMap, class IndexMap>
    struct get_prop_type
    {
        typedef typename mpl::if_<typename is_same<PropertyMap, IndexMap>::type,
                                  get_identity,
                                  get_checked_t>::type extract;
        typedef typename extract::template apply<PropertyMap>::type type;
    };

};


void community_network(GraphInterface& gi, GraphInterface& cgi,
                       boost::any community_property,
                       boost::any condensed_community_property,
                       boost::any vertex_count,
                       boost::any edge_count, boost::any vweight,
                       boost::any eweight)
{
    typedef typename mpl::vector<vweight_map_t, voweight_map_t, no_vweight_map_t>::type
        vweight_properties;
    typedef typename mpl::vector<eweight_map_t, eoweight_map_t, no_eweight_map_t>::type
        eweight_properties;

    if (eweight.empty())
    {
        eweight = no_eweight_map_t(1);
        edge_count = eiweight_map_t(edge_count, edge_scalar_properties());
    }
    else
    {
        try
        {
            eweight = eweight_map_t(eweight, edge_scalar_properties());
            edge_count = eweight_map_t(edge_count, edge_scalar_properties());
        }
        catch (...)
        {
            eweight = eoweight_map_t(eweight, edge_properties());
            edge_count = eoweight_map_t(edge_count, edge_properties());
        }
    }

    if (vweight.empty())
    {
        vweight = no_vweight_map_t(1);
        vertex_count = viweight_map_t(vertex_count, vertex_scalar_properties());
    }
    else
    {
        try
        {
            vweight = vweight_map_t(vweight, vertex_scalar_properties());
            vertex_count = vweight_map_t(vertex_count, vertex_scalar_properties());
        }
        catch (...)
        {
            vweight = voweight_map_t(vweight, vertex_properties());
            vertex_count = voweight_map_t(vertex_count, vertex_properties());
        }
    }

     run_action<>()(gi, bind<void>(get_community_network_dispatch(), _1,
                                   ref(cgi.GetGraph()), cgi.GetVertexIndex(),
                                   cgi.GetEdgeIndex(), _2,
                                   condensed_community_property,
                                   _3, _4, make_pair(vertex_count, edge_count)),
                    vertex_scalar_properties(), vweight_properties(),
                    eweight_properties())
        (community_property, vweight, eweight);
     cgi.ReIndexEdges();
}
