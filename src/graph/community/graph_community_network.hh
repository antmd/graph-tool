// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007-2012 Tiago de Paula Peixoto <tiago@skewed.de>
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

#ifndef GRAPH_COMMUNITY_NETWORK_HH
#define GRAPH_COMMUNITY_NETWORK_HH

#if (GCC_VERSION >= 40400)
#   include <tr1/unordered_map>
#else
#   include <boost/tr1/unordered_map.hpp>
#endif
#include <iostream>
#include <iomanip>

namespace graph_tool
{

using namespace std;
using namespace boost;

// retrieves the network of communities given a community structure

struct get_community_network
{
    template <class Graph, class CommunityGraph, class CommunityMap,
              class CCommunityMap,
              class VertexWeightMap, class EdgeWeightMap, class EdgeIndex,
              class VertexIndex, class VertexProperty, class EdgeProperty>
    void operator()(const Graph& g, CommunityGraph& cg,
                    VertexIndex cvertex_index, EdgeIndex cedge_index,
                    CommunityMap s_map, CCommunityMap cs_map,
                    VertexWeightMap vweight, EdgeWeightMap eweight,
                    VertexProperty vertex_count, EdgeProperty edge_count,
                    bool self_loops) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;
        typedef typename graph_traits<CommunityGraph>::vertex_descriptor
            cvertex_t;
        typedef typename graph_traits<CommunityGraph>::edge_descriptor
            cedge_t;
        typedef typename boost::property_traits<CommunityMap>::value_type
            s_type;
        typedef typename boost::property_traits<VertexProperty>::value_type
            vprop_type;

        tr1::unordered_map<s_type, pair<vector<vertex_t>, vprop_type>, hash<s_type> >
            comms;
        typename graph_traits<Graph>::vertex_iterator v, v_end;
        for (tie(v, v_end) = vertices(g); v != v_end; ++v)
        {
            pair<vector<vertex_t>, vprop_type>& m = comms[get(s_map, *v)];
            m.first.push_back(*v);
            m.second += get(vweight, *v);
        }

        // create vertices
        tr1::unordered_map<s_type, cvertex_t, hash<s_type> >
            comm_vertices;
        for (typeof(comms.begin()) iter = comms.begin(); iter != comms.end();
             ++iter)
        {
            cvertex_t v = add_vertex(cg);
            put(vertex_count, v, iter->second.second);
            comm_vertices[iter->first] = v;
            put_dispatch(cs_map, v, iter->first,
                         typename boost::is_convertible
                            <typename property_traits<CommunityMap>::category,
                             writable_property_map_tag>::type());
        }

        // create edges
        tr1::unordered_map<pair<size_t, size_t>,
                           cedge_t, hash<pair<size_t, size_t> > >
            comm_edges;
        for (typeof(comms.begin()) iter = comms.begin(); iter != comms.end();
             ++iter)
        {
            cvertex_t cs = comm_vertices[iter->first];
            for (size_t i = 0; i < iter->second.first.size(); ++i)
            {
                vertex_t s = iter->second.first[i];
                typename graph_traits<Graph>::out_edge_iterator e, e_end;
                for (tie(e, e_end) = out_edges(s, g); e != e_end; ++e)
                {
                    vertex_t t = target(*e, g);
                    cvertex_t ct = comm_vertices[get(s_map, t)];
                    if (ct == cs && !self_loops)
                        continue;
                    cedge_t ce;
                    if (comm_edges.find(make_pair(cs, ct)) != comm_edges.end())
                        ce = comm_edges[make_pair(cs, ct)];
                    else if (!is_directed::apply<Graph>::type::value &&
                             comm_edges.find(make_pair(ct, cs)) != comm_edges.end())
                        ce = comm_edges[make_pair(ct, cs)];
                    else
                    {
                        ce = add_edge(cs, ct, cg).first;
                        comm_edges[make_pair(cs, ct)] = ce;
                        put(cedge_index, ce, comm_edges.size() - 1);
                    }
                    put(edge_count, ce, get(edge_count, ce) + get(eweight, *e));
                }
            }
        }
    }

    template <class PropertyMap>
    void put_dispatch(PropertyMap cs_map,
                      const typename property_traits<PropertyMap>::key_type& v,
                      const typename property_traits<PropertyMap>::value_type& val,
                      mpl::true_ is_writable) const
    {
        put(cs_map, v, val);
    }

    template <class PropertyMap>
    void put_dispatch(PropertyMap cs_map,
                      const typename property_traits<PropertyMap>::key_type& v,
                      const typename property_traits<PropertyMap>::value_type& val,
                      mpl::false_ is_writable) const
    {
    }

};

} // graph_tool namespace

#endif // GRAPH_COMMUNITY_NETWORK_HH
