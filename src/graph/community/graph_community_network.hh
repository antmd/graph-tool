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

#ifndef GRAPH_COMMUNITY_NETWORK_HH
#define GRAPH_COMMUNITY_NETWORK_HH

#include <unordered_map>

#ifdef HAVE_SPARSEHASH
#include SPARSEHASH_INCLUDE(dense_hash_map)
#endif

#include <iostream>
#include <iomanip>

namespace graph_tool
{

using namespace std;
using namespace boost;

// retrieves the network of communities given a community structure

struct get_community_network_vertices
{
    template <class Graph, class CommunityGraph, class CommunityMap,
              class CCommunityMap, class VertexWeightMap,
              class VertexProperty>
    void operator()(const Graph& g, CommunityGraph& cg, CommunityMap s_map,
                    CCommunityMap cs_map, VertexWeightMap vweight,
                    VertexProperty vertex_count) const
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

#ifdef HAVE_SPARSEHASH
        google::dense_hash_map<s_type, vertex_t, std::hash<s_type> > comms;
        comms.set_empty_key(numeric_limits<s_type>::max());
#else
        std::unordered_map<s_type, vertex_t, std::hash<s_type> > comms;
#endif
        // create vertices
        typename graph_traits<Graph>::vertex_iterator vi, vi_end;
        for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
        {
            s_type s = get(s_map, *vi);
            typeof(comms.begin()) iter = comms.find(s);
            cvertex_t v;
            if (iter == comms.end())
            {
                comms[s] = v = add_vertex(cg);
                put_dispatch(cs_map, v, s,
                             typename boost::is_convertible
                             <typename property_traits<CommunityMap>::category,
                             writable_property_map_tag>::type());
            }
            else
            {
                v = iter->second;
            }
            put(vertex_count, v, get(vertex_count, v) + get(vweight, *vi));
        }
    }

    template <class PropertyMap>
    void put_dispatch(PropertyMap cs_map,
                      const typename property_traits<PropertyMap>::key_type& v,
                      const typename property_traits<PropertyMap>::value_type& val,
                      mpl::true_ /*is_writable*/) const
    {
        put(cs_map, v, val);
    }

    template <class PropertyMap>
    void put_dispatch(PropertyMap,
                      const typename property_traits<PropertyMap>::key_type&,
                      const typename property_traits<PropertyMap>::value_type&,
                      mpl::false_ /*is_writable*/) const
    {
    }

};

struct get_community_network_edges
{
    template <class Graph, class CommunityGraph, class CommunityMap,
              class CCommunityMap, class EdgeWeightMap,
              class EdgeIndex, class EdgeProperty>
    void operator()(const Graph& g, CommunityGraph& cg, EdgeIndex,
                    CommunityMap s_map, CCommunityMap cs_map,
                    EdgeWeightMap eweight, EdgeProperty edge_count,
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

#ifdef HAVE_SPARSEHASH
        typedef google::dense_hash_map<s_type, vertex_t, std::hash<s_type> > comms_t;
        comms_t comms(num_vertices(cg));
        comms.set_empty_key(numeric_limits<s_type>::max());
        typedef google::dense_hash_map<cvertex_t, cedge_t> ecomms_t;
#else
        typedef std::unordered_map<s_type, vertex_t, std::hash<s_type> > comms_t;
        comms_t comms(num_vertices(cg));
        typedef std::unordered_map<cvertex_t, cedge_t> ecomms_t;
#endif

        auto index_map = get(vertex_index_t(), cg);
        unchecked_vector_property_map<ecomms_t, decltype(index_map)>
            comm_edges(index_map, num_vertices(cg));

        typename graph_traits<CommunityGraph>::vertex_iterator v, v_end;
        for (tie(v, v_end) = vertices(cg); v != v_end; ++v)
        {
            comms[cs_map[*v]] = *v;
#ifdef HAVE_SPARSEHASH
            comm_edges[*v].set_empty_key(numeric_limits<cvertex_t>::max());
#endif
        }


        // create edges
        typename graph_traits<Graph>::edge_iterator e, e_end;
        for (tie(e, e_end) = edges(g); e != e_end; ++e)
        {
            cvertex_t cs = comms[get(s_map, source(*e, g))];
            cvertex_t ct = comms[get(s_map, target(*e, g))];
            if (ct == cs && !self_loops)
                continue;
            cedge_t ce;

            typeof(comm_edges[cs].begin()) iter = comm_edges[cs].find(ct);
            if (iter != comm_edges[cs].end())
            {
                ce = iter->second;
            }
            else
            {
                if (!is_directed::apply<Graph>::type::value)
                {
                    iter = comm_edges[ct].find(cs);
                    if (iter != comm_edges[ct].end())
                    {
                        ce = iter->second;
                    }
                    else
                    {
                        ce = add_edge(cs, ct, cg).first;
                        comm_edges[cs][ct] = ce;
                    }
                }
                else
                {
                    ce = add_edge(cs, ct, cg).first;
                    comm_edges[cs][ct] = ce;
                }
            }
            put(edge_count, ce, get(edge_count, ce) + get(eweight, *e));
        }
    }
};


template <class T1, class T2>
inline vector<T1> operator*(const vector<T1>& v, const T2& c)
{
    vector<T1> r(v);
    for (size_t i = 0; i < r.size(); ++i)
        r[i] = v[i] * c;
    return r;
}

template <class T1, class T2>
inline void operator/=(vector<T1>& v, const T2& c)
{
    for (size_t i = 0; i < v.size(); ++i)
        v[i] /= c;
}

template <class T1, class T2>
inline vector<T1> operator*(const vector<T1>& v1, const vector<T2>& v2)
{
    vector<T1> r(max(v1.size(), v2.size()));
    for (size_t i = 0; i < min(v1.size(), v2.size()); ++i)
        r[i] = v1[i] * v2[i];
    return r;
}

template <class T1, class T2>
inline void operator/=(vector<T1>& v1, const vector<T2>& v2)
{
    v1.resize(max(v1.size(), v2.size()));
    for (size_t i = 0; i < v2.size(); ++i)
        v1[i] /= v2[i];
    return v1;
}

template <class T1, class T2>
inline void operator+=(vector<T1>& v1, const vector<T2>& v2)
{
    v1.resize(max(v1.size(), v2.size()));
    for (size_t i = 0; i < v2.size(); ++i)
        v1[i] += v2[i];
}


// retrieves the average property of a community structure

struct get_weighted_vertex_property
{
    template <class Graph, class VertexWeightMap, class Vprop>
    void operator()(const Graph& g, VertexWeightMap vweight, Vprop vprop,
                    Vprop temp) const
    {
        typename graph_traits<Graph>::vertex_iterator vi, vi_end;
        for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
            temp[*vi] = vprop[*vi] * get(vweight, *vi);
    }
};


struct get_vertex_community_property_sum
{
    template <class Graph, class CommunityGraph, class CommunityMap,
              class CCommunityMap, class Vprop>
    void operator()(const Graph& g, CommunityGraph& cg, CommunityMap s_map,
                    CCommunityMap cs_map, Vprop vprop, Vprop cvprop) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;
        typedef typename graph_traits<CommunityGraph>::vertex_descriptor
            cvertex_t;
        typedef typename boost::property_traits<CommunityMap>::value_type
            s_type;

#ifdef HAVE_SPARSEHASH
        google::dense_hash_map<s_type, vertex_t, std::hash<s_type> > comms(num_vertices(cg));
        comms.set_empty_key(numeric_limits<s_type>::max());
#else
        std::unordered_map<s_type, vertex_t, std::hash<s_type> > comms(num_vertices(cg));
#endif
        typename graph_traits<CommunityGraph>::vertex_iterator v, v_end;
        for (tie(v, v_end) = vertices(cg); v != v_end; ++v)
            comms[cs_map[*v]] = *v;

        typename graph_traits<Graph>::vertex_iterator vi, vi_end;
        for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
        {
            s_type s = get(s_map, *vi);
            cvprop[comms[s]] += vprop[*vi];
        }
    }
};

struct get_vertex_community_property_norm
{
    template <class CommunityGraph, class VertexCountMap, class Vprop>
    void operator()(const CommunityGraph& cg, VertexCountMap vertex_count,
                    Vprop cvprop) const
    {
        typename graph_traits<CommunityGraph>::vertex_iterator v, v_end;
        for (tie(v, v_end) = vertices(cg); v != v_end; ++v)
            cvprop[*v] /= vertex_count[*v];
    }
};


struct get_weighted_edge_property
{
    template <class Graph, class EdgeWeightMap, class Eprop>
    void operator()(const Graph& g, EdgeWeightMap eweight, Eprop eprop,
                    Eprop temp) const
    {
        typename graph_traits<Graph>::edge_iterator ei, ei_end;
        for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
            temp[*ei] = eprop[*ei] * get(eweight, *ei);
    }
};

struct get_edge_community_property_sum
{
    template <class Graph, class CommunityGraph, class CommunityMap,
              class CCommunityMap, class Eprop>
    void operator()(const Graph& g, CommunityGraph& cg, CommunityMap s_map,
                    CCommunityMap cs_map, Eprop eprop, Eprop ceprop,
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

#ifdef HAVE_SPARSEHASH
        google::dense_hash_map<s_type, vertex_t, std::hash<s_type> > comms(num_vertices(cg));
        comms.set_empty_key(numeric_limits<s_type>::max());
#else
        std::unordered_map<s_type, vertex_t, std::hash<s_type> > comms(num_vertices(cg));
#endif
        typename graph_traits<CommunityGraph>::vertex_iterator v, v_end;
        for (tie(v, v_end) = vertices(cg); v != v_end; ++v)
            comms[cs_map[*v]] = *v;

#ifdef HAVE_SPARSEHASH
        google::dense_hash_map<pair<size_t, size_t>, cedge_t,
                               boost::hash<pair<size_t, size_t> > >
             comm_edges(num_vertices(cg));
        comm_edges.set_empty_key(make_pair(numeric_limits<size_t>::max(),
                                           numeric_limits<size_t>::max()));
#else
        std::unordered_map<pair<size_t, size_t>, cedge_t,
                           boost::hash<pair<size_t, size_t> > >
            comm_edges(num_vertices(cg));
#endif

        typename graph_traits<CommunityGraph>::edge_iterator e, e_end;
        for (tie(e, e_end) = edges(cg); e != e_end; ++e)
        {
            cvertex_t cs = comms[get(cs_map, source(*e, cg))];
            cvertex_t ct = comms[get(cs_map, target(*e, cg))];
            comm_edges[make_pair(cs, ct)] = *e;
        }

        typename graph_traits<Graph>::edge_iterator ei, ei_end;
        for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
        {
            cvertex_t cs = comms[get(s_map, source(*ei, g))];
            cvertex_t ct = comms[get(s_map, target(*ei, g))];
            if (cs == ct && !self_loops)
                continue;
            ceprop[comm_edges[make_pair(cs, ct)]] += eprop[*ei];
        }
    }
};

struct get_edge_community_property_norm
{
    template <class CommunityGraph, class EdgeCountMap,
              class Eprop>
    void operator()(const CommunityGraph& cg, EdgeCountMap edge_count,
                    Eprop ceprop) const
    {
        typename graph_traits<CommunityGraph>::edge_iterator e, e_end;
        for (tie(e, e_end) = edges(cg); e != e_end; ++e)
            ceprop[*e] /= edge_count[*e];
    }
};

} // graph_tool namespace

#endif // GRAPH_COMMUNITY_NETWORK_HH
