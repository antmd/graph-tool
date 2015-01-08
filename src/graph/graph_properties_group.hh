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

#ifndef GRAPH_PROPERTIES_GROUP_HH
#define GRAPH_PROPERTIES_GROUP_HH

namespace graph_tool
{

template <class Group = boost::mpl::true_, class Edge = boost::mpl::false_>
struct do_group_vector_property
{

    template <class Graph, class VectorPropertyMap, class PropertyMap>
    void operator()(Graph& g, VectorPropertyMap vector_map, PropertyMap map,
                    size_t pos) const
    {
        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i)
        for (i = 0; i < N; ++i)
        {
            typename boost::graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == boost::graph_traits<Graph>::null_vertex())
                continue;
            dispatch_descriptor(g, vector_map, map, v, pos, Edge());
        }
    }

    template <class Graph, class VectorPropertyMap, class PropertyMap,
              class Descriptor>
    inline void dispatch_descriptor(Graph& g, VectorPropertyMap& vector_map,
                                    PropertyMap& map, const Descriptor& v,
                                    size_t pos, boost::mpl::true_) const
    {
        typename boost::graph_traits<Graph>::out_edge_iterator e, e_end;
        for (tie(e,e_end) = out_edges(v, g); e != e_end; ++e)
        {
            typename boost::property_traits<VectorPropertyMap>::value_type& vec =
                vector_map[*e];
            if (vec.size() <= pos)
                vec.resize(pos+1);
            group_or_ungroup(vector_map, map, *e, pos, Group());
        }
    }

    template <class Graph, class VectorPropertyMap, class PropertyMap,
              class Descriptor>
    inline void dispatch_descriptor(Graph&, VectorPropertyMap& vector_map,
                                    PropertyMap& map, const Descriptor& v,
                                    size_t pos, boost::mpl::false_) const
    {
        if (vector_map[v].size() <= pos)
            vector_map[v].resize(pos+1);
        group_or_ungroup(vector_map, map, v, pos, Group());
    }

    template <class VectorPropertyMap, class PropertyMap, class Descriptor>
    inline void group_or_ungroup(VectorPropertyMap& vector_map,
                                 PropertyMap& map, const Descriptor& v,
                                 size_t pos, boost::mpl::true_) const
    {
        convert(get(map,v),  vector_map[v][pos]);
    }

    template <class VectorPropertyMap, class PropertyMap, class Descriptor>
    inline void group_or_ungroup(VectorPropertyMap& vector_map,
                                 PropertyMap& map, const Descriptor& v,
                                 size_t pos, boost::mpl::false_) const
    {
        convert(vector_map[v][pos], map[v]);
    }

    template <class RetVal, class Value>
    inline void convert(const Value& v, RetVal& r)  const
    {
        r = boost::lexical_cast<RetVal>(v);
    }

    template <class RetVal>
    inline void convert(const boost::python::object& v, RetVal& r)  const
    {
        #pragma omp critical
        r = boost::python::extract<RetVal>(v);
    }

    template <class Value>
    inline void convert(const Value& v, boost::python::object& r)  const
    {
        #pragma omp critical
        r = boost::python::object(v);
    }

    template <class Value>
    inline void convert(const Value& v, Value& r)  const
    {
        r = v;
    }

    inline void convert(const boost::python::object& v, boost::python::object& r)  const
    {
       #pragma omp critical
        r = v;
    }

};

} // namespace graph_tool

#endif // GRAPH_PROPERTIES_GROUP_HH
