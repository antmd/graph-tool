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

#include <boost/mpl/contains.hpp>
#include <boost/python/extract.hpp>

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_properties.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

struct graph_copy
{
    template <class GraphDst, class GraphSrc, class DstVertexIndexMap,
              class SrcVertexIndexMap,  class DstEdgeIndexMap,
              class SrcEdgeIndexMap>
    void operator()(GraphDst& dst, GraphSrc& src,
                    DstVertexIndexMap dst_vertex_index,
                    SrcVertexIndexMap src_vertex_index,
                    DstEdgeIndexMap dst_edge_index,
                    SrcEdgeIndexMap src_edge_index) const
    {
        vector<size_t> index_map(num_vertices(src));
        typename graph_traits<GraphSrc>::vertex_iterator v, v_end;
        for (tie(v, v_end) = vertices(src); v != v_end; ++v)
        {
            if (src_vertex_index[*v] >= index_map.size())
                index_map.resize(src_vertex_index[*v]+1);
            typename graph_traits<GraphDst>::vertex_descriptor new_v =
                add_vertex(dst);
            index_map[src_vertex_index[*v]] = dst_vertex_index[new_v];
        }

        typename graph_traits<GraphSrc>::edge_iterator e, e_end;
        for (tie(e, e_end) = edges(src); e != e_end; ++e)
        {
            size_t s = index_map[src_vertex_index[source(*e, src)]];
            size_t t = index_map[src_vertex_index[target(*e, src)]];
            typename graph_traits<GraphDst>::edge_descriptor new_e =
                add_edge(vertex(s,dst), vertex(t,dst), dst).first;
            dst_edge_index[new_e] = src_edge_index[*e];
        }
    }
};

// copy constructor
GraphInterface::GraphInterface(const GraphInterface& gi)
    :_mg(),
     _nedges(gi._nedges),
     _reversed(gi._reversed),
     _directed(gi._directed),
     _vertex_index(get(vertex_index,_mg)),
     _edge_index(get(edge_index_t(),_mg)),
     _free_indexes(gi._free_indexes),
     _max_edge_index(gi._max_edge_index),
     _vertex_filter_map(_vertex_index),
     _vertex_filter_invert(false),
     _vertex_filter_active(false),
     _edge_filter_map(_edge_index),
     _edge_filter_invert(false),
     _edge_filter_active(false)
{

    graph_copy()(_mg, gi._mg, _vertex_index,
                 gi._vertex_index, _edge_index,
                 gi._edge_index);
    // filters will be copied in python
}

//
// Property map copying
// ====================

// handle type convertions

// generic types
template <class Type1, class Type2>
struct convert
{

    Type1 operator()(const Type2& v) const
    {
        return do_convert(v, is_convertible<Type2,Type1>());
    }

    Type1 do_convert(const Type2& v, mpl::bool_<true>) const
    {
        return Type1(v);
    }

    Type1 do_convert(const Type2& v, mpl::bool_<false>) const
    {
        return specific_convert<Type1,Type2>()(v);
    }

    template <class T1, class T2>
    struct specific_convert
    {
        T1 operator()(const T2& v) const
        {
            throw bad_lexical_cast(); // default action
        }
    };

    // specific specializations

    // python::object
    template <class T1>
    struct specific_convert<T1,python::object>
    {
        T1 operator()(const python::object& v) const
        {
            python::extract<Type1> x(v);
            if (x.check())
                return x();
            else
                throw bad_lexical_cast();
        }
    };

    // string
    template <class T1>
    struct specific_convert<T1,string>
    {
        T1 operator()(const string& v) const
        {
            //uint8_t is not char, it is bool!
            if (is_same<T1, uint8_t>::value)
                return convert<T1,int>()(lexical_cast<int>(v));
            else
                return lexical_cast<Type1>(v);
        }
    };

    template <class T2>
    struct specific_convert<string,T2>
    {
        string operator()(const T2& v) const
        {
            //uint8_t is not char, it is bool!
            if (is_same<T2, uint8_t>::value)
                return lexical_cast<string>(convert<int,T2>()(v));
            else
                return lexical_cast<string>(v);
        }
    };

    // vectors
    template <class T1, class T2>
    struct specific_convert<vector<T1>, vector<T2> >
    {
        vector<T1> operator()(const vector<T2>& v) const
        {
            vector<T1> v2(v.size());
            convert<T1,T2> c;
            for (size_t i = 0; i < v.size(); ++i)
                v2[i] = c(v[i]);
            return v2;
        }
    };

};

// python::object to string, to solve ambiguity
template<> template<>
struct convert<string,python::object>::specific_convert<string,python::object>
{
    string operator()(const python::object& v) const
    {
        python::extract<string> x(v);
        if (x.check())
                return x();
        else
            throw bad_lexical_cast();
    }
};


template <class IteratorSel>
struct copy_property
{
    template <class Graph, class PropertySrc,
              class PropertyTgt>
    void operator()(const Graph& tgt, const Graph& src, PropertySrc src_map,
                    PropertyTgt dst_map) const
    {
        typedef typename property_traits<PropertySrc>::value_type val_src;
        typedef typename property_traits<PropertyTgt>::value_type val_tgt;

        try
        {
            convert<val_tgt,val_src> c;
            typename IteratorSel::template apply<Graph>::type vs, vs_end;
            typename IteratorSel::template apply<Graph>::type vt, vt_end;
            tie(vt, vt_end) = IteratorSel::range(tgt);
            for (tie(vs, vs_end) = IteratorSel::range(src); vs != vs_end; ++vs)
            {
                if (vt == vt_end)
                    throw ValueException("Error copying properties: "
                                         "graphs not identical");
                dst_map[*vt] = c(src_map[*vs]);
                ++vt;
            }
        }
        catch (bad_lexical_cast&)
        {
            throw ValueException("property values are not convertible");
        }
    }
};

struct edge_selector
{
    template <class Graph>
    struct apply
    {
        typedef typename graph_traits<Graph>::edge_iterator type;
    };

    template <class Graph>
    static pair<typename apply<Graph>::type,
                typename apply<Graph>::type>
    range(Graph& g)
    {
        return edges(g);
    }
};

struct vertex_selector
{
    template <class Graph>
    struct apply
    {
        typedef typename graph_traits<Graph>::vertex_iterator type;
    };

    template <class Graph>
    static pair<typename apply<Graph>::type,
                typename apply<Graph>::type>
    range(Graph& g)
    {
        return vertices(g);
    }
};

typedef mpl::vector<GraphInterface::multigraph_t> unfiltered;

void GraphInterface::CopyVertexProperty(const GraphInterface& src,
                                        boost::any prop_src,
                                        boost::any prop_tgt)
{
    typedef edge_properties writable_edge_properties;

    run_action<unfiltered>()
        (*this, bind<void>(copy_property<vertex_selector>(),
                           _1, ref(src._mg),  _2, _3),
         vertex_properties(), writable_vertex_properties())
        (prop_src, prop_tgt);
}

void GraphInterface::CopyEdgeProperty(const GraphInterface& src,
                                      boost::any prop_src,
                                      boost::any prop_tgt)
{
    typedef edge_properties writable_edge_properties;

    run_action<unfiltered>()
        (*this, bind<void>(copy_property<edge_selector>(),
                           _1, ref(src._mg), _2, _3),
         edge_properties(), writable_edge_properties())
        (prop_src, prop_tgt);
}
