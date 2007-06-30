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
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/graph/graphviz.hpp>

#include "graph.hh"
#include "histogram.hh"
#include "graph_filtering.hh"
#include "graph_properties.hh"
#include <boost/graph/graphml.hpp>

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

//==============================================================================
// check_value_type()
// this functor will check whether a value is of a specific type, create a 
// corresponding vector_property_map and add the value to it
//==============================================================================
template <class IndexMap>
struct check_value_type
{
    typedef typename IndexMap::key_type key_t;
    check_value_type(IndexMap index_map, const key_t& key, const boost::any& value, dynamic_property_map*& map)
        :_index_map(index_map), _key(key), _value(value), _map(map) {}

    template <class ValueType>
    void operator()(ValueType)
    {
        try
        {
            vector_property_map<ValueType, IndexMap> vector_map(_index_map);
            vector_map[_key] = any_cast<ValueType>(_value);
            _map = new boost::detail::dynamic_property_map_adaptor<vector_property_map<ValueType, IndexMap> >(vector_map);
        }
        catch (bad_any_cast) {}
    }
    IndexMap _index_map;
    const key_t& _key;
    const boost::any& _value;
    dynamic_property_map*& _map;
};

//==============================================================================
// create_dynamic_map()
// this functor will check wether a key is a vertex or edge descriptor, and
// generate the corresponding property map, depending on the value type
//==============================================================================

template <class VertexIndexMap, class EdgeIndexMap>
struct create_dynamic_map
{
    typedef typename VertexIndexMap::key_type vertex_t;
    typedef typename EdgeIndexMap::key_type edge_t;

    create_dynamic_map(VertexIndexMap vertex_map, EdgeIndexMap edge_map) :_vertex_map(vertex_map), _edge_map(edge_map) {}
    auto_ptr<dynamic_property_map> operator()(const string& name, const boost::any& key, const boost::any& value)
    {
        dynamic_property_map* map;
        try
        {
            mpl::for_each<value_types>(check_value_type<VertexIndexMap>(_vertex_map, any_cast<vertex_t>(key), value, map));
        }
        catch (bad_any_cast)
        {
            try 
            {
                mpl::for_each<value_types>(check_value_type<EdgeIndexMap>(_edge_map, any_cast<edge_t>(key), value, map));
            }
            catch (bad_any_cast)
            {
                ConstantPropertyMap<size_t,graph_property_tag> graph_index(0);
                mpl::for_each<value_types>(check_value_type<ConstantPropertyMap<size_t,graph_property_tag> >(graph_index, any_cast<graph_property_tag>(key), value, map));
            }
        }
        return auto_ptr<dynamic_property_map>(map);
    }

    VertexIndexMap _vertex_map;
    EdgeIndexMap _edge_map;
};

//==============================================================================
// GraphEdgeIndexWrap
// this graph wrapper will update the edge index map when edges are added
//==============================================================================

template <class Graph, class EdgeIndexMap>
struct GraphEdgeIndexWrap
{
    GraphEdgeIndexWrap(Graph &g, EdgeIndexMap edge_index_map): _g(g), _edge_index_map(edge_index_map), _n_edges(0) {}
    Graph &_g;
    EdgeIndexMap _edge_index_map;
    size_t _n_edges;

    typedef typename Graph::vertex_property_type vertex_property_type;
    typedef typename Graph::edge_property_type edge_property_type;
    typedef typename Graph::graph_tag graph_tag;
    typedef typename Graph::graph_type graph_type;
};

template <class Graph, class EdgeIndexMap>
inline typename graph_traits<GraphEdgeIndexWrap<Graph,EdgeIndexMap> >::vertex_descriptor
add_vertex(GraphEdgeIndexWrap<Graph,EdgeIndexMap>& g)
{
    return add_vertex(g._g);
}

template <class Graph, class EdgeIndexMap>
inline pair<typename graph_traits<GraphEdgeIndexWrap<Graph,EdgeIndexMap> >::edge_descriptor,bool>
add_edge(typename graph_traits<GraphEdgeIndexWrap<Graph,EdgeIndexMap> >::vertex_descriptor u, 
         typename graph_traits<GraphEdgeIndexWrap<Graph,EdgeIndexMap> >::vertex_descriptor v, 
         GraphEdgeIndexWrap<Graph,EdgeIndexMap>& g)
{
    Graph& orig = g._g; 
    pair<typename graph_traits<GraphEdgeIndexWrap<Graph,EdgeIndexMap> >::edge_descriptor,bool> retval = add_edge(u,v,orig);
    if (retval.second)
        g._edge_index_map[retval.first] = g._n_edges;
    ++g._n_edges;
    return retval;
}

namespace boost {
template <class Graph, class EdgeIndexMap>
class graph_traits<GraphEdgeIndexWrap<Graph,EdgeIndexMap> >: public graph_traits<Graph> {};
}

//==============================================================================
// FakeUndirGraph
// this graph wraps an UndirectedAdaptor, but overrides the underlying
// edge_descriptor type with the original type. This will make the edge property
// maps compatible with the original graph, but will break some things which 
// are not relevant here
//==============================================================================

template <class Graph>
struct FakeUndirGraph: public UndirectedAdaptor<Graph>
{
    FakeUndirGraph(const Graph &g): UndirectedAdaptor<Graph>(g) {}
    FakeUndirGraph(UndirectedAdaptor<Graph> &g): UndirectedAdaptor<Graph>(g) {}
};

template <class Graph>
struct FakeEdgeIterator: public graph_traits<UndirectedAdaptor<Graph> >::edge_iterator
{
    typedef typename graph_traits<FakeUndirGraph<Graph> >::edge_descriptor edge_descriptor;
    typedef typename graph_traits<UndirectedAdaptor<Graph> >::edge_iterator edge_iterator;
    FakeEdgeIterator(){}
    FakeEdgeIterator(edge_iterator e): edge_iterator(e) {}
    edge_descriptor operator*() const
    {
        return edge_descriptor(*edge_iterator(*this));
    }     

};

namespace boost {
template <class Graph>
struct graph_traits<FakeUndirGraph<Graph> >: public graph_traits<UndirectedAdaptor<Graph> > 
{
    typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
    typedef FakeEdgeIterator<Graph> edge_iterator;
};
}


//==============================================================================
// ReadFromFile(file, format)
//==============================================================================


void GraphInterface::ReadFromFile(string file)
{
    bool graphviz = boost::ends_with(file,".dot") || boost::ends_with(file,".dot.gz") || boost::ends_with(file,".dot.bz2");
    if (graphviz)
        ReadFromFile(file, "dot");
    else
        ReadFromFile(file, "xml");
}

void GraphInterface::ReadFromFile(string file, string format)
{
    bool graphviz = false;
    if (format == "dot")
        graphviz = true;
    else if (format != "xml")
        throw GraphException("error reading from file '" + file + "': requested invalid format '" + format + "'");
    try
    {
        boost::iostreams::filtering_stream<boost::iostreams::input> stream;
        std::ifstream file_stream;
        if (file == "-")
            stream.push(std::cin);
        else
        {
            file_stream.open(file.c_str(), std::ios_base::in | std::ios_base::binary);
            file_stream.exceptions(ios_base::badbit | ios_base::failbit);
            if (boost::ends_with(file,".gz"))
                stream.push(boost::iostreams::gzip_decompressor());
            if (boost::ends_with(file,".bz2"))
                stream.push(boost::iostreams::bzip2_decompressor());
            stream.push(file_stream);
        }
        stream.exceptions(ios_base::badbit);

        _properties = dynamic_properties();
        _mg.clear();

        dynamic_properties_copy dp(create_dynamic_map<vertex_index_map_t, edge_index_map_t>(_vertex_index, _edge_index));
        GraphEdgeIndexWrap<multigraph_t,edge_index_map_t> wg(_mg, _edge_index);
        if (_directed)
        {
            if (graphviz)
                read_graphviz(stream, wg, dp, "vertex_name");
            else
                read_graphml(stream, wg, dp);
        }
        else
        {
            FakeUndirGraph<GraphEdgeIndexWrap<multigraph_t,edge_index_map_t> > ug(wg);
            if (graphviz)
                read_graphviz(stream, ug, dp, "vertex_name");
            else
                read_graphml(stream, ug, dp);
        }

        _properties = dp;
    }
    catch (ios_base::failure &e)
    {
        throw GraphException("error reading from file '" + file + "':" + e.what());
    }

};

//==============================================================================
// WriteToFile(file, format)
//==============================================================================
struct write_to_file
{
    template <class Graph, class IndexMap>
    void operator()(ostream& stream, Graph& g, IndexMap index_map, dynamic_properties& dp, bool graphviz) const
    {
        if (graphviz)
        {
            string name;
            try
            {
                find_property_map(dp, "vertex_name", typeid(typename graph_traits<Graph>::vertex_descriptor));
                name = "vertex_name";
            }
            catch (property_not_found)
            {
                name = "vertex_id";
            }
            write_graphviz(stream, g, dp, name);
        }
        else
        {
            write_graphml(stream, g, index_map, dp, true);
        }

    }
};

struct write_to_file_fake_undir: public write_to_file
{
    template <class Graph, class IndexMap>
    void operator()(ostream& stream, Graph& g, IndexMap index_map, dynamic_properties& dp, bool graphviz) const
    {
        typedef typename Graph::original_graph_t graph_t;
        FakeUndirGraph<graph_t> ug(g);
        write_to_file(*this)(stream, ug, index_map, dp, graphviz);
    }
};

struct generate_index
{
    template <class Graph, class IndexMap>
    void operator()(Graph& g, IndexMap index_map) const
    {
        size_t n = 0;
        typename graph_traits<Graph>::vertex_iterator v, v_end;
        for( tie(v, v_end) = vertices(g); v != v_end; ++v)
            index_map[*v] = n++;
    }
};


void GraphInterface::WriteToFile(string file)
{
    bool graphviz = boost::ends_with(file,".dot") || boost::ends_with(file,".dot.gz") || boost::ends_with(file,".dot.bz2");
    if (graphviz)
        WriteToFile(file, "dot");
    else
        WriteToFile(file, "xml");
}

void GraphInterface::WriteToFile(string file, string format)
{
    bool graphviz = false;
    if (format == "dot")
        graphviz = true;
    else if (format != "xml")
        throw GraphException("error writing to file '" + file + "': requested invalid format '" + format + "'");
    try
    {
        boost::iostreams::filtering_stream<boost::iostreams::output> stream;
        std::ofstream file_stream;
        if (file == "-")
            stream.push(std::cout);
        else
        {
            file_stream.open(file.c_str(), std::ios_base::out | std::ios_base::binary);
            file_stream.exceptions(ios_base::badbit | ios_base::failbit);
            if (boost::ends_with(file,".gz"))
                stream.push(boost::iostreams::gzip_compressor());
            if (boost::ends_with(file,".bz2"))
                stream.push(boost::iostreams::bzip2_compressor());
            stream.push(file_stream);
        }
        stream.exceptions(ios_base::badbit | ios_base::failbit);

        dynamic_properties_copy dp = _properties;

        if (IsVertexFilterActive())
        {
            // vertex indexes must be between the [0, HardNumVertices(g)] range
            typedef tr1::unordered_map<graph_traits<multigraph_t>::vertex_descriptor, size_t>  map_t;
            map_t vertex_to_index;
            associative_property_map<map_t> index_map(vertex_to_index);
            check_filter(*this, bind<void>(generate_index(), _1, index_map), reverse_check(), directed_check());
            if (graphviz)
            {
                try
                {
                    find_property_map(dp, "vertex_name", typeid(graph_traits<multigraph_t>::vertex_descriptor));
                }
                catch (property_not_found)
                {
                    dp.property("vertex_id", index_map);
                }
            }
            if (GetDirected())
            {
                check_filter(*this,bind<void>(write_to_file(),var(stream), _1, index_map, var(dp), graphviz),
                             reverse_check(), always_directed());
            }
            else
            {
                check_filter(*this,bind<void>(write_to_file_fake_undir(), var(stream), _1, index_map, var(dp), graphviz),
                             never_reversed(), always_undirected());
            }
        }
        else
        {
            if (graphviz)
            {
                try
                {
                    find_property_map(dp, "vertex_name", typeid(graph_traits<multigraph_t>::vertex_descriptor));
                }
                catch (property_not_found)
                {
                    dp.property("vertex_id", _vertex_index);
                }
            }

            if (GetDirected())
            {
                check_filter(*this,bind<void>(write_to_file(), var(stream), _1, _vertex_index, var(dp), graphviz),
                             reverse_check(), always_directed());
            }
            else
            {
                check_filter(*this,bind<void>(write_to_file_fake_undir(), var(stream), _1, _vertex_index, var(dp), graphviz),
                             never_reversed(), always_undirected());
            }
        }
        stream.reset();
    }
    catch (ios_base::failure &e)
    {
        throw GraphException("error writing to file '" + file + "':" + e.what());
    }

}
