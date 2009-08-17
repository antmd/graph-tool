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

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_properties.hh"
#include "graph_util.hh"

#include <iostream>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/categories.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/python/extract.hpp>
#include <boost/lexical_cast.hpp>

#include "graph_python_interface.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

//
// Data type string representation
// ===============================
//
// String representation of individual data types. We have to take care
// specifically that no information is lost with floating point I/O.

namespace boost
{

template <>
string lexical_cast<string,uint8_t>(const uint8_t& val)
{

    // "chars" should be printed as numbers, since they can be non-printable
    return lexical_cast<std::string>(int(val));
}

template <>
uint8_t lexical_cast<uint8_t,string>(const string& val)
{

    // "chars" should be printed as numbers, since they can be non-printable
    return uint8_t(lexical_cast<int>(val));
}

// double and long double should be printed in hexadecimal format to preserve
// internal representation
template <>
string lexical_cast<string,double>(const double& val)
{
    char* str = 0;
    int retval = asprintf(&str, "%la", val);
    if (retval == -1)
        throw bad_lexical_cast();
    std::string ret = str;
    free(str);
    return ret;
}

template <>
double lexical_cast<double,string>(const string& val)
{
    double ret;
    int nc = sscanf(val.c_str(), "%la", &ret);
    if (nc != 1)
        throw bad_lexical_cast();
    return ret;
}

template <>
string lexical_cast<string,long double>(const long double& val)
{
    char* str = 0;
    int retval = asprintf(&str, "%La", val);
    if (retval == -1)
        throw bad_lexical_cast();
    std::string ret = str;
    free(str);
    return ret;
}

template <>
long double lexical_cast<long double,string>(const string& val)
{
    long double ret;
    int nc = sscanf(val.c_str(), "%La", &ret);
    if (nc != 1)
        throw bad_lexical_cast();
    return ret;
}
}

// std::vector<> stream i/o
namespace std
{
// string vectors need special attention, since separators must be properly
// escaped.
template <>
ostream& operator<<(ostream& out, const vector<string>& vec)
{
    for (size_t i = 0; i < vec.size(); ++i)
    {
        string s = vec[i];
        // escape separators
        boost::replace_all(s, "\\", "\\\\");
        boost::replace_all(s, ", ", ",\\ ");

        out << s;
        if (i < vec.size() - 1)
            out << ", ";
    }
    return out;
}

template <>
istream& operator>>(istream& in, vector<string>& vec)
{
    using namespace boost;
    using namespace boost::algorithm;
    using namespace boost::xpressive;

    vec.clear();
    string data;
    while (in.good())
    {
        string line;
        getline(in, line);
        data += line;
    }

    sregex re = sregex::compile(", ");
    sregex_token_iterator iter(data.begin(), data.end(), re, -1), end;
    for (; iter != end; ++iter)
    {
        vec.push_back(*iter);
        // un-escape separators
        boost::replace_all(vec.back(), ",\\ ", ", ");
        boost::replace_all(vec.back(), "\\\\", "\\");
    }
    return in;
}
} // std namespace


//
// Persistent IO of python::object types. All the magic is done in python,
// through the object_pickler and object_unplickler below
//

namespace graph_tool
{
python::object object_pickler;
python::object object_unpickler;
}

namespace boost
{
template <>
string lexical_cast<string,python::object>(const python::object & o)
{
    stringstream s;
    object_pickler(OStream(s), o);
    return s.str();
    return "";
}

template <>
python::object lexical_cast<python::object,string>(const string& ps)
{
    stringstream s(ps);
    python::object o;
    o = object_unpickler(IStream(s));
    return o;
}
}

// the following source & sink provide iostream access to python file-like
// objects

class python_file_device
{
public:
    typedef char                 char_type;
    typedef iostreams::seekable_device_tag  category;

    python_file_device(python::object file): _file(file) {}
    std::streamsize read(char* s, std::streamsize n)
    {
        python::object pbuf = _file.attr("read")(n);
        string buf = python::extract<string>(pbuf);
        for (size_t i = 0; i < buf.size(); ++i)
            s[i] = buf[i];
        return buf.size();
    }

    std::streamsize write(const char* s, std::streamsize n)
    {
        string buf(s, s+n);
        python::object pbuf(buf);
        _file.attr("write")(pbuf);
        return n;
    }

    iostreams::stream_offset seek(iostreams::stream_offset off,
                                  std::ios_base::seekdir way)
    {
        _file.attr("seek")(off, int(way));
        return python::extract<iostreams::stream_offset>(_file.attr("tell")());
    }

private:
    python::object _file;
};

// Property Maps
// =============

struct get_python_property
{
    template <class ValueType, class IndexMap>
    void operator()(ValueType, IndexMap, dynamic_property_map* map,
                    python::object& pmap) const
    {
        typedef typename property_map_type::apply<ValueType, IndexMap>::type
            map_t;
        try
        {
            pmap = python::object
                (PythonPropertyMap<map_t>
                 (dynamic_cast
                  <boost::detail::dynamic_property_map_adaptor<map_t>&>(*map)
                  .base()));
        } catch (bad_cast&) {}
    }
};

template <class IndexMap>
python::object find_property_map(dynamic_property_map* map, IndexMap)
{
    python::object pmap;
    mpl::for_each<value_types>(lambda::bind<void>(get_python_property(),
                                                  lambda::_1,
                                                  IndexMap(), lambda::var(map),
                                                  lambda::var(pmap)));
    return pmap;
}

// this functor will check whether a value is of a specific type, create a
// corresponding vector_property_map and add the value to it

template <class IndexMap>
struct check_value_type
{
    typedef typename IndexMap::key_type key_t;
    check_value_type(IndexMap index_map, const key_t& key,
                     const boost::any& value, dynamic_property_map*& map)
        :_index_map(index_map), _key(key), _value(value), _map(map) {}

    template <class ValueType>
    void operator()(ValueType)
    {
        try
        {
            typedef typename property_map_type::apply<ValueType, IndexMap>::type
                map_t;
            map_t vector_map(_index_map);
            vector_map[_key] = any_cast<ValueType>(_value);
            _map = new boost::detail::dynamic_property_map_adaptor<map_t>
                (vector_map);
        }
        catch (bad_any_cast) {}
    }
    IndexMap _index_map;
    const key_t& _key;
    const boost::any& _value;
    dynamic_property_map*& _map;
};

// this functor will check wether a key is a vertex or edge descriptor, and
// generate the corresponding property map, depending on the value type

template <class VertexIndexMap, class EdgeIndexMap>
struct create_dynamic_map
{
    typedef typename VertexIndexMap::key_type vertex_t;
    typedef typename EdgeIndexMap::key_type edge_t;

    create_dynamic_map(VertexIndexMap vertex_map, EdgeIndexMap edge_map)
        :_vertex_map(vertex_map), _edge_map(edge_map) {}
    auto_ptr<dynamic_property_map> operator()(const string& name,
                                              const boost::any& key,
                                              const boost::any& value)
    {
        dynamic_property_map* map;
        try
        {
            mpl::for_each<value_types>
                (check_value_type<VertexIndexMap>(_vertex_map,
                                                  any_cast<vertex_t>(key),
                                                  value, map));
        }
        catch (bad_any_cast)
        {
            try
            {
                mpl::for_each<value_types>
                    (check_value_type<EdgeIndexMap>(_edge_map,
                                                    any_cast<edge_t>(key),
                                                    value, map));
            }
            catch (bad_any_cast)
            {
                ConstantPropertyMap<size_t,graph_property_tag> graph_index(0);
                mpl::for_each<value_types>
                    (check_value_type<ConstantPropertyMap<size_t,
                                                          graph_property_tag> >
                     (graph_index, any_cast<graph_property_tag>(key),
                      value, map));
            }
        }
        return auto_ptr<dynamic_property_map>(map);
    }

    VertexIndexMap _vertex_map;
    EdgeIndexMap _edge_map;
};

// this graph wrapper will update the edge index map when edges are added

template <class Graph, class EdgeIndexMap>
struct GraphEdgeIndexWrap
{
    GraphEdgeIndexWrap(Graph &g, EdgeIndexMap edge_index_map)
        : _g(g), _edge_index_map(edge_index_map), _n_edges(0) {}
    Graph &_g;
    EdgeIndexMap _edge_index_map;
    size_t _n_edges;

    typedef typename Graph::vertex_property_type vertex_property_type;
    typedef typename Graph::edge_property_type edge_property_type;
    typedef typename Graph::graph_tag graph_tag;
    typedef typename Graph::graph_type graph_type;
};

template <class Graph, class EdgeIndexMap>
inline
typename graph_traits
    <GraphEdgeIndexWrap<Graph,EdgeIndexMap> >::vertex_descriptor
add_vertex(GraphEdgeIndexWrap<Graph,EdgeIndexMap>& g)
{
    return add_vertex(g._g);
}

template <class Graph, class EdgeIndexMap>
inline
pair<typename graph_traits
     <GraphEdgeIndexWrap<Graph,EdgeIndexMap> >::edge_descriptor,bool>
add_edge(typename graph_traits
         <GraphEdgeIndexWrap<Graph,EdgeIndexMap> >::vertex_descriptor u,
         typename graph_traits
         <GraphEdgeIndexWrap<Graph,EdgeIndexMap> >::vertex_descriptor v,
         GraphEdgeIndexWrap<Graph,EdgeIndexMap>& g)
{
    Graph& orig = g._g;
    pair<typename graph_traits
         <GraphEdgeIndexWrap<Graph,EdgeIndexMap> >::edge_descriptor,
         bool> retval = add_edge(u,v,orig);
    if (retval.second)
        g._edge_index_map[retval.first] = g._n_edges;
    ++g._n_edges;
    return retval;
}

namespace boost {
template <class Graph, class EdgeIndexMap>
class graph_traits<GraphEdgeIndexWrap<Graph,EdgeIndexMap> >
    : public graph_traits<Graph> {};
}

// this graph wraps an UndirectedAdaptor, but overrides the underlying
// edge_descriptor type with the original type. This will make the edge property
// maps compatible with the original graph, but will break some things which
// are not relevant here

template <class Graph>
struct FakeUndirGraph: public UndirectedAdaptor<Graph>
{
    FakeUndirGraph(const Graph &g): UndirectedAdaptor<Graph>(g) {}
    FakeUndirGraph(UndirectedAdaptor<Graph> &g): UndirectedAdaptor<Graph>(g) {}
};

template <class Graph>
struct FakeEdgeIterator:
    public graph_traits<UndirectedAdaptor<Graph> >::edge_iterator
{
    typedef typename graph_traits<FakeUndirGraph<Graph> >::edge_descriptor
        edge_descriptor;
    typedef typename graph_traits<UndirectedAdaptor<Graph> >::edge_iterator
        edge_iterator;
    FakeEdgeIterator(){}
    FakeEdgeIterator(edge_iterator e): edge_iterator(e) {}
    edge_descriptor operator*() const
    {
        return edge_descriptor(*edge_iterator(*this));
    }

};

namespace boost {
template <class Graph>
struct graph_traits<FakeUndirGraph<Graph> >
    : public graph_traits<UndirectedAdaptor<Graph> >
{
    typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
    typedef FakeEdgeIterator<Graph> edge_iterator;
};
}


//==============================================================================
// ReadFromFile(file, pfile, format)
//==============================================================================

void build_stream
    (boost::iostreams::filtering_stream<boost::iostreams::input>& stream,
     const string& file,  python::object& pfile, std::ifstream& file_stream)
{
    stream.reset();
    if (file == "-")
        stream.push(std::cin);
    else
    {
        if (pfile == python::object())
        {
            file_stream.open(file.c_str(), std::ios_base::in |
                             std::ios_base::binary);
            file_stream.exceptions(ios_base::badbit | ios_base::failbit);
            if (boost::ends_with(file,".gz"))
                stream.push(boost::iostreams::gzip_decompressor());
            if (boost::ends_with(file,".bz2"))
                stream.push(boost::iostreams::bzip2_decompressor());
            stream.push(file_stream);
        }
        else
        {
            python_file_device src(pfile);
            stream.push(src);
        }
    }
    stream.exceptions(ios_base::badbit);
}


python::tuple GraphInterface::ReadFromFile(string file, python::object pfile,
                                           string format)
{
    bool graphviz = false;
    if (format == "dot")
        graphviz = true;
    else if (format != "xml")
        throw ValueException("error reading from file '" + file +
                             "': requested invalid format '" + format + "'");
    try
    {
        boost::iostreams::filtering_stream<boost::iostreams::input>
            stream;
        std::ifstream file_stream;
        build_stream(stream, file, pfile, file_stream);

        create_dynamic_map<vertex_index_map_t,edge_index_map_t>
            map_creator(_vertex_index, _edge_index);
        dynamic_properties dp(map_creator);
        _mg.clear();

        GraphEdgeIndexWrap<multigraph_t,edge_index_map_t> wg(_mg, _edge_index);
        _directed = true;
        try
        {
            if (graphviz)
                read_graphviz(stream, wg, dp, "vertex_name");
            else
                read_graphml(stream, wg, dp);
        }
        catch (const undirected_graph_error&)
        {
            _directed = false;
            file_stream.close();
            if (pfile != python::object())
            {
                python_file_device src(pfile);
                src.seek(0, std::ios_base::beg);
            }
            build_stream(stream, file, pfile, file_stream);
            FakeUndirGraph<GraphEdgeIndexWrap<multigraph_t,edge_index_map_t> >
                ug(wg);
            if (graphviz)
                read_graphviz(stream, ug, dp, "vertex_name");
            else
                read_graphml(stream, ug, dp);
        }
        _nedges = num_edges(_mg);
        _max_edge_index = (_nedges > 0) ? _nedges - 1 : 0;

        python::dict vprops, eprops, gprops;
        for(typeof(dp.begin()) iter = dp.begin(); iter != dp.end(); ++iter)
        {
            if (iter->second->key() == typeid(vertex_t))
                vprops[iter->first] = find_property_map(iter->second,
                                                       _vertex_index);
            else if (iter->second->key() == typeid(edge_t))
                eprops[iter->first] = find_property_map(iter->second,
                                                       _edge_index);
            else
                gprops[iter->first] = find_property_map(iter->second,
                                                        _graph_index);
        }
        return python::make_tuple(vprops, eprops, gprops);
    }
    catch (ios_base::failure &e)
    {
        throw IOException("error reading from file '" + file + "':" + e.what());
    }
    catch (parse_error &e)
    {
        throw IOException("error reading from file '" + file + "':" + e.what());
    }
};

template <class IndexMap>
string graphviz_insert_index(dynamic_properties& dp, IndexMap index_map,
                             bool insert = true)
{
    typedef GraphInterface::vertex_t vertex_t;
    bool found = false;
    for(typeof(dp.begin()) iter = dp.begin(); iter != dp.end();
        ++iter)
        if (iter->first == "vertex_name" &&
            iter->second->key() == typeid(vertex_t))
            found = true;
    if (!found && insert)
        dp.property("vertex_id", index_map);
    if (found)
        return "vertex_name";
    else
        return "vertex_id";
}

// writes a graph to a file

struct write_to_file
{
    template <class Graph, class IndexMap>
    void operator()(ostream& stream, Graph& g, IndexMap index_map,
                    dynamic_properties& dp, bool graphviz) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

        if (graphviz)
        {
            string name = graphviz_insert_index(dp, index_map, false);
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
    void operator()(ostream& stream, Graph& g, IndexMap index_map,
                    dynamic_properties& dp, bool graphviz) const
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

void GraphInterface::WriteToFile(string file, python::object pfile,
                                 string format, python::list props)
{
    bool graphviz = false;
    if (format == "dot")
        graphviz = true;
    else if (format != "xml")
        throw ValueException("error writing to file '" + file +
                             "': requested invalid format '" + format + "'");
    try
    {
        boost::iostreams::filtering_stream<boost::iostreams::output> stream;
        std::ofstream file_stream;
        if (file == "-")
            stream.push(std::cout);
        else
        {
            if (pfile == python::object())
            {
                file_stream.open(file.c_str(), std::ios_base::out |
                                 std::ios_base::binary);
                file_stream.exceptions(ios_base::badbit | ios_base::failbit);
                if (boost::ends_with(file,".gz"))
                    stream.push(boost::iostreams::gzip_compressor());
                if (boost::ends_with(file,".bz2"))
                    stream.push(boost::iostreams::bzip2_compressor());
                stream.push(file_stream);
            }
            else
            {
                python_file_device sink(pfile);
                stream.push(sink);
            }
        }
        stream.exceptions(ios_base::badbit | ios_base::failbit);

        dynamic_properties dp;
        for (int i = 0; i < len(props); ++i)
        {
            dynamic_property_map* pmap =
                any_cast<dynamic_property_map*>
                (python::extract<boost::any>
                 (props[i][1].attr("get_dynamic_map")()));
            dp.insert(python::extract<string>(props[i][0]),
                      auto_ptr<dynamic_property_map>(pmap));
        }

        if (IsVertexFilterActive())
        {
            // vertex indexes must be between the [0, HardNumVertices(g)] range
            typedef tr1::unordered_map<vertex_t, size_t>  map_t;
            map_t vertex_to_index;
            associative_property_map<map_t> index_map(vertex_to_index);
            run_action<>()(*this, lambda::bind<void>(generate_index(),
                                                     lambda::_1,
                                                     index_map))();
            if (graphviz)
                graphviz_insert_index(dp, index_map);

            if (GetDirected())
                run_action<detail::always_directed>()
                    (*this, lambda::bind<void>(write_to_file(),
                                               lambda::var(stream),
                                               lambda::_1,
                                               index_map, lambda::var(dp),
                                               graphviz))();
            else
                run_action<detail::never_directed>()
                    (*this,lambda::bind<void>(write_to_file_fake_undir(),
                                              lambda::var(stream),
                                              lambda::_1, index_map,
                                              lambda::var(dp), graphviz))();
        }
        else
        {
            if (graphviz)
                graphviz_insert_index(dp, _vertex_index);

            if (GetDirected())
                run_action<detail::always_directed>()
                    (*this, lambda::bind<void>(write_to_file(),
                                               lambda::var(stream),
                                               lambda::_1, _vertex_index,
                                               lambda::var(dp),
                                               graphviz))();
            else
                run_action<detail::never_directed>()
                    (*this,lambda::bind<void>(write_to_file_fake_undir(),
                                              lambda::var(stream),
                                              lambda::_1, _vertex_index,
                                              lambda::var(dp),
                                              graphviz))();
        }
        stream.reset();
    }
    catch (ios_base::failure &e)
    {
        throw IOException("error writing to file '" + file + "':" + e.what());
    }
}
