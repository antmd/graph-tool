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

#include <boost/python/extract.hpp>

#include <iostream>
#include <boost/iostreams/categories.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/xpressive/xpressive.hpp>

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_properties.hh"
#include "graph_util.hh"

#include "graph_python_interface.hh"
#include "str_repr.hh"

#include "gml.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

// use correct smart pointer type for dynamic properties
#if (BOOST_VERSION / 100 % 1000 >= 44)
    #define DP_SMART_PTR boost::shared_ptr
#else
    #define DP_SMART_PTR std::auto_ptr
#endif

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

#include "graph_io_binary.hh"

// the following source & sink provide iostream access to python file-like
// objects

class python_file_device
{
public:
    typedef char                 char_type;
    typedef iostreams::seekable_device_tag  category;

    python_file_device(boost::python::object file): _file(file) {}
    std::streamsize read(char* s, std::streamsize n)
    {
        boost::python::object pbuf = _file.attr("read")(n);
        string buf = boost::python::extract<string>(pbuf);
        for (size_t i = 0; i < buf.size(); ++i)
            s[i] = buf[i];
        return buf.size();
    }

    std::streamsize write(const char* s, std::streamsize n)
    {
        string buf(s, s+n);
#if (PY_MAJOR_VERSION >= 3)
        // in python 3 we need to construct a 'bytes' instance
        PyObject* bytes = PyBytes_FromStringAndSize(s, n);
        boost::python::handle<> x(bytes);
        boost::python::object pbuf(x);
#else
        boost::python::str pbuf(buf);
#endif
        _file.attr("write")(pbuf);
        return n;
    }

    iostreams::stream_offset seek(iostreams::stream_offset off,
                                  std::ios_base::seekdir way)
    {
        _file.attr("seek")(off, int(way));
        return boost::python::extract<iostreams::stream_offset>(_file.attr("tell")());
    }

private:
    boost::python::object _file;
};

// Property Maps
// =============

struct get_python_property
{
    template <class ValueType, class IndexMap>
    void operator()(ValueType, IndexMap, dynamic_property_map& map,
                    boost::python::object& pmap) const
    {
        typedef typename property_map_type::apply<ValueType, IndexMap>::type
            map_t;
        try
        {
            pmap = boost::python::object
                (PythonPropertyMap<map_t>
                 (dynamic_cast
                  <boost::detail::dynamic_property_map_adaptor<map_t>&>(map)
                  .base()));
        } catch (bad_cast&) {}
    }

    template <class ValueType, class IndexMap>
    void operator()(ValueType, IndexMap, boost::any& map,
                    boost::python::object& pmap) const
    {
        typedef typename property_map_type::apply<ValueType, IndexMap>::type
            map_t;
        try
        {
            PythonPropertyMap<map_t> opmap(any_cast<map_t>(map));
            pmap = boost::python::object(opmap);
        }
        catch (bad_any_cast&) {}
    }

};

template <class IndexMap>
boost::python::object find_property_map(dynamic_property_map& map, IndexMap)
{
    boost::python::object pmap;
    boost::mpl::for_each<value_types>(std::bind(get_python_property(),
                                                std::placeholders::_1, IndexMap(), std::ref(map),
                                                std::ref(pmap)));
    return pmap;
}

template <class IndexMap>
boost::python::object find_property_map(boost::any& map, IndexMap)
{
    boost::python::object pmap;
    boost::mpl::for_each<value_types>(std::bind(get_python_property(),
                                                std::placeholders::_1, IndexMap(), std::ref(map),
                                                std::ref(pmap)));
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
    DP_SMART_PTR<dynamic_property_map> operator()(const string&,
                                                  const boost::any& key,
                                                  const boost::any& value)
    {
        dynamic_property_map* map;
        try
        {
            boost::mpl::for_each<value_types>
                (check_value_type<VertexIndexMap>(_vertex_map,
                                                  any_cast<vertex_t>(key),
                                                  value, map));
        }
        catch (bad_any_cast)
        {
            try
            {
                boost::mpl::for_each<value_types>
                    (check_value_type<EdgeIndexMap>(_edge_map,
                                                    any_cast<edge_t>(key),
                                                    value, map));
            }
            catch (bad_any_cast)
            {
                ConstantPropertyMap<size_t,graph_property_tag> graph_index(0);
                boost::mpl::for_each<value_types>
                    (check_value_type<ConstantPropertyMap<size_t,
                                                          graph_property_tag> >
                     (graph_index, any_cast<graph_property_tag>(key),
                      value, map));
            }
        }
        return DP_SMART_PTR<dynamic_property_map>(map);
    }

    VertexIndexMap _vertex_map;
    EdgeIndexMap _edge_map;
};


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

void build_stream(boost::iostreams::filtering_stream<boost::iostreams::input>& stream,
                  const string& file, boost::python::object& pfile,
                  std::ifstream& file_stream)
{
    stream.reset();
    if (file == "-")
        stream.push(std::cin);
    else
    {
        if (pfile == boost::python::object())
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


boost::python::tuple GraphInterface::ReadFromFile(string file,
                                                  boost::python::object pfile,
                                                  string format,
                                                  boost::python::list ignore_vp,
                                                  boost::python::list ignore_ep,
                                                  boost::python::list ignore_gp)
{
    if (format != "gt" && format != "dot" && format != "xml" && format != "gml")
        throw ValueException("error reading from file '" + file +
                             "': requested invalid format '" + format + "'");
    try
    {
        boost::iostreams::filtering_stream<boost::iostreams::input>
            stream;
        std::ifstream file_stream;
        build_stream(stream, file, pfile, file_stream);

        std::unordered_set<std::string> ivp, iep, igp;
        for (int i = 0; i < len(ignore_vp); ++i)
            ivp.insert(boost::python::extract<string>(ignore_vp[i]));
        for (int i = 0; i < len(ignore_ep); ++i)
            iep.insert(boost::python::extract<string>(ignore_ep[i]));
        for (int i = 0; i < len(ignore_gp); ++i)
            igp.insert(boost::python::extract<string>(ignore_gp[i]));

        create_dynamic_map<vertex_index_map_t,edge_index_map_t>
            map_creator(_vertex_index, _edge_index);
        dynamic_properties dp(map_creator);
        *_mg = multigraph_t();

        if (format == "dot")
            _directed = read_graphviz(stream, *_mg, dp, "vertex_name", true,
                                      ivp, iep, igp);
        else if (format == "xml")
            _directed = read_graphml(stream, *_mg, dp, true, true, true,
                                     ivp, iep, igp);
        else if (format == "gml")
            _directed = read_gml(stream, *_mg, dp, ivp, iep, igp);


        boost::python::dict vprops, eprops, gprops;
        if (format == "gt")
        {
            vector<pair<string, boost::any>> agprops, avprops, aeprops;
            stream.exceptions(ios_base::badbit | ios_base::failbit | ios_base::eofbit);
            _directed = read_graph(stream, *_mg, agprops, avprops, aeprops, igp,
                                   ivp, iep);
            for (auto& p : agprops)
                gprops[p.first] = find_property_map(p.second, _graph_index);
            for (auto& p : avprops)
                vprops[p.first] = find_property_map(p.second, _vertex_index);
            for (auto& p : aeprops)
                eprops[p.first] = find_property_map(p.second, _edge_index);
        }
        else
        {
            for(typeof(dp.begin()) iter = dp.begin(); iter != dp.end(); ++iter)
            {
                if (iter->second->key() == typeid(vertex_t))
                    vprops[iter->first] = find_property_map(*iter->second,
                                                            _vertex_index);
                else if (iter->second->key() == typeid(edge_t))
                    eprops[iter->first] = find_property_map(*iter->second,
                                                            _edge_index);
                else
                    gprops[iter->first] = find_property_map(*iter->second,
                                                            _graph_index);
            }
        }
        return boost::python::make_tuple(vprops, eprops, gprops);
    }
    catch (ios_base::failure &e)
    {
        throw IOException("error reading from file '" + file + "':" + e.what());
    }
    catch (parse_error &e)
    {
        throw IOException("error reading from file '" + file + "':" + e.what());
    }
    catch (gml_parse_error &e)
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
                    dynamic_properties& dp, const string& format) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

        if (format == "dot")
        {
            string name = graphviz_insert_index(dp, index_map, false);
            write_graphviz(stream, g, dp, name);
        }
        else if (format == "xml")
        {
            write_graphml(stream, g, index_map, dp, true);
        }
        else if (format == "gml")
        {
            write_gml(stream, g, index_map, dp);
        }
    }
};

struct write_to_file_fake_undir: public write_to_file
{
    template <class Graph, class IndexMap>
    void operator()(ostream& stream, Graph& g, IndexMap index_map,
                    dynamic_properties& dp, const string& format) const
    {
        typedef typename Graph::original_graph_t graph_t;
        FakeUndirGraph<graph_t> ug(g);
        write_to_file(*this)(stream, ug, index_map, dp, format);
    }
};

struct write_to_binary_file
{
    template <class Graph, class IndexMap>
    void operator()(ostream& stream, Graph& g, IndexMap index_map, size_t N,
                    bool directed, vector<pair<string, boost::any >> & gprops,
                    vector<pair<string, boost::any >> & vprops,
                    vector<pair<string, boost::any >> & eprops) const
    {
        write_graph(g, index_map, N, directed, gprops, vprops, eprops, stream);
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

void GraphInterface::WriteToFile(string file, boost::python::object pfile,
                                 string format, boost::python::list props)
{
    if (format != "gt" && format != "xml" && format != "dot" && format != "gml")
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
            if (pfile == boost::python::object())
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

        if (format == "gt")
        {
            typedef property_map_types::apply<value_types,
                                              GraphInterface::graph_index_map_t>::type
                graph_properties;

            vector<pair<string, boost::any>> agprops, avprops, aeprops;
            for (int i = 0; i < python::len(props); ++i)
            {
                string name = python::extract<string>(props[i][0])();
                boost::any p = python::extract<boost::any>(props[i][1].attr("get_map")())();
                if (belongs<graph_properties>()(p))
                    agprops.push_back(make_pair(name, p));
                if (belongs<vertex_properties>()(p))
                    avprops.push_back(make_pair(name, p));
                if (belongs<edge_properties>()(p))
                    aeprops.push_back(make_pair(name, p));
            }

            bool directed = _directed;
            _directed = true;

            if (IsVertexFilterActive())
            {
                // vertex indexes must be between the [0, HardNumVertices(g)] range
                vector_property_map<size_t> index_map(_vertex_index);
                run_action<>()(*this, std::bind(generate_index(),
                                                std::placeholders::_1,
                                                index_map))();
                run_action<>()(*this, std::bind(write_to_binary_file(),
                                                std::ref(stream),
                                                std::placeholders::_1,
                                                index_map,
                                                GetNumberOfVertices(),
                                                directed,
                                                std::ref(agprops),
                                                std::ref(avprops),
                                                std::ref(aeprops)))();
            }
            else
            {
                run_action<>()(*this, std::bind(write_to_binary_file(),
                                                std::ref(stream),
                                                std::placeholders::_1,
                                                _vertex_index,
                                                GetNumberOfVertices(),
                                                directed,
                                                std::ref(agprops),
                                                std::ref(avprops),
                                                std::ref(aeprops)))();
            }

            _directed = directed;
        }
        else
        {
            dynamic_properties dp;
            for (int i = 0; i < len(props); ++i)
            {
                dynamic_property_map* pmap =
                    any_cast<dynamic_property_map*>
                    (boost::python::extract<boost::any>
                     (props[i][1].attr("get_dynamic_map")()));
                dp.insert(boost::python::extract<string>(props[i][0]),
                          DP_SMART_PTR<dynamic_property_map>(pmap));
            }

            if (IsVertexFilterActive())
            {
                // vertex indexes must be between the [0, HardNumVertices(g)] range
                vector_property_map<size_t> index_map(_vertex_index);
                run_action<>()(*this, boost::bind<void>(generate_index(),
                                                        _1, index_map))();
                if (format == "dot")
                    graphviz_insert_index(dp, index_map);

                if (GetDirected())
                    run_action<detail::always_directed>()
                        (*this, boost::bind<void>(write_to_file(),
                                                  boost::ref(stream), _1,
                                                  index_map, boost::ref(dp),
                                                  format))();
                else
                    run_action<detail::never_directed>()
                        (*this,boost::bind<void>(write_to_file_fake_undir(),
                                                 boost::ref(stream), _1, index_map,
                                                 boost::ref(dp), format))();
            }
            else
            {
                if (format == "dot")
                    graphviz_insert_index(dp, _vertex_index);

                if (GetDirected())
                    run_action<detail::always_directed>()
                        (*this, boost::bind<void>(write_to_file(),
                                                  boost::ref(stream), _1,
                                                  _vertex_index,  boost::ref(dp),
                                                  format))();
                else
                    run_action<detail::never_directed>()
                        (*this,boost::bind<void>(write_to_file_fake_undir(),
                                                 boost::ref(stream), _1,
                                                 _vertex_index, boost::ref(dp),
                                                 format))();
            }
        }
        stream.reset();
    }
    catch (ios_base::failure &e)
    {
        throw IOException("error writing to file '" + file + "':" + e.what());
    }
}
