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

#ifndef GRAPH_IO_BINARY_HH
#define GRAPH_IO_BINARY_HH

#include <iostream>
#include "graph.hh"
#include "graph_properties.hh"
#include "graph_selectors.hh"
#include <unordered_set>

namespace graph_tool
{

const char* _magic = u8"⛾ gt";
size_t _magic_length = 6;
const uint8_t _version = 1;

// deal with endianness

inline bool is_bigendian()
{
    // from: http://esr.ibiblio.org/?p=5095
    return (*(uint16_t *)"\0\xff" < 0x100);
};

template <bool BE, typename T>
void byte_swap(T& p)
{
    if (BE == is_bigendian())
        return;
    char& r = reinterpret_cast<char&>(p);
    std::reverse(&r, &r + sizeof(T));
}


template <typename T>
void write(std::ostream& s, T v)
{
    s.write(reinterpret_cast<const char*>(&v), sizeof(T));
};


template <typename T>
void write(std::ostream& s, const std::vector<T>& v)
{
    uint64_t size = v.size();
    write(s, size);
    s.write(reinterpret_cast<const char*>(v.data()), sizeof(T) * v.size());
};


void write(std::ostream& s, const std::string& v)
{
    uint64_t size = v.size();
    write(s, size);
    s.write(reinterpret_cast<const char*>(v.data()), v.size());
};

void write(std::ostream& s, const std::vector<std::string>& v)
{
    uint64_t size = v.size();
    write(s, size);
    for (auto& x : v)
        write(s, x);
};

void write(std::ostream& s, const boost::python::object& v)
{
    std::string buf = boost::lexical_cast<std::string>(v);
    write(s, buf);
};


template <bool BE, typename T>
void read(std::istream& s, T& v)
{
    s.read(reinterpret_cast<char*>(&v), sizeof(T));
    byte_swap<BE>(v);
};

template <bool BE, typename T>
void skip(std::istream& s, const T&)
{
    s.ignore(sizeof(T));
};

template <bool BE, typename T>
void read(std::istream& s, std::vector<T>& v)
{
    uint64_t size = 0;
    read<BE>(s, size);
    v.resize(size);

    s.read(reinterpret_cast<char*>(v.data()), sizeof(T) * v.size());

    for (auto& x : v)
        byte_swap<BE>(x);
};

template <bool BE, typename T>
void skip(std::istream& s, const std::vector<T>&)
{
    uint64_t size = 0;
    read<BE>(s, size);
    s.ignore(sizeof(T) * size);
};

template <bool BE>
void read(std::istream& s, std::string& v)
{
    uint64_t size = 0;
    read<BE>(s, size);
    v.resize(size);
    s.read(&v[0], v.size());
};


template <bool BE>
void read(std::istream& s, std::vector<std::string>& v)
{
    uint64_t size = 0;
    read<BE>(s, size);
    v.resize(size);

    for (auto& x : v)
        read<BE>(s, x);
};

template <bool BE>
void skip(std::istream& s, const std::string&)
{
    uint64_t size = 0;
    read<BE>(s, size);
    s.ignore(size);
};

template <bool BE>
void read(std::istream& s, boost::python::object& v)
{
    std::string buf;
    read<BE>(s, buf);
    v = boost::lexical_cast<boost::python::object>(buf);
};


template <bool BE>
void skip(std::istream& s, const boost::python::object&)
{
    skip<BE>(s, std::string());
};

template <class Vint, class Graph, class VProp>
void write_adjacency_dispatch(Graph& g, const VProp& vindex, std::ostream& s)
{
    for (auto v : vertices_range(g))
    {
        size_t k = out_degree(v, g);
        std::vector<Vint> us;
        us.reserve(k);
        for (auto e : out_edges_range(v, g))
            us.push_back(vindex[target(e, g)]);
        write(s, us);
    }
}


template <class Graph, class VProp>
void write_adjacency(Graph& g, const VProp& vindex, uint64_t N,
                     bool is_directed, std::ostream& s)
{
    uint8_t directed = is_directed;
    write(s, directed);
    write(s, N);

    if (N <= numeric_limits<uint8_t>::max())
        write_adjacency_dispatch<uint8_t>(g, vindex, s);
    else if (N <= numeric_limits<uint16_t>::max())
        write_adjacency_dispatch<uint16_t>(g, vindex, s);
    else if (N <= numeric_limits<uint32_t>::max())
        write_adjacency_dispatch<uint32_t>(g, vindex, s);
    else
        write_adjacency_dispatch<uint64_t>(g, vindex, s);
}

template <bool BE, class Vint, class Graph>
void read_adjacency_dispatch(Graph& g, size_t N, std::istream& s)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

    for (vertex_t v = 0; v < N; ++v)
    {
        std::vector<Vint> us;
        read<BE>(s, us);
        for (vertex_t u : us)
        {
            if (u >= N)
                throw IOException("error reading graph: vertex index not in range");
            add_edge(v, u, g);
        }
    }
}


template <bool BE, class Graph>
bool read_adjacency(Graph& g, std::istream& s)
{
    uint8_t directed = false;
    read<BE>(s, directed);

    uint64_t N = 0;
    read<BE>(s, N);

    for (size_t i = 0; i < N; ++i)
        add_vertex(g);

    if (N <= numeric_limits<uint8_t>::max())
        read_adjacency_dispatch<BE, uint8_t>(g, N, s);
    else if (N <= numeric_limits<uint16_t>::max())
        read_adjacency_dispatch<BE, uint16_t>(g, N, s);
    else if (N <= numeric_limits<uint32_t>::max())
        read_adjacency_dispatch<BE, uint32_t>(g, N, s);
    else
        read_adjacency_dispatch<BE, uint64_t>(g, N, s);

    return directed;
}

// Property maps

enum class property_type : uint8_t
{
    Graph,
    Vertex,
    Edge
};

using namespace boost;

typedef mpl::joint_view<value_types,
                        mpl::vector<size_t>>::type val_types;

struct graph_range_traits
{
    typedef GraphInterface::graph_index_map_t index_map_t;
    template <class Graph>
    static index_map_t get_index_map(Graph& g) { return index_map_t(0); }

    struct graph_range
    {
        graph_property_tag _g;
        graph_property_tag* begin() { return &_g;     }
        graph_property_tag* end()   { return &_g + 1; }
    };

    template <class Graph>
    static graph_range get_range(Graph&) { return graph_range(); }

    static property_type get_property_id() { return property_type::Graph; }
};

struct vertex_range_traits
{
    typedef GraphInterface::vertex_index_map_t index_map_t;
    template <class Graph>
    static index_map_t get_index_map(Graph& g) { return get(vertex_index, g); }


    template <class Graph>
    IterRange<typename boost::graph_traits<Graph>::vertex_iterator>
    static get_range(Graph& g) { return vertices_range(g); }

    static property_type get_property_id() { return property_type::Vertex; }
};

struct edge_range_traits
{
    typedef GraphInterface::edge_index_map_t index_map_t;
    template <class Graph>
    static index_map_t get_index_map(Graph& g) { return get(edge_index, g); }


    template <class Graph>
    IterRange<typename boost::graph_traits<Graph>::edge_iterator>
    static get_range(Graph& g) { return edges_range(g); }

    static property_type get_property_id() { return property_type::Edge; }
};


template <class RangeTraits>
struct write_property_dispatch
{
    template <class T, class Graph>
    void operator()(T, Graph& g, boost::any& aprop, bool& found,
                    std::ostream& s) const
    {
        try
        {
            typedef typename property_map_type::apply<T, typename RangeTraits::index_map_t>::type pmap_t;
            pmap_t prop = any_cast<pmap_t>(aprop);
            typedef typename mpl::find<val_types, T>::type pos;
            uint8_t val = mpl::distance<typename mpl::begin<val_types>::type, pos>::type::value;
            write(s, val);
            for (auto x : RangeTraits::get_range(g))
                write(s, prop[x]);
            found = true;
        }
        catch (const boost::bad_any_cast&) {}
    }


    template <class Graph>
    void operator()(size_t, Graph& g, boost::any& aprop, bool& found,
                    std::ostream& s) const
    {
        try
        {
            typedef GraphInterface::vertex_index_map_t pmap_t;
            pmap_t prop = any_cast<pmap_t>(aprop);
            typedef typename mpl::find<val_types, int64_t>::type pos;
            uint8_t val = mpl::distance<typename mpl::begin<val_types>::type, pos>::type::value;
            write(s, val);
            int64_t y;
            for (auto x : vertices_range(g))
            {
                y = prop[x];
                write(s, y);
            }
            found = true;
        }
        catch (const boost::bad_any_cast&) {}

        try
        {
            typedef GraphInterface::edge_index_map_t pmap_t;
            pmap_t prop = any_cast<pmap_t>(aprop);
            typedef typename mpl::find<val_types, int64_t>::type pos;
            uint8_t val = mpl::distance<typename mpl::begin<val_types>::type, pos>::type::value;
            write(s, val);
            int64_t y;
            for (auto x : edges_range(g))
            {
                y = prop[x];
                write(s, y);
            }
            found = true;
        }
        catch (const boost::bad_any_cast&) {}
    }

};


template <class RangeTraits, class Graph>
void write_property(Graph& g, std::string& name, boost::any& prop, std::ostream& s)
{
    property_type pt = RangeTraits::get_property_id();
    write(s, pt);
    write(s, name);
    bool found = false;
    mpl::for_each<val_types>(std::bind(write_property_dispatch<RangeTraits>(),
                                       std::placeholders::_1, std::ref(g),
                                       std::ref(prop), std::ref(found),
                                       std::ref(s)));
    if (!found)
        throw GraphException("Error writing graph: unknown property map type (this is a bug)");
}


template <bool BE, class RangeTraits>
struct read_property_dispatch
{
    template <class T, class Graph>
    void operator()(T, Graph& g, boost::any& aprop, uint8_t val, bool ignore,
                    bool& found, std::istream& s) const
    {
        typedef typename mpl::find<val_types, T>::type pos;
        if (mpl::distance<typename mpl::begin<val_types>::type, pos>::type::value == val)
        {
            typedef typename property_map_type::apply<T, typename RangeTraits::index_map_t>::type pmap_t;
            pmap_t prop(RangeTraits::get_index_map(g));
            if (!ignore)
            {
                for (auto x : RangeTraits::get_range(g))
                    read<BE>(s, prop[x]);
                aprop = prop;
            }
            else
            {
                T y;
                for (auto x : RangeTraits::get_range(g))
                {
                    (void)x;
                    skip<BE>(s, y);
                }
            }
            found = true;
        }
    }
};

template <bool BE, class RangeTraits, class Graph>
std::pair<std::string, boost::any>
read_property(Graph& g, const std::unordered_set<std::string>& ignore,
              std::istream& s)
{
    boost::any prop;
    bool found = false;
    std::string name;
    read<BE>(s, name);
    bool skip = ignore.find(name) != ignore.end();
    uint8_t val = 0;
    read<BE>(s, val);
    mpl::for_each<val_types>(std::bind(read_property_dispatch<BE, RangeTraits>(),
                                       std::placeholders::_1, std::ref(g),
                                       std::ref(prop), val, skip, std::ref(found),
                                       std::ref(s)));
    if (!found)
        throw IOException("Error reading graph: invalid property value type index "
                          + boost::lexical_cast<std::string>(val));
    return make_pair(name, prop);
}


template <class Graph, class VProp>
void write_graph(Graph& g, const VProp& vindex, size_t N, bool directed,
                 std::vector<std::pair<std::string, boost::any>>& gprops,
                 std::vector<std::pair<std::string, boost::any>>& vprops,
                 std::vector<std::pair<std::string, boost::any>>& eprops, std::ostream& s)
{
    s.write(_magic, _magic_length);
    write(s, _version);
    uint8_t big_end = is_bigendian();
    write(s, big_end);
    string comment = "graph-tool binary file (http:://graph-tool.skewed.de)"
        " generated by version " VERSION " (commit " GIT_COMMIT ", " GIT_COMMIT_DATE ")";
    comment += " stats: " + lexical_cast<std::string>(N) + " vertices, " +
        lexical_cast<std::string>(num_edges(g)) + " edges, " +
        std::string((directed) ? "directed, " : "undirected, ") +
        lexical_cast<std::string>(gprops.size()) + " graph props, " +
        lexical_cast<std::string>(vprops.size()) + " vertex props, " +
        lexical_cast<std::string>(eprops.size()) + " edge props";
    write(s, comment);

    write_adjacency(g, vindex, N, directed, s);
    uint64_t nprops = gprops.size() + vprops.size() + eprops.size();
    write(s, nprops);
    for (auto& p : gprops)
        write_property<graph_range_traits>(g, p.first, p.second, s);
    for (auto& p : vprops)
        write_property<vertex_range_traits>(g, p.first, p.second, s);
    for (auto& p : eprops)
        write_property<edge_range_traits>(g, p.first, p.second, s);
}

template <bool BE, class Graph>
bool read_graph_dispatch(Graph& g,
                         std::vector<std::pair<std::string, boost::any>>& gprops,
                         std::vector<std::pair<std::string, boost::any>>& vprops,
                         std::vector<std::pair<std::string, boost::any>>& eprops,
                         const std::unordered_set<std::string>& ignore_gp,
                         const std::unordered_set<std::string>& ignore_vp,
                         const std::unordered_set<std::string>& ignore_ep,
                         std::istream& s)
{
    bool directed = read_adjacency<BE>(g, s);
    uint64_t nprops;
    read<BE>(s, nprops);
    for (size_t i = 0; i < nprops; ++i)
    {
        property_type pt;
        read<BE>(s, pt);
        std::pair<std::string, boost::any> p;
        switch (pt)
        {
        case property_type::Graph:
            p = read_property<BE, graph_range_traits>(g, ignore_gp, s);
            if (!p.second.empty())
                gprops.push_back(p);
            break;
        case property_type::Vertex:
            p = read_property<BE, vertex_range_traits>(g, ignore_vp, s);
            if (!p.second.empty())
                vprops.push_back(p);
            break;
        case property_type::Edge:
            p = read_property<BE, edge_range_traits>(g, ignore_ep, s);
            if (!p.second.empty())
                eprops.push_back(p);
            break;
        default:
            throw IOException("Error reading graph: invalid property type " +
                              boost::lexical_cast<std::string>(uint8_t(pt)));
        }
    }
    return directed;
}


template <class Graph>
bool read_graph(std::istream& s, Graph& g,
                std::vector<std::pair<std::string, boost::any>>& gprops,
                std::vector<std::pair<std::string, boost::any>>& vprops,
                std::vector<std::pair<std::string, boost::any>>& eprops,
                const std::unordered_set<std::string>& ignore_gp = std::unordered_set<std::string>(),
                const std::unordered_set<std::string>& ignore_vp = std::unordered_set<std::string>(),
                const std::unordered_set<std::string>& ignore_ep = std::unordered_set<std::string>())
{
    char magic[_magic_length];
    s.read(magic, _magic_length);
    if (strncmp(magic, _magic, _magic_length) != 0)
        throw IOException("Error reading graph: Invalid magic number");
    uint8_t version = 0;
    read<false>(s, version);
    if (version != _version)
        throw IOException("Error reading graph: Invalid format version " +
                          boost::lexical_cast<std::string>(version));
    uint8_t big_end = 0;
    read<false>(s, big_end);
    string comment;
    read<false>(s, comment);

    if (big_end)
        return read_graph_dispatch<true>(g, gprops, vprops, eprops, ignore_gp,
                                         ignore_vp, ignore_ep, s);
    else
        return read_graph_dispatch<false>(g, gprops, vprops, eprops, ignore_gp,
                                          ignore_vp, ignore_ep, s);
}

} // namespace graph_tool

#endif // GRAPH_IO_BINARY_HH
