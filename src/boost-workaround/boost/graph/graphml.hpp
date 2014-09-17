// Copyright (C) 2006  Tiago de Paula Peixoto <tiago@skewed.de>
// Copyright (C) 2004  The Trustees of Indiana University.
//
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
//  Authors: Douglas Gregor
//           Andrew Lumsdaine
//           Tiago de Paula Peixoto

#ifndef BOOST_GRAPH_GRAPHML_HPP
#define BOOST_GRAPH_GRAPHML_HPP

#include <boost/config.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/any.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/graph/graphviz.hpp> // for exceptions
#include <typeinfo>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/find.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/python/object.hpp>
#include <boost/bind.hpp>
#include <string>
#include <exception>
#include <set>
#include <unordered_set>

namespace boost
{

// Base64 Encoding

std::string base64_encode(const std::string& s);
std::string base64_decode(const std::string& s);


/////////////////////////////////////////////////////////////////////////////
// Graph reader exceptions
/////////////////////////////////////////////////////////////////////////////
struct parse_error: public graph_exception
{
    parse_error(const std::string& err)
    {
        error = err;
        statement = "parse error: " + error;
    }
    virtual ~parse_error() throw() {}
    virtual const char* what() const throw() {return statement.c_str();}
    std::string statement;
    std::string error;
};

class mutate_graph
{
public:
    virtual ~mutate_graph() {}
    virtual int is_directed() const = 0;
    virtual void flip_directed(bool) = 0;
    virtual bool get_directed() const = 0;

    virtual boost::any do_add_vertex() = 0;
    virtual std::pair<boost::any,bool> do_add_edge(boost::any source,
                                                   boost::any target) = 0;

    virtual size_t n_vertices() const = 0;

    virtual void
    set_graph_property(const std::string& name, const std::string& value,
                       const std::string& value_type) = 0;

    virtual void
    set_vertex_property(const std::string& name, boost::any vertex,
                        const std::string& value,
                        const std::string& value_type) = 0;

    virtual void
    set_edge_property(const std::string& name, boost::any edge,
                      const std::string& value,
                      const std::string& value_type) = 0;
};

template <typename MutableGraph>
class mutate_graph_impl : public mutate_graph
{
    typedef typename graph_traits<MutableGraph>::vertex_descriptor
        vertex_descriptor;
    typedef typename graph_traits<MutableGraph>::edge_descriptor
        edge_descriptor;

public:
    mutate_graph_impl(MutableGraph& g, dynamic_properties& dp,
                      bool ignore_directedness,
                      std::unordered_set<std::string> ignore_vp,
                      std::unordered_set<std::string> ignore_ep,
                      std::unordered_set<std::string> ignore_gp)
        : m_g(g), m_dp(dp), m_ignore_directedness(ignore_directedness),
          m_is_directed(false), m_ignore_vp(ignore_vp),
          m_ignore_ep(ignore_ep), m_ignore_gp(ignore_gp) { }

    virtual int is_directed() const
    {
        if (m_ignore_directedness)
            return 2;
        return is_convertible
            <typename graph_traits<MutableGraph>::directed_category,
            directed_tag>::value;
    }

    virtual void flip_directed(bool directed)
    {
        if (!m_is_directed)
            m_is_directed = directed;
    }

    virtual bool get_directed() const
    {
        return m_is_directed;
    }

    virtual any do_add_vertex()
    {
        return any(add_vertex(m_g));
    }

    virtual std::pair<any,bool> do_add_edge(any source, any target)
    {
        std::pair<edge_descriptor,bool> retval =
            add_edge(any_cast<vertex_descriptor>(source),
                     any_cast<vertex_descriptor>(target), m_g);
        return std::make_pair(any(retval.first), retval.second);
    }

    virtual size_t n_vertices() const
    {
        return num_vertices(m_g);
    }

    virtual void
    set_graph_property(const std::string& name,
                       const std::string& value, const std::string& value_type)
    {
        if (m_ignore_gp.find(name) != m_ignore_gp.end())
            return;

        bool type_found = false;
        try
        {
            mpl::for_each<value_types>
                (put_property<graph_property_tag,value_types>
                 (name, m_dp, graph_property_tag(), value, value_type,
                  m_type_names, type_found));
        }
        catch (bad_lexical_cast)
        {
            throw parse_error("invalid value \"" + value + "\" for key \"" +
                              name + "\" of type \"" + value_type + "\"");
        }
        if (!type_found)
            throw  parse_error("unrecognized type \"" + value_type +
                               "\" for key " + name + "\"");

    }

    virtual void
    set_vertex_property(const std::string& name, any vertex,
                        const std::string& value, const std::string& value_type)
    {
        if (m_ignore_vp.find(name) != m_ignore_vp.end())
            return;

        bool type_found = false;
        try
        {
            mpl::for_each<value_types>
                (put_property<vertex_descriptor,value_types>
                 (name, m_dp, any_cast<vertex_descriptor>(vertex),
                  value, value_type, m_type_names, type_found));
        }
        catch (bad_lexical_cast)
        {
            throw parse_error("invalid value \"" + value + "\" for key " +
                              name + " of type " + value_type);
        }
        if (!type_found)
            throw  parse_error("unrecognized type \"" + value_type +
                               "\" for key " + name);

    }

    virtual void
    set_edge_property(const std::string& name, any edge,
                      const std::string& value, const std::string& value_type)
    {
        if (m_ignore_ep.find(name) != m_ignore_ep.end())
            return;

        bool type_found = false;
        try
        {
            mpl::for_each<value_types>
                (put_property<edge_descriptor,value_types>
                 (name, m_dp, any_cast<edge_descriptor>(edge),
                  value, value_type, m_type_names, type_found));
        }
        catch (bad_lexical_cast)
        {
            throw parse_error("invalid value \"" + value + "\" for key " +
                              name + " of type " + value_type);
        }
        if (!type_found)
            throw  parse_error("unrecognized type \"" + value_type +
                               "\" for key " + name);
    }

    template <typename Key, typename ValueVector>
    class put_property
    {
    public:
        put_property(const std::string& name, dynamic_properties& dp,
                     const Key& key,
                     const std::string& value, const std::string& value_type,
                     const char** type_names, bool& type_found)
            : m_name(name), m_dp(dp), m_key(key), m_value(value),
              m_value_type(value_type), m_type_names(type_names),
              m_type_found(type_found) {}
        template <class Value>
        void operator()(Value)
        {
            if (m_value_type ==
                m_type_names[mpl::find<ValueVector,Value>::type::pos::value])
            {
                std::string val = m_value;
                if (m_value_type == "boolean")
                {
                    if (val == "true" || val == "True")
                        val = "1";
                    if (val == "false" || val == "False")
                        val = "0";
                }

                if (is_same<Value,uint8_t>::value) // chars are stored as ints
                {
                    int v = lexical_cast<int>(val);
                    put(m_name, m_dp, m_key, uint8_t(v));
                }
                else
                {
                    if (is_same<Value, boost::python::object>::value)
                        val = base64_decode(m_value);

                    put(m_name, m_dp, m_key, lexical_cast<Value>(val));
                }
                m_type_found = true;
            }
        }
    private:
        const std::string& m_name;
        dynamic_properties& m_dp;
        const Key& m_key;
        const std::string& m_value;
        const std::string& m_value_type;
        const  char** m_type_names;
        bool& m_type_found;
    };

protected:
    MutableGraph& m_g;
    dynamic_properties& m_dp;
    bool m_ignore_directedness;
    bool m_is_directed;
    std::unordered_set<std::string> m_ignore_vp;
    std::unordered_set<std::string> m_ignore_ep;
    std::unordered_set<std::string> m_ignore_gp;
    typedef mpl::vector<uint8_t, int16_t, int32_t, int64_t, double, long double,
                        std::vector<uint8_t>, std::vector<int32_t>,
                        std::vector<int64_t>, std::vector<double>,
                        std::vector<long double>, std::vector<std::string>,
                        std::string, python::object>
        value_types;
    static const char* m_type_names[];
};

template<typename MutableGraph>
const char* mutate_graph_impl<MutableGraph>::m_type_names[] =
{"boolean", "short", "int", "long", "float", "double", "vector_boolean", "vector_int",
 "vector_long", "vector_float", "vector_double", "vector_string", "string",
 "python_object"};

void
read_graphml(std::istream& in, mutate_graph& g, bool integer_vertices,
             bool store_ids);

template<typename MutableGraph>
bool
read_graphml(std::istream& in, MutableGraph& g, dynamic_properties& dp,
             bool store_ids, bool integer_vertices, bool ignore_directedness,
             std::unordered_set<std::string> ignore_vp = std::unordered_set<std::string>(),
             std::unordered_set<std::string> ignore_ep = std::unordered_set<std::string>(),
             std::unordered_set<std::string> ignore_gp = std::unordered_set<std::string>())
{
    mutate_graph_impl<MutableGraph> mg(g, dp, ignore_directedness, ignore_vp,
                                       ignore_ep, ignore_gp);
    read_graphml(in, mg, integer_vertices, store_ids);
    return mg.get_directed();
}

template <typename Types>
class get_type_name
{
public:
    get_type_name(const std::type_info& type, const char** type_names,
                  std::string& type_name)
        : m_type(type), m_type_names(type_names), m_type_name(type_name) {}
    template <typename Type>
    void operator()(Type)
    {
        if (typeid(Type) == m_type)
            m_type_name = m_type_names[mpl::find<Types,Type>::type::pos::value];
    }
private:
    const std::type_info &m_type;
    const char** m_type_names;
    std::string &m_type_name;
};

// we will rely on lexical_cast to convert the values to strings. Thus it is
// possible to deal with correct representation issues by specializing the
// lexical_cast<> template.
struct get_string
{
    template <typename ValueType>
    void operator()(const boost::any& val, std::string& sval, ValueType) const
    {
        const ValueType* v = any_cast<ValueType>(&val);
        if (v != 0)
        {
            sval = lexical_cast<std::string>(*v);
            if (is_same<ValueType, boost::python::object>::value)
                sval = base64_encode(sval);
        }
    }
};

template <typename ValueTypes, typename Descriptor>
std::string print_value(dynamic_property_map& pmap, Descriptor v)
{
    std::string val;
    boost::any oval = pmap.get(v);
    mpl::for_each<ValueTypes>(boost::bind<void>(get_string(), boost::ref(oval),
                                                boost::ref(val), _1));
    return val;
}

std::string protect_xml_string(const std::string& s);

template <typename Graph, typename VertexIndexMap>
void
write_graphml(std::ostream& out, const Graph& g, VertexIndexMap vertex_index,
              const dynamic_properties& dp, bool ordered_vertices=false)
{
    typedef typename graph_traits<Graph>::directed_category directed_category;
    typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;

    BOOST_STATIC_CONSTANT(bool, graph_is_directed =
                          (is_convertible<directed_category*,
                                          directed_tag*>::value));

    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        << "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"\n"
           "         xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
           "         xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns"
           " http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n\n";

    typedef mpl::vector<bool, uint8_t, int8_t, uint16_t, int16_t, uint32_t,
                        int32_t, uint64_t, int64_t, float, double, long double,
                        std::vector<uint8_t>, std::vector<int32_t>,
                        std::vector<int64_t>, std::vector<double>,
                        std::vector<long double>, std::vector<std::string>,
                        std::string, python::object> value_types;
    const char* type_names[] = {"boolean", "boolean", "boolean", "short",
                                "short", "int", "int", "long", "long", "float",
                                "float", "double", "vector_boolean",
                                "vector_int", "vector_long", "vector_float",
                                "vector_double", "vector_string", "string",
                                "python_object"};

    std::map<std::string, std::string> graph_key_ids;
    std::map<std::string, std::string> vertex_key_ids;
    std::map<std::string, std::string> edge_key_ids;
    int key_count = 0;

    out << "  <!-- property keys -->\n";

    bool has_vertex_ids = false;
    bool has_edge_ids = false;

    // Output keys
    for (dynamic_properties::const_iterator i = dp.begin(); i != dp.end(); ++i)
    {
        if (i->first == "_graphml_vertex_id")
        {
            has_vertex_ids = true;
            continue;
        }

        if (i->first == "_graphml_edge_id")
        {
            has_edge_ids = true;
            continue;
        }

        std::string key_id = "key" + lexical_cast<std::string>(key_count++);
        if (i->second->key() == typeid(graph_property_tag))
            graph_key_ids[i->first] = key_id;
        else if (i->second->key() == typeid(vertex_descriptor))
            vertex_key_ids[i->first] = key_id;
        else if (i->second->key() == typeid(edge_descriptor))
            edge_key_ids[i->first] = key_id;
        else
            continue;
        std::string type_name = "string";
        mpl::for_each<value_types>
            (get_type_name<value_types>(i->second->value(), type_names,
                                        type_name));
        out << "  <key id=\"" << protect_xml_string(key_id) << "\" for=\""
            << (i->second->key() == typeid(graph_property_tag) ? "graph" :
                (i->second->key() == typeid(vertex_descriptor) ? "node" :
                 "edge")) << "\""
            << " attr.name=\"" << protect_xml_string(i->first) << "\""
            << " attr.type=\"" << protect_xml_string(type_name) << "\""
            << " />\n";
    }

    bool canonical_vertices = ordered_vertices && !has_vertex_ids;
    bool canonical_edges = !has_edge_ids;

    out << "\n  <graph id=\"G\" edgedefault=\""
        << (graph_is_directed ? "directed" : "undirected") << "\""
        << " parse.nodeids=\"" << (canonical_vertices ? "canonical" : "free")
        << "\""
        << " parse.edgeids=\"" << (canonical_edges ? "canonical" : "free")
        << "\" parse.order=\"nodesfirst\">\n\n";

    out << "   <!-- graph properties -->\n";
    // Output graph data
    for (dynamic_properties::const_iterator i = dp.begin(); i != dp.end(); ++i)
    {
        if (i->second->key() == typeid(graph_property_tag))
        {
            std::string val = protect_xml_string
                (print_value<value_types>(*i->second, graph_property_tag()));
            if (val.empty())
                continue;
            out << "   <data key=\""
                << protect_xml_string(graph_key_ids[i->first]) << "\">"
                << val << "</data>\n";
        }
    }

    out << "\n   <!-- vertices -->\n";

    typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
    vertex_iterator v, v_end;
    for (tie(v, v_end) = vertices(g); v != v_end; ++v)
    {
        out << "    <node id=\"";
        if (has_vertex_ids)
            out << protect_xml_string(get("_graphml_vertex_id", dp, *v));
        else
            out << "n" << get(vertex_index, *v);
        out << "\">\n";

        // Output data
        for (dynamic_properties::const_iterator i = dp.begin(); i != dp.end();
             ++i)
        {
            if (i->first == "_graphml_vertex_id")
                continue;

            if (i->second->key() == typeid(vertex_descriptor))
            {
                std::string val = protect_xml_string
                    (print_value<value_types>(*i->second, *v));
                if (val.empty())
                    continue;
                out << "      <data key=\""
                    << protect_xml_string(vertex_key_ids[i->first]) << "\">"
                    << val << "</data>\n";
            }
        }
        out << "    </node>\n";
    }

    out << "\n   <!-- edges -->\n";
    typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
    edge_iterator e, e_end;
    typename graph_traits<Graph>::edges_size_type edge_count = 0;
    for (tie(e, e_end) = edges(g); e != e_end; ++e)
    {
        out << "    <edge id=\"";
        if (has_edge_ids)
            out << protect_xml_string(get("_graphml_edge_id", dp, *e));
        else
            out << "e" << edge_count;
        edge_count++;

        out << "\" source=\"";
        if (has_vertex_ids)
            out << protect_xml_string(get("_graphml_vertex_id", dp,
                                           source(*e, g)));
        else
            out << "n" << get(vertex_index, source(*e, g));
        out << "\" target=\"";
        if (has_vertex_ids)
            out << protect_xml_string(get("_graphml_vertex_id", dp,
                                          target(*e, g)));
        else
            out << "n" << get(vertex_index, target(*e, g));
        out<< "\">\n";

        // Output data
        for (dynamic_properties::const_iterator i = dp.begin(); i != dp.end();
             ++i)
        {
            if (i->first == "_graphml_edge_id")
                continue;

            if (i->second->key() == typeid(edge_descriptor))
            {
                std::string val = protect_xml_string
                    (print_value<value_types>(*i->second, *e));
                if (val.empty())
                    continue;
                out << "      <data key=\""
                    << protect_xml_string(edge_key_ids[i->first]) << "\">"
                    << val << "</data>\n";
            }
        }
        out << "    </edge>\n";
    }

    out << "\n  </graph>\n"
        << "</graphml>\n";
}


template <typename Graph>
void
write_graphml(std::ostream& out, const Graph& g, const dynamic_properties& dp,
              bool ordered_vertices=false)
{
    write_graphml(out, g, get(vertex_index, g), dp, ordered_vertices);
}

} // boost namespace

#endif // BOOST_GRAPH_GRAPHML_HPP
