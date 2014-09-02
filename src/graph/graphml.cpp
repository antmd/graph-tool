// Copyright (C) 2006  Tiago de Paula Peixoto <tiago@skewed.de>
// Copyright (C) 2004  The Trustees of Indiana University.
//
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Douglas Gregor
//           Andrew Lumsdaine
//           Tiago de Paula Peixoto

#include <boost/python.hpp>
#include <boost/variant.hpp>
#include <expat.h>
#include <boost/graph/graphml.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/archive/iterators/xml_escape.hpp>
#include <boost/archive/iterators/ostream_iterator.hpp>
#include <boost/archive/iterators/base64_from_binary.hpp>
#include <boost/archive/iterators/binary_from_base64.hpp>
#include <boost/archive/iterators/insert_linebreaks.hpp>
#include <boost/archive/iterators/transform_width.hpp>
#include <boost/archive/iterators/ostream_iterator.hpp>
#include <sstream>

using namespace boost;
namespace boost
{

std::string base64_encode(const std::string& s)
{
    static const std::string base64_padding[] = {"", "==","="};
    namespace bai = boost::archive::iterators;
    std::stringstream os;
    typedef bai::base64_from_binary<bai::transform_width<const char *, 6, 8> >
        base64_enc;
    std::copy(base64_enc(s.c_str()), base64_enc(s.c_str() + s.size()),
              std::ostream_iterator<char>(os));
    os << base64_padding[s.size() % 3];
    return os.str();
}

std::string base64_decode(const std::string& s)
{
    namespace bai = boost::archive::iterators;
    std::stringstream os;

    typedef bai::transform_width<bai::binary_from_base64<const char *>, 8, 6> base64_dec;

    unsigned int size = s.size();

    // Remove the padding characters, cf. https://svn.boost.org/trac/boost/ticket/5629
    if (size && s[size - 1] == '=')
    {
        --size;
        if (size && s[size - 1] == '=')
            --size;
    }
    if (size == 0)
        return std::string();

    std::copy(base64_dec(s.data()), base64_dec(s.data() + size),
              std::ostream_iterator<char>(os));

    return os.str();
}


std::string protect_xml_string(const std::string& os)
{
    using namespace boost::archive::iterators;
    std::stringstream s;
    std::copy(xml_escape<const char*>(os.c_str()),
              xml_escape<const char*>(os.c_str()+os.size()),
              ostream_iterator<char>(s));
    return s.str();
}
}

class graphml_reader
{
public:
    graphml_reader(mutate_graph& g, bool integer_vertices, bool store_ids)
        : m_g(g),
          m_canonical_vertices(false),
          m_canonical_edges(false),
          m_integer_vertices(integer_vertices),
          m_store_ids(store_ids) { }

    void run(std::istream& in)
    {
        const int buffer_size = 4096;
        m_parser = XML_ParserCreateNS(0,'|');
        XML_SetElementHandler(m_parser, &on_start_element, &on_end_element);
        XML_SetCharacterDataHandler(m_parser, &on_character_data);
        XML_SetUserData(m_parser, this);
        char buffer[buffer_size];

        bool okay = true;
        do
        {
            in.read(buffer, buffer_size);
            okay = XML_Parse(m_parser, buffer, in.gcount(), in.gcount() == 0);
        }
        while (okay && in.good());

        if (!okay)
        {
            std::stringstream s;
            s << "on line " << XML_GetCurrentLineNumber(m_parser)
              <<", column " << XML_GetCurrentColumnNumber(m_parser)
              << ": " << XML_ErrorString(XML_GetErrorCode(m_parser));
            throw parse_error(s.str());
        }
        XML_ParserFree(m_parser);
    }

private:
    /// The kinds of keys. Not all of these are supported
    enum key_kind {
        graph_key,
        node_key,
        edge_key,
        hyperedge_key,
        port_key,
        endpoint_key,
        all_key
    };

    enum desc_kind
    {
        M_VERTEX_DESCRIPTOR,
        M_EDGE_DESCRIPTOR,
        M_GRAPH_DESCRIPTOR
    };


    static void
    on_start_element(void* user_data, const XML_Char *c_name,
                     const XML_Char **atts)
    {
        graphml_reader* self = static_cast<graphml_reader*>(user_data);

        std::string name(c_name);
        replace_first(name, "http://graphml.graphdrawing.org/xmlns|", "");

        if (name == "edge")
        {
            std::string id;
            std::string source, target;
            while (*atts)
            {
                std::string name = *atts++;
                std::string value = *atts++;

                if (name == "id") id = value;
                else if (name == "source") source = value;
                else if (name == "target") target = value;
                else if (name == "directed")
                {
                    bool edge_is_directed = (value == "directed");
                    if (self->m_g.is_directed() != 2 &&
                        edge_is_directed != self->m_g.is_directed())
                    {
                        if (edge_is_directed)
                            throw directed_graph_error();
                        else
                            throw undirected_graph_error();
                    }
                    self->m_g.flip_directed(edge_is_directed);
                }
            }

            self->m_active_descriptor = self->handle_edge(id, source, target);
            self->m_descriptor_kind = M_EDGE_DESCRIPTOR;
        }
        else if (name == "node")
        {
            std::string id;

            while (*atts)
            {
                std::string name = *atts++;
                std::string value = *atts++;

                if (name == "id") id = value;
            }

            self->m_active_descriptor = self->handle_vertex(id);
            self->m_descriptor_kind = M_VERTEX_DESCRIPTOR;
        }
        else if (name == "data")
        {
            while (*atts)
            {
                std::string name = *atts++;
                std::string value = *atts++;

                if (name == "key") self->m_active_key = value;
            }
        }
        else if (name == "key")
        {
            std::string id;
            std::string key_name;
            std::string key_type;
            key_kind kind = all_key;

            while (*atts)
            {
                std::string name = *atts++;
                std::string value = *atts++;

                if (name == "id") id = value;
                else if (name == "attr.name") key_name = value;
                else if (name == "attr.type") key_type = value;
                else if (name == "for")
                {
                    if (value == "graph") kind = graph_key;
                    else if (value == "node") kind = node_key;
                    else if (value == "edge") kind = edge_key;
                    else if (value == "hyperedge") kind = hyperedge_key;
                    else if (value == "port") kind = port_key;
                    else if (value == "endpoint") kind = endpoint_key;
                    else if (value == "all") kind = all_key;
                    else
                    {
                        std::stringstream s;
                        s << "on line "
                          << XML_GetCurrentLineNumber(self->m_parser)
                          << ", column "
                          << XML_GetCurrentColumnNumber(self->m_parser)
                          << ": unrecognized key kind '" << value << "'";
                        throw parse_error(s.str());
                    }
                }
            }

            self->m_keys[id] = kind;
            self->m_key_name[id] = key_name;
            self->m_key_type[id] = key_type;
            self->m_active_key = id;
        }
        else if (name == "graph")
        {
            while (*atts)
            {
                std::string name = *atts++;
                std::string value = *atts++;

                if (name == "edgedefault")
                {
                    bool edge_is_directed = (value == "directed");
                    if ( self->m_g.is_directed() != 2 &&
                         edge_is_directed != self->m_g.is_directed())
                    {
                        if (edge_is_directed)
                            throw directed_graph_error();
                        else
                            throw undirected_graph_error();
                    }
                    self->m_g.flip_directed(edge_is_directed);
                }
                else if (name == "parse.nodeids")
                {
                    self->m_canonical_vertices = (value == "canonical");
                }
                else if (name == "parse.edgeids")
                {
                    self->m_canonical_edges = (value == "canonical");
                }
            }
            self->m_active_descriptor = any();
            self->m_descriptor_kind = M_GRAPH_DESCRIPTOR;
        }

        self->m_character_data.clear();
    }

    static void
    on_end_element(void* user_data, const XML_Char *c_name)
    {
        graphml_reader* self = static_cast<graphml_reader*>(user_data);

        std::string name(c_name);
        replace_first(name, "http://graphml.graphdrawing.org/xmlns|", "");

        if (name == "data")
        {
            switch (self->m_descriptor_kind)
            {
            case M_VERTEX_DESCRIPTOR:
                self->handle_vertex_property(self->m_active_key, self->m_active_descriptor,
                                             self->m_character_data);
                break;
            case M_EDGE_DESCRIPTOR:
                self->handle_edge_property(self->m_active_key, self->m_active_descriptor,
                                           self->m_character_data);
                break;
            case M_GRAPH_DESCRIPTOR:
                self->handle_graph_property(self->m_active_key, self->m_active_descriptor,
                                             self->m_character_data);
                break;
            }
        }
        else if (name == "default")
        {
            self->m_key_default[self->m_active_key] = self->m_character_data;
        }
    }

    static void
    on_character_data(void* user_data, const XML_Char* s, int len)
    {
        graphml_reader* self = static_cast<graphml_reader*>(user_data);
        self->m_character_data.append(s, len);
    }

    any
    handle_vertex(const std::string& v)
    {
        bool is_new = false;

        if (m_canonical_vertices)
        {
            size_t id;

            //strip leading "n" from name
            try
            {
                id = lexical_cast<size_t>(std::string(v,1));
            }
            catch (bad_lexical_cast)
            {
                std::stringstream s;
                s << "on line " << XML_GetCurrentLineNumber(m_parser)
                  << ", column " << XML_GetCurrentColumnNumber(m_parser)
                  << ": invalid vertex: " << v;
                throw parse_error(s.str());
            }

            if (m_integer_vertices)
            {
                is_new = (m_g.n_vertices() <= id);
                for (size_t i = m_g.n_vertices(); i <= id; ++i)
                    m_g.do_add_vertex();
            }
            else
            {
                while(id >= m_canonical_vertex.size())
                {
                    m_canonical_vertex.push_back(m_g.do_add_vertex());
                    is_new = true;
                }
            }
        }
        else
        {
            if (m_vertex.find(v) == m_vertex.end())
            {
                m_vertex[v] = m_g.do_add_vertex();
                is_new = true;
            }
        }

        any vd = get_vertex_descriptor(v);
        if (is_new)
        {
            std::map<std::string, std::string>::iterator iter;
            for (iter = m_key_default.begin(); iter != m_key_default.end();
                 ++iter)
            {
                if (m_keys[iter->first] == node_key)
                    handle_vertex_property(iter->first, vd, iter->second);
            }
            if (m_store_ids && !m_canonical_vertices)
                m_g.set_vertex_property("_graphml_vertex_id",
                                        vd, v, "string");
        }
        return vd;
    }

    any
    get_vertex_descriptor(const std::string& v)
    {
        if (m_canonical_vertices)
        {
            //strip leading "n" from name
            size_t id = lexical_cast<size_t>(std::string(v,1));
            if (m_integer_vertices)
                return id;
            return m_canonical_vertex[id];
        }
        else
        {
            return m_vertex[v];
        }
    }

    any
    handle_edge(const std::string& id, const std::string& u,
                const std::string& v)
    {
        handle_vertex(u);
        handle_vertex(v);

        any source, target;
        source = get_vertex_descriptor(u);
        target = get_vertex_descriptor(v);

        any edge;
        bool added;
        tie(edge, added) = m_g.do_add_edge(source, target);
        if (!added)
            throw bad_parallel_edge(u, v);

        std::map<std::string, std::string>::iterator iter;
        for (iter = m_key_default.begin(); iter != m_key_default.end(); ++iter)
        {
            if (m_keys[iter->first] == edge_key)
                handle_edge_property(iter->first, edge, iter->second);
        }

        if (m_store_ids && !m_canonical_edges)
            m_g.set_edge_property("_graphml_edge_id", edge, id, "string");

        return edge;
    }

    void handle_edge_property(const std::string& key_id,
                              const any& descriptor,
                              const std::string& value)
    {
        try
        {
            m_g.set_edge_property(m_key_name[key_id], descriptor,
                              value, m_key_type[key_id]);
        }
        catch (parse_error &e)
        {
            std::stringstream s;
            s << "on line " << XML_GetCurrentLineNumber(m_parser)
              << ", column " << XML_GetCurrentColumnNumber(m_parser)
              << ": " << e.error;
            throw parse_error(s.str());
        }
    }

    void handle_vertex_property(const std::string& key_id,
                                const any& descriptor,
                                const std::string& value)
    {
        try
        {
            m_g.set_vertex_property(m_key_name[key_id], descriptor,
                                    value, m_key_type[key_id]);
        }
        catch (parse_error &e)
        {
            std::stringstream s;
            s << "on line " << XML_GetCurrentLineNumber(m_parser)
              << ", column " << XML_GetCurrentColumnNumber(m_parser)
              << ": " << e.error;
            throw parse_error(s.str());
        }
    }

    void handle_graph_property(const std::string& key_id,
                               const any&,
                               const std::string& value)
    {
        try
        {
            m_g.set_graph_property(m_key_name[key_id], value,
                                   m_key_type[key_id]);
        }
        catch (parse_error &e)
        {
            std::stringstream s;
            s << "on line " << XML_GetCurrentLineNumber(m_parser)
              << ", column " << XML_GetCurrentColumnNumber(m_parser)
              << ": " << e.error;
            throw parse_error(s.str());
        }
    }

    mutate_graph& m_g;
    std::map<std::string, key_kind> m_keys;
    std::map<std::string, std::string> m_key_name;
    std::map<std::string, std::string> m_key_type;
    std::map<std::string, std::string> m_key_default;
    std::map<std::string, any> m_vertex;
    std::vector<any> m_canonical_vertex;

    any m_active_descriptor;
    desc_kind m_descriptor_kind;

    std::string m_active_key;
    std::string m_character_data;
    bool m_canonical_vertices;
    bool m_canonical_edges;
    bool m_integer_vertices;
    bool m_store_ids;
    bool m_ignore_directedness;
    XML_Parser m_parser;
};

namespace boost
{
void
read_graphml(std::istream& in, mutate_graph& g, bool integer_vertices, bool store_ids)
{
    graphml_reader reader(g, integer_vertices, store_ids);
    reader.run(in);
}
}
