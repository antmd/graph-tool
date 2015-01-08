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

#ifndef GML_HH
#define GML_HH

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/variant/recursive_variant.hpp>
#include <boost/variant/get.hpp>
#include <boost/spirit/include/support_istream_iterator.hpp>
#include <boost/foreach.hpp>
#include <boost/type_traits.hpp>

#include <boost/algorithm/string/replace.hpp>

#include <boost/property_map/dynamic_property_map.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/bind/bind.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <boost/python.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <unordered_map>

namespace graph_tool{

using namespace std;
using namespace boost;

class gml_parse_error: public std::exception
{
public:
    gml_parse_error(const string& w): _what(w) {}
    ~gml_parse_error() throw() {}
    virtual const char* what() const throw() {return _what.c_str();}

private:
    std::string _what;
};


template <class Graph>
class gml_state
{
public:
    gml_state(Graph& g, dynamic_properties& dp,
              const std::unordered_set<std::string>& ignore_vp = std::unordered_set<std::string>(),
              const std::unordered_set<std::string>& ignore_ep = std::unordered_set<std::string>(),
              const std::unordered_set<std::string>& ignore_gp = std::unordered_set<std::string>())
        : _g(g), _dp(dp), _directed(false), _ignore_vp(ignore_vp),
          _ignore_ep(ignore_ep), _ignore_gp(ignore_gp) {}

    typedef boost::variant<std::string, int, double> val_t;

    // key / value mechanics
    void push_key(const std::string& key)
    {
        _stack.push_back(make_pair(key, prop_list_t()));
    }

    void push_value(const val_t& value)
    {
        if (_stack.empty())
            return;
        std::string k = _stack.back().first;
        _stack.pop_back();
        if (!_stack.empty())
            _stack.back().second[k] = value;
    }

    // actual parsing
    void finish_list()
    {
        if (_stack.empty())
            return;

        std::string &k = _stack.back().first;
        if (k == "node")
        {
            int id;
            if (_stack.back().second.find("id") == _stack.back().second.end())
                throw gml_parse_error("node does not have an id");
            try
            {
                id = boost::get<double>(_stack.back().second["id"]);
            }
            catch (bad_get)
            {
                throw gml_parse_error("invalid node id");
            }

            typename graph_traits<Graph>::vertex_descriptor v = get_vertex(id);

            // put properties
            for (typeof(_stack.back().second.begin()) iter = _stack.back().second.begin();
                 iter != _stack.back().second.end(); ++iter)
            {
                if (iter->first == "id")
                    continue;
                if (_ignore_vp.find(iter->first) != _ignore_vp.end())
                    continue;
                try
                {
                    put(iter->first, _dp, v, boost::get<string>(iter->second));
                }
                catch (bad_get)
                {
                    put(iter->first, _dp, v, boost::get<double>(iter->second));
                }
            }
        }
        else if (k == "edge")
        {
            int source, target;
            if (_stack.back().second.find("source") ==
                _stack.back().second.end() ||
                _stack.back().second.find("target") ==
                _stack.back().second.end())
                throw gml_parse_error("edge does not have source and target ids");
            try
            {
                source = boost::get<double>(_stack.back().second["source"]);
                target = boost::get<double>(_stack.back().second["target"]);
            }
            catch (bad_get)
            {
                throw gml_parse_error("invalid source and target ids");
            }

            typename graph_traits<Graph>::vertex_descriptor s, t;
            s = get_vertex(source);
            t = get_vertex(target);

            typename graph_traits<Graph>::edge_descriptor e =
                add_edge(s, t, _g).first;

            // put properties
            for (typeof(_stack.back().second.begin()) iter = _stack.back().second.begin();
                 iter != _stack.back().second.end(); ++iter)
            {
                if (iter->first == "id" || iter->first == "source" || iter->first == "target")
                    continue;
                if (_ignore_ep.find(iter->first) != _ignore_ep.end())
                    continue;
                try
                {
                    put(iter->first, _dp, e, boost::get<string>(iter->second));
                }
                catch (bad_get)
                {
                    put(iter->first, _dp, e, boost::get<double>(iter->second));
                }
            }

        }
        else if (k == "graph")
        {
            // put properties
            for (typeof(_stack.back().second.begin()) iter = _stack.back().second.begin();
                 iter != _stack.back().second.end(); ++iter)
            {
                if (iter->first == "directed")
                    _directed = boost::get<double>(iter->second);
                if (_ignore_gp.find(iter->first) != _ignore_gp.end())
                    continue;
                try
                {
                    put(iter->first, _dp, graph_property_tag(), boost::get<string>(iter->second));
                }
                catch (bad_get)
                {
                    put(iter->first, _dp, graph_property_tag(), boost::get<double>(iter->second));
                }
            }

        }
        _stack.pop_back();
    }

    typename graph_traits<Graph>::vertex_descriptor get_vertex(size_t index)
    {
        if (_vmap.find(index) == _vmap.end())
            _vmap[index] = add_vertex(_g);
        return _vmap[index];
    }


    bool is_directed()
    {
        return _directed;
    }


private:
    Graph& _g;
    dynamic_properties& _dp;
    bool _directed;
    std::unordered_map<int, typename graph_traits<Graph>::vertex_descriptor> _vmap;

    // the stack holds the keys, and its properties (but omits nested lists)
    typedef std::unordered_map<std::string, val_t> prop_list_t;
    vector<pair<std::string,  prop_list_t> > _stack;

    const std::unordered_set<std::string>& _ignore_vp;
    const std::unordered_set<std::string>& _ignore_ep;
    const std::unordered_set<std::string>& _ignore_gp;
};


template <class Iterator, class Graph, class Skipper>
struct gml : spirit::qi::grammar<Iterator, void(), Skipper>
{
    gml(Graph& g, dynamic_properties& dp,
        const std::unordered_set<std::string>& ignore_vp = std::unordered_set<std::string>(),
        const std::unordered_set<std::string>& ignore_ep = std::unordered_set<std::string>(),
        const std::unordered_set<std::string>& ignore_gp = std::unordered_set<std::string>())
        : gml::base_type(start), _state(g, dp, ignore_vp, ignore_ep, ignore_gp)
    {
        using namespace spirit;
        using spirit::ascii::char_;

        unesc_str = spirit::lexeme['"' >> *(unesc_char | (spirit::qi::char_ - "\"") | "\\x" >> qi::hex) >> '"'];
        unesc_char.add("\\a", '\a')("\\b", '\b')("\\f", '\f')("\\n", '\n')
            ("\\r", '\r')("\\t", '\t')("\\v", '\v')("\\\\", '\\')
            ("\\\'", '\'')("\\\"", '\"');
        key_identifier %= spirit::lexeme[((+spirit::qi::alnum) >> *spirit::qi::alnum)];
        key = key_identifier
            [boost::bind(&gml_state<Graph>::push_key, &_state, ::_1)];
        value_identifier %= (spirit::lexeme[spirit::qi::double_] | unesc_str);
        value %= value_identifier
            [boost::bind(&gml_state<Graph>::push_value, &_state, ::_1)];
        list_identifier = *(key >> (value | "[" >> list >> "]"));
        list = list_identifier
            [boost::bind(&gml_state<Graph>::finish_list, &_state)];
        start = list;
    }

    typedef boost::variant<std::string, double> val_t;

    spirit::qi::rule<Iterator, std::string(), Skipper> unesc_str;
    spirit::qi::symbols<char const, char const> unesc_char;
    spirit::qi::rule<Iterator, std::string(), Skipper> key, key_identifier;
    spirit::qi::rule<Iterator, val_t(), Skipper> value, value_identifier;
    spirit::qi::rule<Iterator, void(), Skipper> list, list_identifier;
    spirit::qi::rule<Iterator, void(), Skipper> start;

    gml_state<Graph> _state;
};

template <class Iterator, class Graph, class Skipper>
bool parse_grammar(Iterator begin, Iterator end, Graph& g,
                   dynamic_properties& dp, Skipper skip,
                   const std::unordered_set<std::string>& ignore_vp = std::unordered_set<std::string>(),
                   const std::unordered_set<std::string>& ignore_ep = std::unordered_set<std::string>(),
                   const std::unordered_set<std::string>& ignore_gp = std::unordered_set<std::string>())
{
    using namespace spirit;
    gml<spirit::istream_iterator, Graph, Skipper> parser(g, dp, ignore_vp,
                                                         ignore_ep, ignore_gp);
    bool ok = qi::phrase_parse(begin, end, parser, skip);
    if (!ok)
        throw gml_parse_error("invalid syntax");
    return parser._state.is_directed();
}


template <class Graph>
bool read_gml(istream& in, Graph& g, dynamic_properties& dp,
              const std::unordered_set<std::string>& ignore_vp = std::unordered_set<std::string>(),
              const std::unordered_set<std::string>& ignore_ep = std::unordered_set<std::string>(),
              const std::unordered_set<std::string>& ignore_gp = std::unordered_set<std::string>())
{
    using namespace spirit;

    in >> std::noskipws;
    spirit::istream_iterator begin(in);
    spirit::istream_iterator end;

    bool directed =
        parse_grammar(begin, end, g, dp,
                      (ascii::space |'#' >> *(ascii::char_ - qi::eol) >> qi::eol),
                      ignore_vp, ignore_ep, ignore_gp);

    return directed;
}

struct get_str
{
    template <typename ValueType>
    void operator()(const boost::any& val, std::string& sval, ValueType) const
    {
        try
        {
            ValueType v = any_cast<ValueType>(val);
            if (std::is_same<ValueType, python::object>::value)
            {
                sval = lexical_cast<string>(v);
            }
            else
            {
                stringstream s;
                s << v;
                sval = s.str();
            }

            if (!std::is_scalar<ValueType>::value)
            {
                replace_all(sval, "\"", "\\\"");
                sval = "\"" + sval + "\"";
            }
        }
        catch (bad_any_cast)
        {
        }
    }
};

template <typename ValueTypes, typename Descriptor>
std::string print_val(dynamic_property_map& pmap, const Descriptor& v)
{
    std::string val;
    boost::any oval = pmap.get(v);
    mpl::for_each<ValueTypes>(bind<void>(get_str(), boost::ref(oval),
                                         boost::ref(val), _1));
    return val;
}


template <typename Graph, typename VertexIndexMap>
void write_gml(std::ostream& out, const Graph& g, VertexIndexMap vertex_index,
               const dynamic_properties& dp)
{
    typedef typename graph_traits<Graph>::directed_category directed_category;
    typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;

    typedef mpl::vector<bool, uint8_t, int8_t, uint32_t, int32_t,
                        uint64_t, int64_t, float, double, long double,
                        std::vector<uint8_t>, std::vector<int32_t>,
                        std::vector<int64_t>, std::vector<double>,
                        std::vector<long double>, std::vector<std::string>,
                        std::string, python::object> value_types;

    BOOST_STATIC_CONSTANT(bool, graph_is_directed =
                          (std::is_convertible<directed_category*,
                                               directed_tag*>::value));

    out << "graph [" << endl;

    if (graph_is_directed)
        out << "   directed " << 1 << endl;

    for (dynamic_properties::const_iterator i = dp.begin(); i != dp.end();
         ++i)
    {
        if (i->second->key() == typeid(graph_property_tag))
        {
            std::string val = print_val<value_types>(*i->second,
                                                     graph_property_tag());
            if (val.empty())
                continue;
            out << "   " << i->first << " " << val << endl;
        }
    }

    typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
    vertex_iterator v, v_end;
    for (tie(v, v_end) = vertices(g); v != v_end; ++v)
    {
        out << "   node [" << endl;
        out << "      id " << get(vertex_index, *v) << endl;

        for (dynamic_properties::const_iterator i = dp.begin(); i != dp.end();
             ++i)
        {
            if (i->second->key() == typeid(vertex_descriptor))
            {
                std::string val = print_val<value_types>(*i->second, *v);
                if (val.empty())
                    continue;
                out << "      " << i->first << " " << val << endl;
            }
        }
        out << "   ]" << endl;
    }

    typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
    edge_iterator e, e_end;
    typename graph_traits<Graph>::edges_size_type edge_count = 0;
    for (tie(e, e_end) = edges(g); e != e_end; ++e)
    {
        out << "   edge [" << endl;
        out << "      id " << edge_count++ << endl;
        out << "      source " << get(vertex_index, source(*e, g)) << endl;
        out << "      target " << get(vertex_index, target(*e, g)) << endl;

        for (dynamic_properties::const_iterator i = dp.begin(); i != dp.end();
             ++i)
        {
            if (i->second->key() == typeid(edge_descriptor))
            {
                std::string val = print_val<value_types>(*i->second, *e);
                if (val.empty())
                    continue;
                out << "      " << i->first << " " << val << endl;
            }
        }
        out << "   ]" << endl;
    }
    out << "]" << endl;
}


template <typename Graph>
void write_gml(std::ostream& out, const Graph& g, const dynamic_properties& dp)
{
    write_gml(out, g, get(vertex_index, g), dp);
}

} // namespace graph_tool

#endif // GML_HH
