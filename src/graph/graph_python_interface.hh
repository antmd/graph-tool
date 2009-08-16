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

#ifndef PYTHON_INTERFACE_HH
#define PYTHON_INTERFACE_HH

#include <boost/lambda/bind.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/functional/hash.hpp>
#include <boost/iterator/iterator_facade.hpp>

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_selectors.hh"
#include "numpy_bind.hh"

#include <boost/python.hpp>
#include <boost/python/type_id.hpp>

// this file includes a simple python interface for the internally kept
// graph. It defines a PythonVertex, PythonEdge and PythonIterator template
// classes, which contain the proper member functions for graph traversal. These
// types are then specialized for each version of the adapted graph (directed,
// undirected, filtered, reversed).

namespace graph_tool
{
using namespace boost;

// generic iterator adaptor which can be used to iterate vertices, edges,
// out_edges and in_edges through python
template <class Descriptor, class Iterator>
class PythonIterator
{
public:
    PythonIterator(GraphInterface& gi, std::pair<Iterator,Iterator> e)
        : _gi(gi), _e(e) {}
    Descriptor Next()
    {
        if (_e.first == _e.second)
            python::objects::stop_iteration_error();
        Descriptor e(_gi, *_e.first);
        ++_e.first;
        return e;
    }
private:
    GraphInterface& _gi;
    std::pair<Iterator,Iterator> _e;
};


// forward declaration of PythonEdge
template <class Graph>
class PythonEdge;

// below are classes related to the PythonVertex type
class PythonVertex
{
public:
    PythonVertex(GraphInterface& gi, GraphInterface::vertex_t v):
        _gi(gi), _v(v), _valid(true)
    {
        CheckValid();
    }

    bool IsValid() const
    {
        return _valid &&
            (_v != graph_traits<GraphInterface::multigraph_t>::null_vertex()) &&
            (_v < num_vertices(_gi._mg));
    }

    void SetValid(bool valid)
    {
        _valid = valid;
    }

    void CheckValid() const
    {
        if (!IsValid())
            throw ValueException("invalid vertex descriptor");
    }

    GraphInterface::vertex_t GetDescriptor() const
    {
        return _v;
    }

    template <class DegSelector>
    struct get_degree
    {
        template<class Graph>
        void operator()(const Graph& g,
                        typename graph_traits<Graph>::vertex_descriptor v,
                        size_t& deg) const
        {
            deg = DegSelector()(v, g);
        }
    };

    size_t GetInDegree() const
    {
        CheckValid();
        size_t in_deg;
        run_action<>()(_gi,lambda::bind<void>(get_degree<in_degreeS>(),
                                              lambda::_1, _v,
                                              lambda::var(in_deg)))();
        return in_deg;
    }

    size_t GetOutDegree() const
    {
        CheckValid();
        size_t out_deg;
        run_action<>()(_gi,lambda::bind<void>(get_degree<out_degreeS>(),
                                              lambda::_1, _v,
                                              lambda::var(out_deg)))();
        return out_deg;
    }

    // provide iterator support for out_edges

    struct get_out_edges
    {
        template<class Graph>
        void operator()(const Graph& g, GraphInterface& gi,
                        typename graph_traits<Graph>::vertex_descriptor v,
                        python::object& iter) const
        {
            typedef typename graph_traits<Graph>::out_edge_iterator
                out_edge_iterator;
            iter = python::object(PythonIterator<PythonEdge<Graph>,
                                                 out_edge_iterator>
                                  (gi, out_edges(v, g)));
        }
    };

    python::object
    OutEdges() const
    {
        CheckValid();
        python::object iter;
        run_action<>()(_gi, lambda::bind<void>(get_out_edges(), lambda::_1,
                                               lambda::var(_gi), _v,
                                               lambda::var(iter)))();
        return iter;
    }

    struct get_in_edges
    {
        template<class Graph>
        void operator()(const Graph& g, GraphInterface& gi,
                        typename graph_traits<Graph>::vertex_descriptor v,
                        python::object& iter) const
        {
            typedef typename in_edge_iteratorS<Graph>::type
                in_edge_iterator;
            iter = python::object
                (PythonIterator<PythonEdge<Graph>,in_edge_iterator>
                 (gi, in_edge_iteratorS<Graph>::get_edges(v, g)));
        }
    };

    python::object
    InEdges() const
    {
        CheckValid();
        python::object iter;
        run_action<>()(_gi, lambda::bind<void>(get_in_edges(), lambda::_1,
                                               lambda::var(_gi), _v,
                                               lambda::var(iter)))();
        return iter;
    }

    std::string GetString() const
    {
        CheckValid();
        return lexical_cast<std::string>(_v);
    }

    size_t GetHash() const
    {
        return hash<size_t>()(_gi._vertex_index[_v]);
    }

    size_t GetIndex() const
    {
        return _gi._vertex_index[_v];
    }

    bool operator==(const PythonVertex& other) const
    {
        CheckValid();
        other.CheckValid();
        return other._v == _v;
    }

    bool operator!=(const PythonVertex& other) const
    {
        CheckValid();
        other.CheckValid();
        return other._v != _v;
    }

private:
    GraphInterface& _gi;
    GraphInterface::vertex_t _v;
    bool _valid;
};

// below are classes related to the PythonEdge type

template <class Graph>
class PythonEdge
{
public:
    typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
    PythonEdge(GraphInterface& gi, edge_descriptor e)
        : _gi(gi), _e(e), _valid(true)
    {
        CheckValid();
    }

    bool IsValid() const
    {
        return _valid &&
            (_gi._vertex_index[target(_e, _gi._mg)] < num_vertices(_gi._mg)) &&
            (_gi._vertex_index[source(_e, _gi._mg)] < num_vertices(_gi._mg)) &&
            (target(_e, _gi._mg) != graph_traits<Graph>::null_vertex()) &&
            (source(_e, _gi._mg) != graph_traits<Graph>::null_vertex());
    }

    void SetValid(bool valid)
    {
        _valid = valid;
    }

    void CheckValid() const
    {
        if (!_valid)
            throw ValueException("invalid edge descriptor");
    }

    GraphInterface::edge_t GetDescriptor() const
    {
        return _e;
    }

    struct get_source
    {
        template<class GraphType>
        void operator()(const GraphType& g, GraphInterface& gi,
                        const edge_descriptor& edge, python::object& vertex)
            const
        {
            vertex = python::object(PythonVertex(gi, source(edge, g)));
        }
    };

    python::object GetSource() const
    {
        CheckValid();
        python::object v;
        run_action<>()(_gi, lambda::bind<void>(get_source(), lambda::_1,
                                               lambda::var(_gi),
                                               lambda::var(_e),
                                               lambda::var(v)))();
        return v;
    }

    struct get_target
    {
        template<class GraphType>
        void operator()(const GraphType& g, GraphInterface& gi,
                        const edge_descriptor& edge,
                        python::object& vertex) const
        {
            vertex = python::object(PythonVertex(gi, target(edge, g)));
        }
    };

    python::object GetTarget() const
    {
        CheckValid();
        python::object v;
        run_action<>()(_gi, lambda::bind<void>(get_target(), lambda::_1,
                                               lambda::var(_gi),
                                               lambda::var(_e),
                                               lambda::var(v)))();
        return v;
    }

    std::string GetString() const
    {
        CheckValid();
        PythonVertex src = python::extract<PythonVertex>(GetSource());
        PythonVertex tgt = python::extract<PythonVertex>(GetTarget());
        return "(" + src.GetString() + "," + tgt.GetString() + ")";
    }

    size_t GetHash() const
    {
        CheckValid();
        return hash<size_t>()(_gi._edge_index[_e]);
    }

    bool operator==(const PythonEdge& other) const
    {
        CheckValid();
        return other._e == _e;
    }

    bool operator!=(const PythonEdge& other) const
    {
        CheckValid();
        return other._e != _e;
    }

private:
    GraphInterface& _gi;
    edge_descriptor _e;
    bool _valid;
};

// metafunction to determine wether or not to return copies or internal
// references to property types
struct return_reference
{
    template <class ValueType>
    struct apply
    {
        // return actual references only for non-string and non-python::object
        // classes

        typedef typename mpl::if_<
            typename mpl::and_<
                is_class<ValueType>,
                typename mpl::and_<
                    typename mpl::not_<is_same<ValueType,
                                               string> >::type,
                    typename mpl::not_<is_same<ValueType,
                                               python::object> >::type>::type
                >::type,
            mpl::bool_<true>,
            mpl::bool_<false> >::type type;
    };
};

template <class PropertyMap>
class PythonPropertyMap
{
public:
    PythonPropertyMap(const PropertyMap& pmap)
        : _pmap(pmap) {}

    typedef typename property_traits<PropertyMap>::value_type value_type;

    typedef typename mpl::if_<
        typename return_reference::apply<value_type>::type,
        value_type&,
        value_type>::type reference;

    template <class PythonDescriptor>
    reference GetValue(const PythonDescriptor& key)
    {
        key.CheckValid();
        return get(_pmap, key.GetDescriptor());
    }

    template <class PythonDescriptor>
    void SetValue(const PythonDescriptor& key, const value_type& val)
    {
        set_value(key, val,
                  is_convertible<
                  typename property_traits<PropertyMap>::category,
                  writable_property_map_tag>());
    }

    template <class PythonDescriptor>
    void set_value(const PythonDescriptor& key, const value_type& val,
                   true_type)
    {
        key.CheckValid();
        put(_pmap, key.GetDescriptor(), val);
    }

    template <class PythonDescriptor>
    void set_value(const PythonDescriptor& key, const value_type& val,
                   false_type)
    {
        throw ValueException("property is read-only");
    }

    size_t GetHash() const
    {
        return hash<size_t>()(size_t(this));
    }

    std::string GetType() const
    {
        using python::detail::gcc_demangle;
        if (is_same<typename mpl::find<value_types,value_type>::type,
                    typename mpl::end<value_types>::type>::value)
            return gcc_demangle(typeid(value_type).name());
        else
            return type_names[mpl::find<value_types,
                                        value_type>::type::pos::value];
    }

    boost::any GetMap() const
    {
        return _pmap;
    }

    boost::any GetDynamicMap() const
    {
        return (dynamic_property_map*)
            (new boost::detail::dynamic_property_map_adaptor<PropertyMap>
             (_pmap));
    }

    python::object GetArray(size_t size)
    {
        typedef typename mpl::or_<
              typename mpl::or_<
                   is_same<PropertyMap,
                           GraphInterface::vertex_index_map_t>,
                   is_same<PropertyMap,
                           GraphInterface::edge_index_map_t> >::type,
              typename mpl::not_<
                  typename mpl::has_key<numpy_types, value_type>::type >
            ::type>::type isnt_vector_map;
        return get_array(_pmap, size, isnt_vector_map());
    }

    python::object get_array(PropertyMap pmap, size_t size, mpl::bool_<false>)
    {
        _pmap.reserve(size);
        return wrap_vector_not_owned(_pmap.get_storage());
    }

    python::object get_array(PropertyMap pmap, size_t size, mpl::bool_<true>)
    {
        return python::object();
    }

    bool IsWritable() const
    {
        return is_convertible<
            typename property_traits<PropertyMap>::category,
            readable_property_map_tag>::value;
    }

private:
    PropertyMap _pmap; // hold an internal copy, since it's cheap
};


//
// Create new properties
//

struct new_property_map
{
    template <class ValueType, class IndexMap>
    void operator()(ValueType, IndexMap index, const string& type_name,
                    python::object& new_prop, bool& found) const
    {
        size_t i = mpl::find<value_types,ValueType>::type::pos::value;
        if (type_name == type_names[i])
        {
            typedef typename property_map_type::apply<ValueType, IndexMap>::type
                map_t;
            map_t prop(index);
            new_prop = python::object(PythonPropertyMap<map_t>(prop));
            found = true;
        }
    }
};

template <class IndexMap>
python::object new_property(const string& type, IndexMap index_map)
{
    python::object prop;
    bool found = false;
    mpl::for_each<value_types>(lambda::bind<void>(new_property_map(),
                                                  lambda::_1, index_map,
                                                  lambda::var(type),
                                                  lambda::var(prop),
                                                  lambda::var(found)));
    if (!found)
        throw ValueException("Invalid property type: " + type);
    return prop;
}

} //graph_tool namespace

#endif
