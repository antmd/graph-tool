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

#include <boost/python.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/functional/hash.hpp>
#include <boost/iterator/iterator_facade.hpp>

#include "graph_selectors.hh"

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
template <class Graph, class Descriptor, class Iterator>
class PythonIterator
{
public:
    PythonIterator(const Graph& g, std::pair<Iterator,Iterator> e)
        : _g(g), _e(e) {}
    Descriptor Next()
    {
        if (_e.first == _e.second)
            python::objects::stop_iteration_error();
        Descriptor e(_g, *_e.first);
        ++_e.first;
        return e;
    }
private:
    const Graph& _g;
    std::pair<Iterator,Iterator> _e;
};


// forward declaration of PythonEdge
template <class Graph>
class PythonEdge;

// below are classes related to the PythonVertex type
template <class Graph>
    class PythonVertex
{
public:
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

    PythonVertex(const Graph& g, vertex_descriptor v):
        _g(g), _v(v) {}

    PythonVertex(const PythonVertex& v): _g(v._g)
    {
        _v = v._v;
    }

    const vertex_descriptor& GetDescriptor() const
    {
        return _v;
    }

    size_t GetInDegree() const
    {
        return in_degreeS()(_v, _g);
    }

    size_t GetOutDegree() const
    {
        return out_degreeS()(_v, _g);
    }

    // provide iterator support for out_edges

    typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iterator;
    PythonIterator<Graph,PythonEdge<Graph>,out_edge_iterator>
    OutEdges() const
    {
        return PythonIterator<Graph,PythonEdge<Graph>,out_edge_iterator>
            (_g, out_edges(_v, _g));
    }

    typedef typename in_edge_iteratorS<Graph>::type in_edge_iterator;
    PythonIterator<Graph,PythonEdge<Graph>, in_edge_iterator>
    InEdges() const
    {
        return PythonIterator<Graph,PythonEdge<Graph>,in_edge_iterator>
            (_g, in_edge_iteratorS<Graph>::in_edges(_v, _g));
    }

    std::string GetString() const
    {
        return lexical_cast<std::string>(_v);
    }

    size_t GetHash() const
    {
        return hash<vertex_descriptor>()(_v);
    }

    bool operator==(const PythonVertex& other) const
    {
        return other._v == _v;
    }

    bool operator!=(const PythonVertex& other) const
    {
        return other._v != _v;
    }

private:
    const Graph& _g;
    vertex_descriptor _v;
};

// below are classes related to the PythonEdge type

template <class Graph>
class PythonEdge
{
public:
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

    PythonEdge(const Graph& g, edge_descriptor e): _g(g), _e(e) {}
    PythonEdge(const PythonEdge& e): _g(e._g)
    {
        _e = e._e;
    }

    const edge_descriptor& GetDescriptor() const
    {
        return _e;
    }

    PythonVertex<Graph> GetSource() const
    {
        return PythonVertex<Graph>(_g, source(_e, _g));
    }

    PythonVertex<Graph> GetTarget() const
    {
        return PythonVertex<Graph>(_g, target(_e, _g));
    }

    std::string GetString() const
    {
        return "("+GetSource().GetString() + "," +
            GetTarget().GetString() + ")";
    }

    size_t GetHash() const
    {
        vertex_descriptor s,t;
        s = source(_e, _g);
        t = target(_e, _g);
        return hash<std::pair<vertex_descriptor, vertex_descriptor> >()
            (make_pair(s,t));
    }

    bool operator==(const PythonEdge& other) const
    {
        return other._e == _e;
    }

    bool operator!=(const PythonEdge& other) const
    {
        return other._e != _e;
    }

private:
    const Graph& _g;
    edge_descriptor _e;
};

template <class ValueType>
class PythonPropertyMap
{
public:
    PythonPropertyMap(const std::string& name, dynamic_property_map& pmap)
        : _name(name), _pmap(pmap) {}

    template <class PythonDescriptor>
    ValueType GetValue(const PythonDescriptor& key) const
    {
        any val = _pmap.get(key.GetDescriptor());
        return any_cast<ValueType>(val);
    }

    template <class PythonDescriptor>
    void SetValue(const PythonDescriptor& key, const ValueType& val)
    {
        _pmap.put(key.GetDescriptor(), val);
    }

    size_t GetHash() const
    {
        return hash<std::string>()(_name + this->GetType());
    }

    std::string GetType() const
    {
        return type_names[mpl::find<value_types,ValueType>::type::pos::value];
    }

private:
    const std::string& _name;
    dynamic_property_map& _pmap;
};


} //graph_tool namespace

#endif
