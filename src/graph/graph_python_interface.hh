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
#include <boost/lambda/bind.hpp>
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
template <class Descriptor, class Iterator>
class PythonIterator
{
public:
    PythonIterator(const GraphInterface& gi, std::pair<Iterator,Iterator> e)
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
    const GraphInterface& _gi;
    std::pair<Iterator,Iterator> _e;
};


// forward declaration of PythonEdge
template <class Graph>
class PythonEdge;

// below are classes related to the PythonVertex type
class PythonVertex
{
public:
    PythonVertex(const GraphInterface& gi, GraphInterface::vertex_t v):
        _gi(gi), _v(v) {}

    PythonVertex(const PythonVertex& v): _gi(v._gi)
    {

        _v = v._v;
    }

    GraphInterface::vertex_t GetDescriptor() const
    {
        return _v;
    }

    template <class DegSelector>
    struct get_degree
    {
        template<class Graph>
        void operator()(const Graph* gp,
                        typename graph_traits<Graph>::vertex_descriptor v,
                        size_t& deg) const
        {
            deg = DegSelector()(v, *gp);
        }
    };
    

    size_t GetInDegree() const
    {
        size_t in_deg;
        run_action<>()(_gi,lambda::bind<void>(get_degree<in_degreeS>(),
                                              lambda::_1, _v,
                                              lambda::var(in_deg)))();
        return in_deg;
    }

    size_t GetOutDegree() const
    {
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
        void operator()(const Graph* gp, const GraphInterface& gi,
                        typename graph_traits<Graph>::vertex_descriptor v,
                        python::object& iter) const
        {
            typedef typename graph_traits<Graph>::out_edge_iterator
                out_edge_iterator;
            iter = python::object(PythonIterator<PythonEdge<Graph>,
                                                 out_edge_iterator>
                                  (gi, out_edges(v, *gp)));
        }
    };

    python::object
    OutEdges() const
    {
        python::object iter;
        run_action<>()(_gi, lambda::bind<void>(get_out_edges(), lambda::_1,
                                               lambda::var(_gi), _v,
                                               lambda::var(iter)))();
        return iter;
    }

    struct get_in_edges
    {
        template<class Graph>
        void operator()(const Graph* gp, const GraphInterface& gi,
                        typename graph_traits<Graph>::vertex_descriptor v,
                        python::object& iter) const
        {
            typedef typename in_edge_iteratorS<Graph>::type
                in_edge_iterator;
            iter = python::object
                (PythonIterator<PythonEdge<Graph>,in_edge_iterator>
                 (gi, in_edge_iteratorS<Graph>::get_edges(v, *gp)));
        }
    };

    python::object
    InEdges() const
    {
        python::object iter;
        run_action<>()(_gi, lambda::bind<void>(get_in_edges(), lambda::_1, 
                                               lambda::var(_gi), _v,
                                               lambda::var(iter)))();
        return iter;
    }
    
    std::string GetString() const
    {
        return lexical_cast<std::string>(_v);
    }

    size_t GetHash() const
    {
        return hash<GraphInterface::vertex_t>()(_v);
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
    const GraphInterface& _gi;
    GraphInterface::vertex_t _v;
};

// below are classes related to the PythonEdge type

template <class Graph>
class PythonEdge
{
public:
    typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
    PythonEdge(const GraphInterface& gi, edge_descriptor e)
        : _gi(gi), _e(e) {}
    PythonEdge(const PythonEdge& e): _gi(e._gi)
    {
        _e = e._e;
    }

    const GraphInterface::edge_t& GetDescriptor() const
    {
        return _e;
    }

    struct get_source
    {
        template<class GraphType>
        void operator()(const GraphType* gp, const GraphInterface& gi,
                        const edge_descriptor& edge, python::object& vertex)
            const
        {
            vertex = python::object(PythonVertex(gi, source(edge, *gp)));
        }
    };
    
    python::object GetSource() const
    {
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
        void operator()(const GraphType* gp, const GraphInterface& gi,
                        const edge_descriptor& edge, 
                        python::object& vertex) const
        {
            vertex = python::object(PythonVertex(gi, target(edge, *gp)));
        }
    };
    
    python::object GetTarget() const
    {
        python::object v;
        run_action<>()(_gi, lambda::bind<void>(get_target(), lambda::_1, 
                                               lambda::var(_gi),
                                               lambda::var(_e),
                                               lambda::var(v)))();
        return v;        
    }

    std::string GetString() const
    {
        PythonVertex& src = python::extract<PythonVertex&>(GetSource());
        PythonVertex& tgt = python::extract<PythonVertex&>(GetTarget());
        return "(" + src.GetString() + "," + tgt.GetString() + ")";
    }

    size_t GetHash() const
    {
        return _gi.GetEdgeHash(_e);;
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
    const GraphInterface& _gi;
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
