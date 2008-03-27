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

#ifndef GRAPH_HH
#define GRAPH_HH

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <boost/vector_property_map.hpp>
#include <boost/dynamic_property_map.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/variant.hpp>
#include <boost/python/object.hpp>
#include <boost/python/dict.hpp>
#include <boost/mpl/vector.hpp>
#include "config.h"

namespace graph_tool
{
using namespace std;
using namespace boost;

// GraphInterface
// this class is the main interface to the internally kept graph. This is how
// the external world will manipulate the graph. All the algorithms should be
// registered here. This class will be exported to python in graph_bind.hh

namespace detail
{
// Generic graph_action functor. See graph_filtering.hh for details.
template <class Action, class GraphViews,
          class TR1=boost::mpl::vector<>, class TR2=boost::mpl::vector<>,
          class TR3=boost::mpl::vector<>, class TR4=boost::mpl::vector<> >
struct graph_action;
}

// default visibility is necessary for related typeinfo objects to work across
// DSO boundaries
#pragma GCC visibility push(default)
class GraphInterface
{
public:
    GraphInterface();
    GraphInterface(const GraphInterface& gi);
    ~GraphInterface();

    // useful enums

    typedef enum
    {
        IN_DEGREE,
        OUT_DEGREE,
        TOTAL_DEGREE
    } degree_t;

    typedef boost::variant<degree_t, string> deg_t; // general "degree" type,
                                                    // i.e., either a degree_t
                                                    // above or a string
                                                    // representing a scalar
                                                    // vertex property

    //
    // Basic manipulation
    //

    size_t GetNumberOfVertices() const;
    size_t GetNumberOfEdges() const;
    void SetDirected(bool directed) {_directed = directed;}
    bool GetDirected() const {return _directed;}
    void SetReversed(bool reversed) {_reversed = reversed;}
    bool GetReversed() const {return _reversed;}

    // graph filtering
    void SetVertexFilterProperty(string property, bool invert);
    pair<string, bool> GetVertexFilterProperty() const;
    bool IsVertexFilterActive() const;
    void SetEdgeFilterProperty(string property, bool invert);
    pair<string, bool> GetEdgeFilterProperty() const;
    bool IsEdgeFilterActive() const;

    // graph modification
    void RemoveVertexProperty(string property);
    void RemoveEdgeProperty(string property);
    void RemoveGraphProperty(string property);
    void InsertEdgeIndexProperty(string property);
    void InsertVertexIndexProperty(string property);
    void AddVertexProperty(string property, string type);
    void AddEdgeProperty(string property, string type);
    void AddGraphProperty(string property, string type);
    void ReIndexEdges();
    void PurgeVertices(); // removes filtered vertices
    void PurgeEdges();    // removes filtered edges
    void Clear();

    // i/o
    void WriteToFile(string s, python::object pf, string format);
    void ReadFromFile(string s, python::object pf, string format);

    //
    // python interface
    //
    python::object Vertices() const;
    python::object Vertex(size_t i) const;
    python::object Edges() const;

    python::object AddVertex();
    void           RemoveVertex(const python::object& v);
    python::object AddEdge(const python::object& s, const python::object& t);
    void           RemoveEdge(const python::object& e);

    python::dict GetVertexProperties() const;
    python::dict GetEdgeProperties() const;
    python::dict GetGraphProperties() const;

    // used for graph properties
    graph_property_tag GetDescriptor() const { return graph_property_tag(); }
    bool CheckValid() const {return true;}

    void ExportPythonInterface() const;

    // signal handling
    void InitSignalHandling();

    //
    // Internal types
    //

    // the following defines the edges' internal properties
    typedef property<edge_index_t, size_t> EdgeProperty;

    // this is the main graph type
    typedef adjacency_list <vecS, // edges
                            vecS, // vertices
                            bidirectionalS,
                            no_property,
                            EdgeProperty>  multigraph_t;
    typedef graph_traits<multigraph_t>::vertex_descriptor vertex_t;
    typedef graph_traits<multigraph_t>::edge_descriptor edge_t;

    typedef property_map<multigraph_t,vertex_index_t>::type vertex_index_map_t;
    typedef property_map<multigraph_t,edge_index_t>::type edge_index_map_t;

private:
    // Gets the encapsulated graph view. See graph_filtering.cc for details
    boost::any GetGraphView() const;

    // Generic graph_action functor. See graph_filtering.hh for details.
    template <class Action,
              class TR1=boost::mpl::vector<>, class TR2=boost::mpl::vector<>,
              class TR3=boost::mpl::vector<>, class TR4=boost::mpl::vector<> >
    friend struct detail::graph_action;

    // Arbitrary code execution
    template <class Action>
    friend void RunAction(GraphInterface &g, const Action& a);

    friend boost::any degree_selector(deg_t deg, const GraphInterface&gi);

    friend boost::any vertex_prop(const string& name, const GraphInterface& gi);
    friend boost::any edge_prop(const string& name, const GraphInterface& gi);
    friend boost::any graph_prop(const string& name, const GraphInterface& gi);

    // python interface
    friend class PythonVertex;
    template <class Graph>
    friend class PythonEdge;

    // this is the main graph
    multigraph_t _mg;

    // this will hold an instance of the graph views at run time
    vector<boost::any> _graph_views;

    // reverse and directed states
    bool _reversed;
    bool _directed;

    // vertex index map
    vertex_index_map_t _vertex_index;

    // edge index map
    edge_index_map_t _edge_index;

    // graph properties
    dynamic_properties _properties;

    // vertex filter
    typedef vector_property_map<uint8_t,vertex_index_map_t> vertex_filter_t;
    vertex_filter_t _vertex_filter_map;
    string _vertex_filter_property;
    bool _vertex_filter_invert;
    bool _vertex_filter_active;

    // edge filter
    typedef vector_property_map<uint8_t,edge_index_map_t> edge_filter_t;
    edge_filter_t _edge_filter_map;
    string _edge_filter_property;
    bool _edge_filter_invert;
    bool _edge_filter_active;
};
#pragma GCC visibility pop

// This is the main exception which will be thrown the outside world, when
// things go wrong

#pragma GCC visibility push(default)
class GraphException : public exception
{
    string _error;
public:
    GraphException(const string& error) {_error = error;}
    virtual ~GraphException() throw () {}
    virtual const char * what () const throw () {return _error.c_str();}
protected:
    virtual void SetError(const string& error) {_error = error;}
};
#pragma GCC visibility pop

} //namespace graph_tool

#endif



