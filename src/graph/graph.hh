
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

#include <boost/graph/filtered_graph.hpp>

#include <boost/vector_property_map.hpp>
#include <boost/dynamic_property_map.hpp>
#include <boost/variant.hpp>
#include <boost/python/object.hpp>
#include <boost/python/dict.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include "histogram.hh"
#include "config.h"
#include "graph_properties.hh"

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
    ~GraphInterface();

    // this enum specifies all the different types of degree
    enum degree_t
    {
        IN_DEGREE,
        OUT_DEGREE,
        TOTAL_DEGREE, // in + out
    };

    typedef variant<degree_t,string> deg_t; // useful when function also expects
                                            // a scalar vertex property

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
    bool IsVertexFilterActive() const;
    void SetEdgeFilterProperty(string property, bool invert);
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
    void RandomRewire(std::string strat, bool self_loops, bool parallel_edges,
                      size_t seed);

    // i/o
    void WriteToFile(string s);
    void WriteToFile(string s, string format);
    void ReadFromFile(string s);
    void ReadFromFile(string s, string format);

    //
    // Algorithms
    // Below are all the algorithms that operate somehow on the graph
    //

    // basic statistics
    hist_t GetVertexHistogram(deg_t degree) const;
    hist_t GetEdgeHistogram(string property) const;

    // correlations
    hist2d_t   GetCombinedVertexHistogram(deg_t degree1, deg_t degree2) const;
    avg_corr_t GetAverageCombinedVertexCorrelation(deg_t degree1,
                                                   deg_t degree2)const;
    hist2d_t   GetVertexCorrelationHistogram(deg_t degree1, deg_t degree2,
                                             string weight) const;
    hist3d_t   GetEdgeVertexCorrelationHistogram(deg_t deg1, string scalar,
                                                 deg_t deg2) const;
    avg_corr_t GetAverageNearestNeighboursCorrelation(deg_t origin_degree,
                                                      deg_t neighbour_degree,
                                                      string weight) const;

    // vertex mixing
    pair<double,double> GetAssortativityCoefficient(deg_t deg) const;
    pair<double,double> GetScalarAssortativityCoefficient(deg_t deg) const;

    // clustering
    void SetLocalClusteringToProperty(string property);
    pair<double,double> GetGlobalClustering();
    void SetExtendedClusteringToProperty(string property_prefix,
                                         size_t max_depth);

    // other
    void   LabelComponents(string property);
    void   LabelParallelEdges(string property);
    hist_t GetDistanceHistogram(string weight) const;
    hist_t GetSampledDistanceHistogram(string weight, size_t samples,
                                       size_t seed) const;
    double GetReciprocity() const;
    void   GetMinimumSpanningTree(string weight, string property);
    void   GetLineGraph(string out_file, string format);
    void   GetBetweenness(string weight, string edge_betweenness,
                          string vertex_betweenness);
    double GetCentralPointDominance(string vertex_betweenness);

    // community structure
    enum comm_corr_t // null model correlation type
    {
        ERDOS_REYNI,
        UNCORRELATED,
        CORRELATED
    };

    void   GetCommunityStructure(double gamma, comm_corr_t corr, size_t n_iter,
                                 double Tmin, double Tmax, size_t Nseeds,
                                 size_t seed, bool verbose, string history_file,
                                 string weight, string property);
    double GetModularity(string weight, string property);
    // TODO: this should return a GraphInterface type
    void   GetCommunityNetwork(string property, string size_property,
                               string out_file, string format) const;

    // Graph generation
    void GenerateCorrelatedConfigurationalModel
        (size_t N, python::object ppjk, python::object pceil_pjk,
         python::object pinv_ceil_pjk, double ceil_pjk_bound,
         python::object pcorr, python::object pceil_corr,
         python::object pinv_ceil_corr, double ceil_corr_bound,
         bool undirected_corr, size_t seed, bool verbose);

    // graph layout
    void ComputeGraphLayoutGursoy(string prop, string weight, string topology,
                                  size_t iter = 0, size_t seed = 4357);
    void ComputeGraphLayoutSpringBlock(string prop, string weight, string type,
                                       size_t iter = 0, 
                                       bool progressive = false,
                                       size_t seed = 4357);

    // python interface
    python::object Vertices() const;
    python::object Vertex(size_t i) const;
    python::object Edges() const;

    python::object AddVertex();
    void RemoveVertex(python::object v);
    python::object AddEdge(python::object s, python::object t);
    void RemoveEdge(python::object e);

    python::dict GetVertexProperties() const;
    python::dict GetEdgeProperties() const;
    python::dict GetGraphProperties() const;

    // used for graph properties
    graph_property_tag GetDescriptor() const { return graph_property_tag(); }

    void ExportPythonInterface() const;

    // signal handling
    void InitSignalHandling();

    // arbitrary code execution, for run-time code integration
    template <class Action, class Args>
    void RunAction(const Action &a, const Args& args)
    {
        using namespace boost::lambda;
        run_action(*this, bind<void>(a, boost::lambda::_1, _vertex_index,
                                     _edge_index, var(_properties), var(args)));
    }

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

    size_t GetEdgeHash(const edge_t& e) const;

private:
    // Gets the encapsulated graph view. See graph_filtering.cc for details
    boost::any GetGraphView() const;

    // Generic graph_action functor. See graph_filtering.hh for details.
    template <class Action,
              class TR1=boost::mpl::vector<>, class TR2=boost::mpl::vector<>,
              class TR3=boost::mpl::vector<>, class TR4=boost::mpl::vector<> >
    friend struct detail::graph_action;

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
    typedef vector_property_map<bool,vertex_index_map_t> vertex_filter_t;
    vertex_filter_t _vertex_filter_map;
    bool _vertex_filter_invert;
    bool _vertex_filter_active;

    // edge filter
    typedef vector_property_map<bool,edge_index_map_t> edge_filter_t;
    edge_filter_t _edge_filter_map;
    bool _edge_filter_invert;
    bool _edge_filter_active;
};
#pragma GCC visibility pop

pair<GraphInterface::degree_t,string>
get_degree_type(GraphInterface::deg_t degree);

// This is the main exception which will be thrown the outside world, when
// things go wrong

#pragma GCC visibility push(default)
class GraphException : public exception
{
    string _error;
public:
    GraphException(string error) {_error = error;}
    virtual ~GraphException() throw () {}
    virtual const char * what () const throw () {return _error.c_str();}
};
#pragma GCC visibility pop

} //namespace graph_tool

#endif



