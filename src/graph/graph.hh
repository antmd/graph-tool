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
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#ifndef GRAPH_HH
#define GRAPH_HH

#include <tr1/unordered_map>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/vector_property_map.hpp>
#include <boost/dynamic_property_map.hpp>
#include <boost/variant.hpp>
#include <boost/python/object.hpp>
#include "histogram.hh"
#include "config.h"
#include "graph_properties.hh"

namespace graph_tool
{

//==============================================================================
// GraphInterface
// this structure is an interface to the internally kept graph
//==============================================================================

class GraphInterface
{
public:
    GraphInterface(); 
    ~GraphInterface();

    enum degree_t
    {
        IN_DEGREE,
        OUT_DEGREE,
        TOTAL_DEGREE,
        SCALAR
    };

    // histogram types
    typedef std::tr1::unordered_map<double,size_t> hist_t;
    typedef std::tr1::unordered_map<std::pair<double,double>,double, boost::hash<std::pair<double,double> > > hist2d_t;
    typedef std::tr1::unordered_map<boost::tuple<double,double,double>,double_t, boost::hash<boost::tuple<double,double,double> > > hist3d_t;
    typedef std::tr1::unordered_map<double,std::pair<double,double> > avg_corr_t;

    typedef boost::variant<degree_t,std::string> deg_t;

    //graph generation
    typedef boost::function<double (size_t j, size_t k)> pjk_t;
    typedef boost::function<std::pair<size_t,size_t> (double r1, double r2)> inv_ceil_t;
    typedef boost::function<double (size_t jl, size_t kl, size_t j, size_t k)> corr_t;
    typedef boost::function<std::pair<size_t,size_t>  (double r1, double r2, size_t j, size_t k)> inv_corr_t;
    void GenerateCorrelatedConfigurationalModel(size_t N, pjk_t p, pjk_t ceil, inv_ceil_t inv_ceil, double ceil_pjk_bound, corr_t corr, corr_t ceil_corr, inv_corr_t inv_ceil_corr, double ceil_corr_bound, bool undirected_corr, size_t seed, bool verbose);

    // basic stats
    size_t GetNumberOfVertices() const;
    size_t GetNumberOfEdges() const;
    hist_t GetVertexHistogram(deg_t degree) const;
    hist_t GetEdgeHistogram(std::string property) const;

    //correlations
    hist2d_t   GetCombinedVertexHistogram(deg_t degree1, deg_t degree2) const;    
    avg_corr_t GetAverageCombinedVertexCorrelation(deg_t degree1, deg_t degree2) const;
    hist2d_t   GetVertexCorrelationHistogram(deg_t degree1, deg_t degree2, std::string weight) const;
    hist3d_t   GetEdgeVertexCorrelationHistogram(deg_t deg1, std::string scalar, deg_t deg2) const;
    avg_corr_t GetAverageNearestNeighboursCorrelation(deg_t origin_degree, deg_t neighbour_degree, std::string weight) const;

    // mixing
    std::pair<double,double> GetAssortativityCoefficient(deg_t deg) const;
    std::pair<double,double> GetScalarAssortativityCoefficient(deg_t deg) const;

    //clustering
    void SetLocalClusteringToProperty(std::string property);
    std::pair<double,double> GetGlobalClustering();
    void SetExtendedClusteringToProperty(std::string property_prefix, size_t max_depth);

    // other
    void   LabelComponents(std::string property);
    void   LabelParallelEdges(std::string property);
    hist_t GetDistanceHistogram(std::string weight) const;
    hist_t GetSampledDistanceHistogram(std::string weight, size_t samples, size_t seed) const;
    double GetReciprocity() const;
    void   GetMinimumSpanningTree(std::string weight, std::string property);
    void   GetLineGraph(std::string out_file, std::string format);
    void   GetBetweenness(std::string weight, std::string edge_betweenness, std::string vertex_betweenness);
    double GetCentralPointDominance(std::string vertex_betweenness);
    
    // community structure
    enum comm_corr_t
    {
        ERDOS_REYNI,
        UNCORRELATED,
        CORRELATED
    };

    void   GetCommunityStructure(double gamma, comm_corr_t corr, size_t n_iter, double Tmin, double Tmax, size_t Nseeds, size_t seed, bool verbose, std::string history_file, std::string weight, std::string property);
    double GetModularity(std::string weight, std::string property);
    void   GetCommunityNetwork(std::string property, std::string out_file, std::string format) const;

    // filtering
    void SetDirected(bool directed) {_directed = directed;}
    bool GetDirected() const {return _directed;}

    void SetReversed(bool reversed) {_reversed = reversed;}
    bool GetReversed() const {return _reversed;}

    void SetVertexFilterProperty(std::string property);
    std::string GetVertexFilterProperty() const {return _vertex_filter_property;}
    void SetVertexFilterRange(std::pair<double,double> allowed_range);
    std::pair<double, double> GetVertexFilterRange() const {return _vertex_range;}
    bool IsVertexFilterActive() const;

    void SetEdgeFilterProperty(std::string property);
    std::string GetEdgeFilterProperty() const {return _edge_filter_property;}
    void SetEdgeFilterRange(std::pair<double,double> allowed_range) {_edge_range = allowed_range;}
    std::pair<double,double> GetEdgeFilterRange() const {return _edge_range;}
    bool IsEdgeFilterActive() const;

    void SetGenericVertexFilter(boost::python::object filter);
    void SetGenericEdgeFilter(boost::python::object filter);

    // modification
    void RemoveVertexProperty(std::string property);
    void RemoveEdgeProperty(std::string property);
    void RemoveGraphProperty(std::string property);
    void InsertEdgeIndexProperty(std::string property);
    void InsertVertexIndexProperty(std::string property);
    void EditVertexProperty(std::string property, std::string type, boost::python::object op);
    void EditEdgeProperty(std::string property, std::string type, boost::python::object op);
    void EditGraphProperty(std::string property, std::string type, boost::python::object op);
    void ListProperties() const;

    // layout
    void ComputeGraphLayoutGursoy(size_t iter = 0, size_t seed = 4357);
    void ComputeGraphLayoutSpringBlock(size_t iter = 0, size_t seed = 4357);

    // i/o
    void WriteToFile(std::string s);
    void WriteToFile(std::string s, std::string format);
    void ReadFromFile(std::string s);
    void ReadFromFile(std::string s, std::string format);

    // signal handling
    void InitSignalHandling();

    // the following defines the edges' internal properties
    typedef boost::property<boost::edge_index_t, size_t> EdgeProperty;

    // this is the main graph type
    typedef boost::adjacency_list <boost::vecS, // edges
                                   boost::vecS, // vertices
                                   boost::bidirectionalS,
                                   boost::no_property,
                                   EdgeProperty >  multigraph_t;

private:
    template <class Action, class ReverseCheck, class DirectedCheck> 
    friend void check_filter(const GraphInterface &g, Action a, ReverseCheck, DirectedCheck);
    template <class Graph, class Action, class ReverseCheck, class DirectedCheck> 
    friend void check_python_filter(const Graph& g, const GraphInterface &gi, Action a, bool& found, ReverseCheck, DirectedCheck);

    friend class scalarS;

    multigraph_t _mg;

    bool _reversed;
    bool _directed;

    // vertex index map
    typedef boost::property_map<multigraph_t, boost::vertex_index_t>::type vertex_index_map_t;
    vertex_index_map_t _vertex_index;

    // edge index map
    typedef boost::property_map<multigraph_t, boost::edge_index_t>::type edge_index_map_t;
    edge_index_map_t _edge_index;

    // graph properties
    boost::dynamic_properties _properties;

    // vertex filter
    std::string _vertex_filter_property;
    typedef boost::variant<boost::vector_property_map<double, vertex_index_map_t>,
                           HashedDescriptorMap<vertex_index_map_t,double>,
                           boost::vector_property_map<size_t, vertex_index_map_t>,
                           HashedDescriptorMap<vertex_index_map_t, size_t>,
                           vertex_index_map_t,
                           DynamicPropertyMapWrap<double, boost::graph_traits<multigraph_t>::vertex_descriptor> > vertex_filter_map_t;
    vertex_filter_map_t _vertex_filter_map;
    std::pair<double,double> _vertex_range;
    boost::python::object _vertex_python_filter;

    // edge filter
    std::string _edge_filter_property;
    typedef boost::variant<boost::vector_property_map<double, edge_index_map_t>,
                           HashedDescriptorMap<edge_index_map_t, double>,                            
                           boost::vector_property_map<size_t, edge_index_map_t>,
                           HashedDescriptorMap<edge_index_map_t, size_t>,
                           edge_index_map_t,
                           DynamicPropertyMapWrap<double, boost::graph_traits<multigraph_t>::edge_descriptor> > edge_filter_map_t;
    edge_filter_map_t _edge_filter_map;
    std::pair<double,double> _edge_range;
    boost::python::object _edge_python_filter;
    
};

std::pair<GraphInterface::degree_t,std::string> get_degree_type(GraphInterface::deg_t degree);

//==============================================================================
// GraphException
//==============================================================================

class GraphException : public std::exception
{
    std::string _error;
public:
    GraphException(std::string error) {_error = error;}
    virtual ~GraphException() throw () {}
    virtual const char * what () const throw () {return _error.c_str();}
};

} //namespace graph_tool

//#include "node_graph_io.hh"

#endif



