// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006  Tiago de Paula Peixoto <tiago@forked.de>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
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

#include <algorithm>
#include <tr1/unordered_set>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/gursoy_atun_layout.hpp>
#include <boost/graph/fruchterman_reingold.hpp>
#include <boost/graph/random_layout.hpp>
#include <boost/random.hpp>
#include <boost/python/make_function.hpp>

#include <unistd.h>    /* standard unix functions, like getpid()         */
#include <sys/types.h> /* various type definitions, like pid_t           */
#include <signal.h>    /* signal name macros, and the signal() prototype */

#include "graph.hh"
#include "histogram.hh"
#include "graph_filtering.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;


pair<GraphInterface::degree_t,string> graph_tool::get_degree_type(GraphInterface::deg_t degree)
{
    GraphInterface::degree_t deg;
    string name;
    try 
    {
	deg = boost::get<GraphInterface::degree_t>(degree);
    }
    catch (bad_get)
    {
	name = boost::get<std::string>(degree);
	deg = GraphInterface::SCALAR;
    }
    return make_pair(deg,name);
}


//==============================================================================
// GraphInterface()
//==============================================================================
GraphInterface::GraphInterface()
:_mg(), 
 _reversed(false),
 _directed(true),
 _vertex_index(get(vertex_index,_mg)),
 _edge_index(get(edge_index_t(),_mg)),
 _vertex_range(make_pair(numeric_limits<double>::min(), numeric_limits<double>::max())),
 _edge_range(make_pair(numeric_limits<double>::min(), numeric_limits<double>::max()))
{

}

//==============================================================================
// ~GraphInterface()
//==============================================================================
GraphInterface::~GraphInterface()
{
}

//==============================================================================
// SetVertexFilter()
//==============================================================================

python::object python_range_filter(python::object properties, string filter_property, pair<double,double> range) 
{
    bool accept = true;
    if (python::extract<double>(properties[filter_property]) < range.first || python::extract<double>(properties[filter_property]) > range.second)
	accept = false;
    return python::object(accept);
}

void GraphInterface::SetVertexFilterProperty(string property)
{
    _vertex_filter_property = property;
    
    if (property != "")
    {
	try
	{
	    _vertex_filter_map = get_static_property_map<vertex_filter_map_t>(find_property_map(_properties, property, 
												typeid(graph_traits<multigraph_t>::vertex_descriptor)));
	}
	catch (bad_cast)
	{
	    // set generic vertex filter instead
	    function<python::object(python::object)> filter = bind<python::object>(python_range_filter, _1, property, _vertex_range);
	    SetGenericVertexFilter(python::make_function(filter, python::default_call_policies(), mpl::vector<python::object,python::object>::type()));
	}
	catch (property_not_found) 
	{
	    throw GraphException("property " + property + " not found");
	}
    }
}

bool GraphInterface::IsVertexFilterActive() const { return _vertex_filter_property != "" || _vertex_python_filter != python::object(); }

void GraphInterface::SetEdgeFilterProperty(string property) 
{
    _edge_filter_property = property;
    
    if (property != "")
    {
	try
	{
	    _edge_filter_map = get_static_property_map<edge_filter_map_t>(find_property_map(_properties, property, 
											    typeid(graph_traits<multigraph_t>::edge_descriptor)));
	}
	catch (bad_cast)
	{
	    // set generic edge filter instead
	    function<python::object(python::object)> filter = bind<python::object>(python_range_filter, _1, property, _edge_range);
	    SetGenericEdgeFilter(python::make_function(filter, python::default_call_policies(), mpl::vector<python::object,python::object>::type()));
	}
	catch (property_not_found) 
	{
	    throw GraphException("property " + property + " not found");
	}
    }
}

bool GraphInterface::IsEdgeFilterActive() const {return _edge_filter_property != "" || _edge_python_filter != python::object();}


//==============================================================================
// GetNumberOfVertices()
//==============================================================================
size_t GraphInterface::GetNumberOfVertices() const
{
    size_t n = 0;
    if (IsVertexFilterActive())
	check_filter(*this,var(n)=bind<size_t>(HardNumVertices(),_1),reverse_check(),directed_check()); 
    else
        check_filter(*this,var(n)=bind<size_t>(SoftNumVertices(),_1),reverse_check(),directed_check());
    return n;
}

//==============================================================================
// GetNumberOfEdges()
//==============================================================================
size_t GraphInterface::GetNumberOfEdges() const
{
    size_t n = 0;
    if (IsEdgeFilterActive() || IsVertexFilterActive())
	check_filter(*this,var(n)=bind<size_t>(HardNumEdges(),_1),reverse_check(),directed_check()); 
    else
        check_filter(*this,var(n)=bind<size_t>(SoftNumEdges(),_1),reverse_check(),directed_check());
    return n;
}

//==============================================================================
// get_vertex_histogram
// generalized functor to obtain histogram of different types of "degrees"
//==============================================================================
template <class DegreeSelector>
struct get_vertex_histogram
{
    get_vertex_histogram(DegreeSelector& deg): _degree(deg) {}
    template <class Graph, class Hist>
    void operator()(const Graph& g, Hist& hist) const
    {
	typename graph_traits<Graph>::vertex_iterator v, v_begin, v_end;
        tie(v_begin, v_end) = vertices(g);
        for(v = v_begin; v != v_end; ++v)
            hist[_degree(*v,g)]++;
    }
    DegreeSelector& _degree;
};

struct choose_vertex_histogram
{
    choose_vertex_histogram(const GraphInterface& g, GraphInterface::deg_t deg, GraphInterface::hist_t& hist)
	: _g(g), _hist(hist) 
    {
	tie(_deg, _deg_name) = get_degree_type(deg);
    }
    template <class DegreeSelector>
    void operator()(DegreeSelector)
    {
	if (mpl::at<degree_selector_index,DegreeSelector>::type::value == _deg)
	{
	    DegreeSelector selector(_deg_name, _g);
	    check_filter(_g, bind<void>(get_vertex_histogram<DegreeSelector>(selector), _1, var(_hist)),reverse_check(),directed_check());
	}
    }
    const GraphInterface &_g;
    GraphInterface::hist_t &_hist;
    GraphInterface::degree_t _deg;
    string _deg_name;    
};

//==============================================================================
// GetVertexHistogram(deg_type)
//==============================================================================
GraphInterface::hist_t GraphInterface::GetVertexHistogram(GraphInterface::deg_t deg) const
{
    hist_t hist;
    try 
    {
	mpl::for_each<mpl::vector<in_degreeS, out_degreeS, total_degreeS, scalarS> >(choose_vertex_histogram(*this, deg, hist));
    }
    catch (dynamic_get_failure &e)
    {
	throw GraphException("error getting scalar property: " + string(e.what()));
    }

    return hist;
}

//==============================================================================
// get_edge_histogram
// generalized functor to obtain histogram of edge properties
//==============================================================================
struct get_edge_histogram
{
    get_edge_histogram(scalarS& prop): _prop(prop) {}
    template <class Graph, class Hist>
    void operator()(const Graph& g, Hist& hist) const
    {
	typename graph_traits<Graph>::edge_iterator e, e_begin, e_end;
        tie(e_begin, e_end) = edges(g);
        for(e = e_begin; e != e_end; ++e)
            hist[_prop(*e,g)]++;
    }
    scalarS& _prop;
};

//==============================================================================
// GetEdgeHistogram(property)
//==============================================================================
GraphInterface::hist_t GraphInterface::GetEdgeHistogram(string property) const
{
    hist_t hist;
    try 
    {
	scalarS prop(property, *this);
	check_filter(*this, bind<void>(get_edge_histogram(prop), _1, var(hist)),reverse_check(),directed_check());
    }
    catch (dynamic_get_failure &e)
    {
	throw GraphException("error getting scalar property: " + string(e.what()));
    }


    return hist;
}

//==============================================================================
// GetComponentSizeHistogram()
//==============================================================================

struct strong_component
{
    template <class Graph, class ComponentMap>
    void operator()(const Graph &g, ComponentMap comp_map) const
    {
	strong_components(g, comp_map);
    }
};

struct connected_component
{
    template <class Graph, class ComponentMap>
    void operator()(const Graph &g, ComponentMap comp_map) const
    {
	connected_components(g, comp_map);
    }
};

template <class ComponentSelector>
struct get_component_size_histogram
{
    template <class Graph, class IndexMap, class Hist>
    void operator()(const Graph &g, IndexMap index_map, Hist &hist) const
    {   
	typedef HashedDescriptorMap<IndexMap, size_t> comp_map_t;
	comp_map_t comp_map(index_map);

	_components(g, comp_map);

	tr1::unordered_map<size_t, size_t> comp_count;
	typename graph_traits<Graph>::vertex_iterator v, v_begin, v_end;
	tie(v_begin, v_end) = vertices(g);
	for(v = v_begin; v != v_end; ++v)
	    comp_count[comp_map[*v]]++;
	for(typeof(comp_count.begin()) iter = comp_count.begin(); iter != comp_count.end(); ++iter)
	    hist[iter->second]++;
    }
    ComponentSelector _components;
};

GraphInterface::hist_t GraphInterface::GetComponentSizeHistogram() const
{

    hist_t hist;
    if (_directed)
	check_filter(*this, bind<void>(get_component_size_histogram<strong_component>(), _1, _vertex_index, var(hist)),
		     reverse_check(),always_directed()); 
    else
	check_filter(*this, bind<void>(get_component_size_histogram<connected_component>(), _1, _vertex_index, var(hist)),
		     reverse_check(),always_undirected());
    return hist;
}

//==============================================================================
// RemoveVertexProperty(property)
//==============================================================================
void GraphInterface::RemoveVertexProperty(string property)
{
    dynamic_properties_copy dp;
    try
    {
	dynamic_property_map& prop_map = find_property_map(_properties, property, typeid(graph_traits<multigraph_t>::vertex_descriptor));
	for (typeof(_properties.begin()) iter = _properties.begin(); iter != _properties.end(); ++iter)
	{
	    if (iter->second != &prop_map)
		dp.insert(iter->first, auto_ptr<dynamic_property_map>(iter->second));
	}
    }
    catch (property_not_found)
    {
	throw GraphException("property '" + property + "' not found");
    }
    _properties = dp;
}

//==============================================================================
// RemoveEdgeProperty(property)
//==============================================================================
void GraphInterface::RemoveEdgeProperty(string property)
{
    dynamic_properties_copy dp;
    try
    {
	dynamic_property_map& prop_map = find_property_map(_properties, property, typeid(graph_traits<multigraph_t>::edge_descriptor));
	for (typeof(_properties.begin()) iter = _properties.begin(); iter != _properties.end(); ++iter)
	{
	    if (iter->second != &prop_map)
		dp.insert(iter->first, auto_ptr<dynamic_property_map>(iter->second));
	}
    }
    catch (property_not_found)
    {
	throw GraphException("property '" + property + "' not found");
    }
    _properties = dp;
}

//==============================================================================
// InsertEdgeIndexProperty(property)
//==============================================================================
void GraphInterface::InsertEdgeIndexProperty(string property)
{
    _properties.property(property, _edge_index);
}

//==============================================================================
// InsertVertexIndexProperty(property)
//==============================================================================
void GraphInterface::InsertVertexIndexProperty(string property)
{
    _properties.property(property, _vertex_index);
}

//==============================================================================
// RemoveParallelEdges(property)
//==============================================================================
struct remove_parallel_edges
{
    template <class Graph>
    void operator()(Graph &g) const
    {
	typename graph_traits<Graph>::edge_iterator e, e_end;
	for (tie(e, e_end) = edges(g); e != e_end; ++e)
	{
	    bool finished = false;
	    while (!finished)
	    {
		finished = true;
		typename graph_traits<Graph>::out_edge_iterator oe, oe_end;
		typename graph_traits<Graph>::vertex_descriptor s,t;
		s = source(*e,g);
		t = target(*e,g);
		for (tie(oe,oe_end) = out_edges(s,g); oe != oe_end; ++oe)
		    if (*oe != *e && target(*oe,g) == t) // is a parallel edge
		    {
			remove_edge(*oe, g);
			finished = false;
			break;
		    }
	    }
	}
    }
};

void GraphInterface::RemoveParallelEdges()
{
    if (_directed)
    {
	remove_parallel_edges()(_mg);
    }
    else
    {
	UndirectedAdaptor<multigraph_t> ug(_mg);
	remove_parallel_edges()(ug);
    }
}

//==============================================================================
// ComputeGraphLayoutGursoy(iter, seed)
//==============================================================================

struct compute_gursoy
{
    template <class Graph, class PosMap, class IndexMap>
    void operator()(Graph &g, size_t iter, size_t seed, PosMap pos, IndexMap index_map) const
    {
	mt19937 rng(seed);
	size_t n = HardNumVertices()(g);
	
	vector_property_map<square_topology<mt19937>::point_type, IndexMap> position_map(index_map);
	if (iter == 0)
	    iter = n;
	square_topology<mt19937> topology(rng, n);
	gursoy_atun_layout(g, topology, position_map, iterations(iter).
			   diameter_range(make_pair(sqrt(double(n)), 1.0)).
			   learning_constant_range(make_pair(0.8, 0.2)).
			   vertex_index_map(index_map));
	typename graph_traits<Graph>::vertex_iterator v, v_begin, v_end;
	tie(v_begin, v_end) = vertices(g);
	for(v = v_begin; v != v_end; ++v)
	{
	    pos[*v].x = position_map[*v][0];
	    pos[*v].y = position_map[*v][1];
	}
    }
};

typedef struct { double x; double y; } pos_t;
ostream& operator<<(ostream &o, const pos_t &p ) { o << p.x << "," << p.y; return o;}
istream& operator>>(istream &o, pos_t &p ) { char c; o >> p.x >> c >> p.y; return o;}

void GraphInterface::ComputeGraphLayoutGursoy(size_t iter, size_t seed)
{       
    // vertex postion map
    typedef HashedDescriptorMap<vertex_index_map_t, pos_t> pos_map_t;
    pos_map_t pos_map(_vertex_index);    

    check_filter(*this,bind<void>(compute_gursoy(),_1,iter,seed,var(pos_map),var(_vertex_index)),reverse_check(),directed_check());

    _properties.property("pos", pos_map);
}

//==============================================================================
// ComputeGraphLayoutSpringBlock(iter,seed)
//==============================================================================
struct compute_spring_block
{

    template <class Graph, class PosMap, class IndexMap>
    void operator()(Graph &g, size_t iter, size_t seed, PosMap pos, IndexMap index_map) const
    {
  
	mt19937 rng(seed);
	size_t n = HardNumVertices()(g);

        if (iter == 0)
            iter = 100;
        double radius = n*iter;
        random_graph_layout(g, pos, -radius/2, radius/2, -radius/2, radius/2, rng);
        fruchterman_reingold_force_directed_layout(g, pos, radius, radius, 
                                                   cooling(linear_cooling<double>(iter)).
                                                   vertex_index_map(index_map));
    }
    
};

void GraphInterface::ComputeGraphLayoutSpringBlock(size_t iter, size_t seed)
{
    // vertex postion map
    typedef HashedDescriptorMap<vertex_index_map_t,pos_t> pos_map_t;
    pos_map_t pos_map(_vertex_index);    

    check_filter(*this,bind<void>(compute_spring_block(),_1,iter,seed,var(pos_map),var(_vertex_index)),reverse_check(),directed_check()); 

    _properties.property("pos", pos_map);
}

//==============================================================================
// InitSignalHandling()
//==============================================================================

void catch_sig_stop(int sig_num)
{
    std::cerr << "graph-tool: received signal ";
    switch (sig_num)
    {
    case SIGINT:
        std::cerr << "SIGINT (Interrupt).";
        break;
    case SIGTERM:
        std::cerr << "SIGTERM (Termination).";
        break;
    case SIGQUIT:
        std::cerr << "SIGQUIT (Terminal quit).";
        break;
    case SIGHUP:
        std::cerr << "SIGHUP (Hangup).";
        break;
    case SIGPWR:
        std::cerr << "SIGPWR (Power failure restart).";
        break;
    case SIGSEGV:
        std::cerr << "SIGSEGV (Segmentation fault). There's a bug somewhere in the program. Go fix it.";
        break;
    case SIGBUS:
        std::cerr << "SIGBUS (BUS error). The bus is broken, I guess...";
        break;
    case SIGFPE:
        std::cerr << "SIGFPE (Floating-point exception). INFs and NANs are wreaking havoc.";
        break;
    case SIGILL:
        std::cerr << "SIGILL (Illegal instruction). Did you compile it right?";
        break;
    case SIGXCPU:
        std::cerr << "SIGXCPU (CPU limit exceeded). Time's over.";
        break;
    case SIGXFSZ:
        std::cerr << "SIGXFSZ (File size limit exceeded). The fascist sysadmin doesn't let you write more.";
        break;
    }
    std::cerr << " Bailing out cowardly and calling abort()." << std::endl;
    abort();
}


void GraphInterface::InitSignalHandling()
{
    signal(SIGINT, catch_sig_stop);
    signal(SIGTERM, catch_sig_stop);
    signal(SIGQUIT, catch_sig_stop);
    signal(SIGHUP, catch_sig_stop);
    signal(SIGPWR, catch_sig_stop);
    signal(SIGSEGV, catch_sig_stop);
    signal(SIGBUS, catch_sig_stop);
    signal(SIGFPE, catch_sig_stop);
    signal(SIGILL, catch_sig_stop);
    signal(SIGXCPU, catch_sig_stop);
    signal(SIGXFSZ, catch_sig_stop);
}
