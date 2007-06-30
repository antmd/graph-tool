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
#include "shared_map.hh"

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

void GraphInterface::SetVertexFilterProperty(string property)
{
#ifndef NO_RANGE_FILTERING
    _vertex_filter_property = property;
    
    if (property != "")
    {
        try 
        {
            dynamic_property_map& pmap = find_property_map(_properties, property, typeid(graph_traits<multigraph_t>::vertex_descriptor));

            if (get_static_property_map<vector_property_map<double,vertex_index_map_t> >(&pmap))
                _vertex_filter_map = get_static_property_map<vector_property_map<double,vertex_index_map_t> >(pmap);
            else if (get_static_property_map<HashedDescriptorMap<vertex_index_map_t,double> >(&pmap))
                _vertex_filter_map = get_static_property_map<HashedDescriptorMap<vertex_index_map_t,double> >(pmap);
            else if (get_static_property_map<vector_property_map<size_t,vertex_index_map_t> >(&pmap))
                _vertex_filter_map = get_static_property_map<vector_property_map<size_t,vertex_index_map_t> >(pmap);
            else if (get_static_property_map<HashedDescriptorMap<vertex_index_map_t,size_t> >(&pmap))
                _vertex_filter_map = get_static_property_map<HashedDescriptorMap<vertex_index_map_t,size_t> >(pmap);
            else if (get_static_property_map<vertex_index_map_t>(&pmap))
                _vertex_filter_map = get_static_property_map<vertex_index_map_t>(pmap);
            else 
                _vertex_filter_map = DynamicPropertyMapWrap<double, graph_traits<multigraph_t>::vertex_descriptor>(pmap);
        }
        catch (property_not_found) 
        {
            throw GraphException("property " + property + " not found");
        }
    }
#else
    if (property != "")
        throw GraphException("support for graph range filtering was not enabled during compilation.");
#endif
}

bool GraphInterface::IsVertexFilterActive() const { return _vertex_filter_property != "" || _vertex_python_filter != python::object(); }

void GraphInterface::SetEdgeFilterProperty(string property) 
{
#ifndef NO_RANGE_FILTERING
    _edge_filter_property = property;
    
    if (property != "")
    {
        try
        {
            dynamic_property_map& pmap = find_property_map(_properties, property, typeid(graph_traits<multigraph_t>::edge_descriptor));
            
            if (get_static_property_map<vector_property_map<double,edge_index_map_t> >(&pmap))
                _edge_filter_map = get_static_property_map<vector_property_map<double,edge_index_map_t> >(pmap);
            else if (get_static_property_map<HashedDescriptorMap<edge_index_map_t,double> >(&pmap))
                _edge_filter_map = get_static_property_map<HashedDescriptorMap<edge_index_map_t,double> >(pmap);
            else if (get_static_property_map<vector_property_map<size_t,edge_index_map_t> >(&pmap))
                _edge_filter_map = get_static_property_map<vector_property_map<size_t,edge_index_map_t> >(pmap);
            else if (get_static_property_map<HashedDescriptorMap<edge_index_map_t,size_t> >(&pmap))
                _edge_filter_map = get_static_property_map<HashedDescriptorMap<edge_index_map_t,size_t> >(pmap);
            else if (get_static_property_map<edge_index_map_t>(&pmap))
                _edge_filter_map = get_static_property_map<edge_index_map_t>(pmap);
            else 
                _edge_filter_map = DynamicPropertyMapWrap<double, graph_traits<multigraph_t>::edge_descriptor>(pmap);
        }
        catch (property_not_found) 
        {
            throw GraphException("property " + property + " not found");
        }
    }

#else
    if (property != "")
        throw GraphException("support for graph range filtering was not enabled during compilation.");
#endif
}

bool GraphInterface::IsEdgeFilterActive() const {return _edge_filter_property != "" || _edge_python_filter != python::object();}

void GraphInterface::SetVertexFilterRange(std::pair<double,double> allowed_range)
{
#ifndef NO_RANGE_FILTERING
    _vertex_range = allowed_range;
#else
    throw GraphException("support for graph range filtering was not enabled during compilation.");
#endif
}


void GraphInterface::SetGenericVertexFilter(boost::python::object filter) 
{
#ifndef NO_PYTHON_FILTERING
    _vertex_python_filter = filter;
#else
    if (filter != python::object())
        throw GraphException("support for graph python filtering was not enabled during compilation.");
#endif    
}

void GraphInterface::SetGenericEdgeFilter(boost::python::object filter) 
{ 
#ifndef NO_PYTHON_FILTERING
    _edge_python_filter = filter;
#else
    if (filter != python::object())
        throw GraphException("support for graph python filtering was not enabled during compilation.");
#endif    
}


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
        SharedMap<Hist> s_hist(hist);
        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) firstprivate(s_hist) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;            
            s_hist[_degree(v,g)]++;
        }
        s_hist.Gather();
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
        SharedMap<Hist> s_hist(hist);

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) firstprivate(s_hist) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            typename graph_traits<Graph>::out_edge_iterator e, e_begin, e_end;
            tie(e_begin,e_end) = out_edges(v,g);
            for(e = e_begin; e != e_end; ++e)
                if(is_convertible<typename graph_traits<Graph>::directed_category, undirected_tag>::value)
                    s_hist[_prop(*e,g)] += 0.5;
                else
                    s_hist[_prop(*e,g)]++;
        }
        s_hist.Gather();
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
// LabelComponents(property)
//==============================================================================

struct label_components
{
    template <class Graph, class CompMap>
    void operator()(const Graph &g, CompMap comp_map) const
    {   
        get_components(g, comp_map, typename is_convertible<typename graph_traits<Graph>::directed_category, directed_tag>::type());
    }

    template <class Graph, class CompMap>
    void get_components(const Graph& g, CompMap comp_map, boost::true_type is_directed) const
    {
        strong_components(g, comp_map);
    }

    template <class Graph, class CompMap>
    void get_components(const Graph& g, CompMap comp_map, boost::false_type is_directed) const
    {
        connected_components(g, comp_map);
    }
    
};

void GraphInterface::LabelComponents(string prop)
{
    typedef HashedDescriptorMap<vertex_index_map_t, size_t> comp_map_t;
    static comp_map_t comp_map(_vertex_index);

    check_filter(*this, bind<void>(label_components(), _1, comp_map), reverse_check(), directed_check());

    try
    {
        find_property_map(_properties, prop, typeid(graph_traits<multigraph_t>::vertex_descriptor));
        RemoveVertexProperty(prop);
    }
    catch (property_not_found) {}

    _properties.property(prop, comp_map);
}

//==============================================================================
// label_parallel_edges
// label parallel edges in order
//==============================================================================
struct label_parallel_edges
{
    template <class Graph, class EdgeIndexMap, class ParallelMap>
    void operator()(const Graph& g, EdgeIndexMap edge_index, ParallelMap parallel) const
    {
        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            tr1::unordered_set<typename graph_traits<Graph>::edge_descriptor,DescriptorHash<EdgeIndexMap> > p_edges(0, DescriptorHash<EdgeIndexMap>(edge_index));
            typename graph_traits<Graph>::out_edge_iterator e1, e2, e_end;
            for (tie(e1, e_end) = out_edges(v, g); e1 != e_end; ++e1)
            {
                if (p_edges.find(*e1) != p_edges.end())
                    continue;
                size_t n = 0;
                put(parallel, *e1, n);
                for (tie(e2, e_end) = out_edges(v, g); e2 != e_end; ++e2)
                    if (*e2 != *e1 && target(*e1, g) == target(*e2, g))
                    {
                        put(parallel, *e2, ++n);
                        p_edges.insert(*e2);
                    }
            }
        }
    }
};

//==============================================================================
// LabelParallelEdges(property)
//==============================================================================
void GraphInterface::LabelParallelEdges(string property)
{
    try
    {
        DynamicPropertyMapWrap<size_t,graph_traits<multigraph_t>::edge_descriptor> parallel_map(find_property_map(_properties, property, typeid(graph_traits<multigraph_t>::edge_descriptor)));
        check_filter(*this, bind<void>(label_parallel_edges(), _1, _edge_index, parallel_map), reverse_check(), directed_check());
    }
    catch (property_not_found) 
    {
        typedef HashedDescriptorMap<edge_index_map_t,size_t> parallel_map_t;
        parallel_map_t parallel_map(_edge_index);
        check_filter(*this, bind<void>(label_parallel_edges(), _1, _edge_index, parallel_map), reverse_check(), directed_check());
        _properties.property(property, parallel_map);
    }
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
// ComputeGraphLayoutGursoy(iter, seed)
//==============================================================================

struct compute_gursoy
{
    template <class Graph, class PosMap, class IndexMap>
    void operator()(Graph &g, size_t iter, size_t seed, PosMap pos, IndexMap index_map) const
    {
        mt19937 rng(static_cast<mt19937::result_type>(seed));
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
  
        mt19937 rng(static_cast<mt19937::result_type>(seed));
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
