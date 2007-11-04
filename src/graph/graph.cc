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

#include <algorithm>
#include <tr1/unordered_set>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/connected_components.hpp>
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

// this will return the degree type contained in 'degree' or the name of the
// scalar property.
pair<GraphInterface::degree_t,string>
graph_tool::get_degree_type(GraphInterface::deg_t degree)
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

// this is the constructor for the graph interface
GraphInterface::GraphInterface()
    :_mg(),
     _reversed(false),
     _directed(true),
     _vertex_index(get(vertex_index,_mg)),
     _edge_index(get(edge_index_t(),_mg)),
     _vertex_filter_map(_vertex_index),
     _vertex_range(make_pair(0.0, numeric_limits<double>::max())),
     _vertex_range_include(make_pair(true, true)),
     _vertex_range_invert(false),
     _edge_filter_map(_edge_index),
     _edge_range(make_pair(0.0, numeric_limits<double>::max())),
     _edge_range_include(make_pair(true, true)),
     _edge_range_invert(false)
{

}

// the destructor
GraphInterface::~GraphInterface()
{
}

// these test whether or not the vertex and filter are active, by checking the
// corresponding property value
bool GraphInterface::IsVertexFilterActive() const
{ return _vertex_filter_property != ""; }

bool GraphInterface::IsEdgeFilterActive() const
{ return _edge_filter_property != ""; }

// this function will get a dynamic property map, and return its static type, if
// it matches a couple of possibilities that have 'double' as value
// type. Otherwise it wraps it around a DynamicPropertyMapWrap
template <class FilterMap, class Descriptor, class IndexMap>
FilterMap choose_filter_map(string property, dynamic_properties& dp)
{
    FilterMap filter_map;
    try
    {
        dynamic_property_map& pmap =
            find_property_map(dp, property, typeid(Descriptor));

        try
        {
            return get_static_property_map
                <vector_property_map<double,IndexMap> >(pmap);
        }
        catch(bad_cast){}

        try
        {
            return get_static_property_map
                <HashedDescriptorMap<IndexMap,double> >(pmap);
        }
        catch(bad_cast){}

        try
        {
            return get_static_property_map
                <vector_property_map<size_t,IndexMap> >(pmap);
        }
        catch(bad_cast){}

        try
        {
            return get_static_property_map
                <HashedDescriptorMap<IndexMap,size_t> >(pmap);
        }
        catch(bad_cast){}

        try
        {
            return get_static_property_map<IndexMap>(pmap);
        }
        catch(bad_cast){}

        return DynamicPropertyMapWrap<double,Descriptor>(pmap);
    }
    catch (property_not_found)
    {
        throw GraphException("property " + property + " not found");
    }

}

// this will set enable the vertex filter, according to the given property. If
// property is an empty string, this will disable the filter
void GraphInterface::SetVertexFilterProperty(string property)
{
#ifndef NO_RANGE_FILTERING
    _vertex_filter_property = property;

    if (property != "")
    {
        _vertex_filter_map = choose_filter_map<vertex_filter_map_t,
                                               vertex_t,vertex_index_map_t>
            (property, _properties);
    }
    else
    {
        _vertex_filter_map = _vertex_index;
        _vertex_range = make_pair(0.0, numeric_limits<double>::max());
        _vertex_range_include = make_pair(true, true);
        _vertex_range_invert = false;
    }
#else
    if (property != "")
        throw GraphException("support for graph range filtering was not enabled"
                             " during compilation.");
#endif
}

// this will set enable the edge filter, according to the given property. If
// property is an empty string, this will disable the filter
void GraphInterface::SetEdgeFilterProperty(string property)
{
#ifndef NO_RANGE_FILTERING
    _edge_filter_property = property;

    if (property != "")
    {
        _edge_filter_map = choose_filter_map<edge_filter_map_t,
                                               edge_t,edge_index_map_t>
            (property, _properties);
    }
    else
    {
        _edge_filter_map = _edge_index;
        _edge_range = make_pair(0.0, numeric_limits<double>::max());
        _edge_range_include = make_pair(true, true);
        _edge_range_invert = false;
    }
#else
    if (property != "")
        throw GraphException("support for graph range filtering was not enabled"
                             " during compilation.");
#endif
}


void GraphInterface::SetVertexFilterRange(pair<double,double> allowed_range,
                                          pair<bool,bool> include, bool invert)
{
#ifndef NO_RANGE_FILTERING
    if (isinf(allowed_range.first) == 1)
        allowed_range.first = numeric_limits<double>::max();
    else if (isinf(allowed_range.first) == -1)
        allowed_range.first = -numeric_limits<double>::max();
    if (isinf(allowed_range.second) == 1)
        allowed_range.second = numeric_limits<double>::max();
    else if (isinf(allowed_range.second) == -1)
        allowed_range.second = -numeric_limits<double>::max();

    _vertex_range = allowed_range;
    _vertex_range_include = include;
    _vertex_range_invert = invert;
#else
    throw GraphException("support for graph range filtering was not enabled"
                         " during compilation.");
#endif
}

void GraphInterface::SetEdgeFilterRange(pair<double,double> allowed_range,
                                        pair<bool,bool> include, bool invert)
{
#ifndef NO_RANGE_FILTERING
    if (isinf(allowed_range.first) == 1)
        allowed_range.first = numeric_limits<double>::max();
    else if (isinf(allowed_range.first) == -1)
        allowed_range.first = -numeric_limits<double>::max();
    if (isinf(allowed_range.second) == 1)
        allowed_range.second = numeric_limits<double>::max();
    else if (isinf(allowed_range.second) == -1)
        allowed_range.second = -numeric_limits<double>::max();

    _edge_range = allowed_range;
    _edge_range_include = include;
    _edge_range_invert = invert;
#else
    throw GraphException("support for graph range filtering was not enabled"
                         " during compilation.");
#endif
}


// this function will reindex all the edges, in the order in which they are
// found, taking care of preserving all the properties
template<class Graph, class EdgeIndex>
void reindex_edges(Graph &g, EdgeIndex edge_index, dynamic_properties& dp)
{
    typedef typename graph_traits<Graph>::edge_descriptor edge_t;

    size_t n_edges = num_edges(g);
    if (n_edges == 0)
        return;
    vector<pair<edge_t,bool> > edge_map
        (n_edges, make_pair(edge_t(), false));
    typename graph_traits<Graph>::vertex_iterator v, v_end;
    typename graph_traits<Graph>::out_edge_iterator e, e_end;
    for (tie(v, v_end) = vertices(g); v != v_end; ++v)
        for (tie(e, e_end) = out_edges(*v, g); e != e_end; ++e)
        {
            size_t index = edge_index[*e];
            if (index >= edge_map.size())
                edge_map.resize(index+1);
            edge_map[index] = make_pair(*e, true);
        }
    size_t new_index = 0;
    for (tie(v, v_end) = vertices(g); v != v_end; ++v)
        for (tie(e, e_end) = out_edges(*v, g); e != e_end; ++e)
        {
            edge_t old_edge = edge_map[new_index].first;
            if (edge_map[new_index].second)
            {
                // another edge exists with the same index; indexes
                // must be swapped, as well as the properties
                edge_index[old_edge] = edge_index[*e];
                edge_map[edge_index[*e]] = make_pair(old_edge, true);
                edge_index[*e] = new_index;
                edge_map[new_index] = make_pair(*e, true);
                for (typeof(dp.begin()) p = dp.begin(); p != dp.end(); ++p)
                    if (p->second->key() == typeid(edge_t))
                    {
                        boost::any temp = p->second->get(*e);
                        p->second->put(*e, p->second->get(old_edge));
                        p->second->put(old_edge, temp);
                    }
            }
            else
            {
                // no other edge has this index; it must be then
                // assigned for this edge, and the properties must be
                // copied over
                size_t old_index = edge_index[*e];
                for (typeof(dp.begin()) p = dp.begin(); p != dp.end(); ++p)
                    if (p->second->key() == typeid(edge_t))
                    {
                        edge_index[*e] = old_index;
                        boost::any val = p->second->get(*e);
                        edge_index[*e] = new_index;
                        p->second->put(*e, val);
                    }
            }
            ++new_index;
        }
}

// this will definitively remove all the edges from the graph, which are being
// currently filtered out. This will also disable the edge filter
void GraphInterface::PurgeEdges()
{
    if (!IsEdgeFilterActive())
        return;

    RangeFilter<edge_filter_map_t> filter(_edge_filter_map, _edge_range,
                                          _edge_range_include,
                                          _edge_range_invert);
    graph_traits<multigraph_t>::vertex_iterator v, v_end;
    graph_traits<multigraph_t>::out_edge_iterator e, e_end;
    vector<graph_traits<multigraph_t>::edge_descriptor> deleted_edges;
    for (tie(v, v_end) = vertices(_mg); v != v_end; ++v)
    {
        for (tie(e, e_end) = out_edges(*v, _mg); e != e_end; ++e)
            if (!filter(*e))
                deleted_edges.push_back(*e);
        for (typeof(deleted_edges.begin()) iter = deleted_edges.begin();
             iter != deleted_edges.end(); ++iter)
            remove_edge(*iter, _mg);
        deleted_edges.clear();
    }
    reindex_edges(_mg, _edge_index, _properties);
}

// this will definitively remove all the verticess from the graph, which are being
// currently filtered out. This will also disable the vertex filter
void GraphInterface::PurgeVertices()
{
    if (!IsVertexFilterActive())
        return;

    RangeFilter<vertex_filter_map_t> filter(_vertex_filter_map, _vertex_range,
                                            _vertex_range_include,
                                            _vertex_range_invert);
    size_t N = num_vertices(_mg);
    vector<bool> deleted(N, false);
    for (size_t i = 0; i < N; ++i)
        deleted[i] = !filter(vertex(i, _mg));

    //migrate properties
    size_t last = 0;
    for (size_t i = 0; i < N; ++i)
    {
        graph_traits<multigraph_t>::vertex_descriptor v = vertex(i, _mg);
        if (filter(v))
        {
            for (typeof(_properties.begin()) p = _properties.begin();
                 p != _properties.end(); ++p)
                if (p->second->key() == typeid(vertex_t))
                    try
                    {
                        p->second->put(vertex(last, _mg), p->second->get(v));
                    }
                    catch (dynamic_const_put_error) {} // index prop. is const
            last++;
        }
    }

    //remove vertices
    for (size_t i = N-1; i >= 0 && i < N; --i)
    {
        if (deleted[i])
        {
            graph_traits<multigraph_t>::vertex_descriptor v = vertex(i, _mg);
            clear_vertex(v, _mg);
            remove_vertex(v, _mg);
        }
    }
    reindex_edges(_mg, _edge_index, _properties);
}

// this will get the number of vertices, either the "soft" O(1) way, or the hard
// O(V) way, which is necessary if the graph is filtered
size_t GraphInterface::GetNumberOfVertices() const
{
    size_t n = 0;
    if (IsVertexFilterActive())
        check_filter(*this,var(n)=bind<size_t>(HardNumVertices(),_1),
                     reverse_check(),directed_check());
    else
        check_filter(*this,var(n)=bind<size_t>(SoftNumVertices(),_1),
                     reverse_check(),directed_check());
    return n;
}

// this will get the number of edges, either the "soft" O(E) way, or the hard
// O(E) way, which is necessary if the graph is filtered. Both cases are of
// linear complexity, since num_edges() is O(E) in Boost's adjacency_list
size_t GraphInterface::GetNumberOfEdges() const
{
    size_t n = 0;
    if (IsEdgeFilterActive() || IsVertexFilterActive())
        check_filter(*this,var(n)=bind<size_t>(HardNumEdges(),_1),
                     reverse_check(),directed_check());
    else
        check_filter(*this,var(n)=bind<size_t>(SoftNumEdges(),_1),
                     reverse_check(),directed_check());
    return n;
}


// generalized functor to obtain histogram of different types of "degrees"
template <class DegreeSelector>
struct get_vertex_histogram
{
    get_vertex_histogram(DegreeSelector& deg): _degree(deg) {}
    template <class Graph, class Hist>
    void operator()(const Graph& g, Hist& hist) const
    {
        SharedMap<Hist> s_hist(hist);
        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) \
            firstprivate(s_hist) schedule(dynamic)
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

// used to cycle through the different types of histograms
struct choose_vertex_histogram
{
    choose_vertex_histogram(const GraphInterface& g, GraphInterface::deg_t deg,
                            GraphInterface::hist_t& hist)
        : _g(g), _hist(hist)
    {
        tie(_deg, _deg_name) = get_degree_type(deg);
    }
    template <class DegreeSelector>
    void operator()(DegreeSelector)
    {
        if (mpl::at<degree_selector_index,DegreeSelector>::type::value == _deg)
        {
            DegreeSelector selector(_deg_name, _g, true);
            check_filter(_g, bind<void>(get_vertex_histogram<DegreeSelector>
                                        (selector), _1, var(_hist)),
                         reverse_check(),directed_check());
        }
    }
    const GraphInterface &_g;
    GraphInterface::hist_t &_hist;
    GraphInterface::degree_t _deg;
    string _deg_name;
};


// this will return the vertex histogram of degrees or scalar properties
GraphInterface::hist_t
GraphInterface::GetVertexHistogram(GraphInterface::deg_t deg) const
{
    hist_t hist;
    try
    {
        mpl::for_each<mpl::vector<in_degreeS,
                                  out_degreeS,
                                  total_degreeS,
                                  scalarS> >
            (choose_vertex_histogram(*this, deg, hist));
    }
    catch (dynamic_get_failure &e)
    {
        throw GraphException("error getting scalar property: " +
                             string(e.what()));
    }

    return hist;
}

// generalized functor to obtain histogram of edge properties
struct get_edge_histogram
{
    get_edge_histogram(scalarS& prop): _prop(prop) {}
    template <class Graph, class Hist>
    void operator()(const Graph& g, Hist& hist) const
    {
        typedef typename graph_traits<Graph>::directed_category
            directed_category;
        SharedMap<Hist> s_hist(hist);

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) \
            firstprivate(s_hist)  schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            typename graph_traits<Graph>::out_edge_iterator e, e_begin, e_end;
            tie(e_begin,e_end) = out_edges(v,g);
            for(e = e_begin; e != e_end; ++e)
                if(is_convertible<directed_category,undirected_tag>::value)
                    s_hist[_prop(*e,g)] += 0.5;
                else
                    s_hist[_prop(*e,g)]++;
        }
        s_hist.Gather();
    }
    scalarS& _prop;
};



// returns the histogram of a given edge property
GraphInterface::hist_t GraphInterface::GetEdgeHistogram(string property) const
{
    hist_t hist;
    try
    {
        scalarS prop(property, *this, false);
        check_filter(*this, bind<void>(get_edge_histogram(prop), _1, var(hist)),
                     reverse_check(),directed_check());
    }
    catch (dynamic_get_failure &e)
    {
        throw GraphException("error getting scalar property: " +
                             string(e.what()));
    }

    return hist;
}

// this will label the components of a graph to a given vertex property, from
// [0, number of components - 1]. If the graph is directed the "strong
// components" are used.
struct label_components
{
    template <class Graph, class CompMap>
    void operator()(const Graph &g, CompMap comp_map) const
    {
        typedef typename graph_traits<Graph>::directed_category
            directed_category;
        get_components(g, comp_map,
                       typename is_convertible<directed_category,
                                               directed_tag>::type());
    }

    template <class Graph, class CompMap>
    void get_components(const Graph& g, CompMap comp_map,
                        boost::true_type is_directed) const
    {
        strong_components(g, comp_map);
    }

    template <class Graph, class CompMap>
    void get_components(const Graph& g, CompMap comp_map,
                        boost::false_type is_directed) const
    {
        connected_components(g, comp_map);
    }

};

void GraphInterface::LabelComponents(string prop)
{
    typedef vector_property_map<size_t, vertex_index_map_t> comp_map_t;
    comp_map_t comp_map(_vertex_index);

    check_filter(*this, bind<void>(label_components(), _1, comp_map),
                 reverse_check(), directed_check());

    try
    {
        find_property_map(_properties, prop, typeid(vertex_t));
        RemoveVertexProperty(prop);
    }
    catch (property_not_found) {}

    _properties.property(prop, comp_map);
}

// label parallel edges in the order they are found, starting from 0
struct label_parallel_edges
{
    template <class Graph, class EdgeIndexMap, class ParallelMap>
    void operator()(const Graph& g, EdgeIndexMap edge_index,
                    ParallelMap parallel) const
    {
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            tr1::unordered_set<edge_t,DescriptorHash<EdgeIndexMap> >
                p_edges(0, DescriptorHash<EdgeIndexMap>(edge_index));

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

void GraphInterface::LabelParallelEdges(string property)
{
    try
    {
        DynamicPropertyMapWrap<size_t,edge_t>
            parallel_map(find_property_map(_properties, property,
                                           typeid(edge_t)));
        check_filter(*this, bind<void>(label_parallel_edges(), _1, _edge_index,
                                       parallel_map),
                     reverse_check(), directed_check());
    }
    catch (property_not_found)
    {
        typedef vector_property_map<size_t,edge_index_map_t> parallel_map_t;
        parallel_map_t parallel_map(_edge_index);
        check_filter(*this, bind<void>(label_parallel_edges(), _1, _edge_index,
                                       parallel_map),
                     reverse_check(), directed_check());
        _properties.property(property, parallel_map);
    }
}


// this inserts the edge index map as a property
void GraphInterface::InsertEdgeIndexProperty(string property)
{
    _properties.property(property, _edge_index);
}

// this inserts the vertex index map as a property
void GraphInterface::InsertVertexIndexProperty(string property)
{
    _properties.property(property, _vertex_index);
}


// signal handling

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
        std::cerr << "SIGSEGV (Segmentation fault). "
                  << "There's a bug somewhere in the program. Go fix it.";
        break;
    case SIGBUS:
        std::cerr << "SIGBUS (BUS error). The bus is broken, I guess...";
        break;
    case SIGFPE:
        std::cerr << "SIGFPE (Floating-point exception). "
                  << "INFs and NANs are wreaking havoc.";
        break;
    case SIGILL:
        std::cerr << "SIGILL (Illegal instruction). Did you compile it right?";
        break;
    case SIGXCPU:
        std::cerr << "SIGXCPU (CPU limit exceeded). Time's over.";
        break;
    case SIGXFSZ:
        std::cerr << "SIGXFSZ (File size limit exceeded). "
                  << "The fascist sysadmin doesn't let you write more.";
        break;
    }
    std::cerr << " Bailing out cowardly and calling abort()." << std::endl;
    abort();
}

// initialize signal handling

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
