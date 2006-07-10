#include <ext/hash_map>
#include <map>
#include <iostream>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_concepts.hpp>
#include "../filtered_graph.hpp"
#include "../graph_filter.hh"
#include "../graph_adaptor.hh"

using namespace boost;


template <class Graph>
void print_graph(Graph &g)
{
    for(typename graph_traits<Graph>::vertex_iterator iter = vertices(g).first; 
        iter != vertices(g).second; ++iter)
    {
        typename graph_traits<Graph>::vertex_descriptor v = *iter;
        std::cout << g[v].index << ": ";

        for(typename graph_traits<Graph>::out_edge_iterator iter2 = 
                                                           out_edges(v,g).first;
            iter2 != out_edges(v,g).second; ++iter2)
            std::cout << g[target(*iter2,g)].index
                      << "(" << g[*iter2].hidden << ") ";
        std::cout << std::endl;
    }
}

template <class Graph>
void print_graph_adjacency(Graph &g)
{
    for(typename graph_traits<Graph>::vertex_iterator iter = vertices(g).first; 
        iter != vertices(g).second; ++iter)
    {
        typename graph_traits<Graph>::vertex_descriptor v = *iter;
        std::cout << g[v].index << ": ";
        for(typename graph_traits<Graph>::adjacency_iterator 
                                       iter2 = adjacent_vertices(*iter,g).first;
            iter2 != adjacent_vertices(*iter,g).second; ++iter2)
            std::cout << g[*iter2].index <<" ";
        std::cout << std::endl;
    }
}


struct VertexProperties
{ 
    size_t index;
    bool hidden;
};

struct EdgeProperties
{
    size_t index;
    bool hidden;

};

template <class Graph, class VertexOrEdge>
class HiddenFilter
{
    public:
    HiddenFilter(): _g(0) {}
    HiddenFilter(const Graph &g):_g(&g) {}
    
    bool operator() (VertexOrEdge u) const 
    {       
        return (*_g)[u].hidden;
    }
private:
    Graph const *_g;
};


int main()
{
    
    typedef  adjacency_list
                < boost::vecS, // boost::vecS, boost::hash_setS
                  boost::setS, //boost::slistS, 
                  boost::bidirectionalS,
                  VertexProperties,
                  EdgeProperties > Graph;

    typedef graph_traits<Graph>::vertex_descriptor Vertex;
    typedef graph_traits<Graph>::edge_descriptor Edge;
    
    typedef HiddenFilter<Graph, Vertex> vpred_t;
    typedef HiddenFilter<Graph, Edge> epred_t;

    function_requires<VertexListGraphConcept<UncoupledAdaptor<Graph, epred_t, vpred_t> > >();
    function_requires<EdgeListGraphConcept<UncoupledAdaptor<Graph, epred_t, vpred_t> > >();  
    function_requires<VertexMutableGraphConcept<UncoupledAdaptor<Graph, epred_t, vpred_t> > >();
    function_requires<EdgeMutableGraphConcept<UncoupledAdaptor<Graph, epred_t, vpred_t> > >();
    function_requires<MutableIncidenceGraphConcept<UncoupledAdaptor<Graph, epred_t, vpred_t> > >();
    function_requires<AdjacencyGraphConcept<UncoupledAdaptor<Graph, epred_t, vpred_t> > >();
    function_requires<VertexMutablePropertyGraphConcept<UncoupledAdaptor<Graph, epred_t, vpred_t> > >();
    function_requires<EdgeMutablePropertyGraphConcept<UncoupledAdaptor<Graph, epred_t, vpred_t> > >();
    
    Graph g;
    
    
    Vertex v1 = add_vertex(g);
    g[v1].index = 1;
    g[v1].hidden = true;
    Vertex v2 = add_vertex(g);
    g[v2].index = 2;
    g[v2].hidden = false;
    Vertex v3 = add_vertex(g);
    g[v3].index = 3;
    g[v3].hidden = true;
    Vertex v4 = add_vertex(g);
    g[v4].index = 4;
    g[v4].hidden = false;
    Vertex v5 = add_vertex(g);
    g[v5].index = 5;
    g[v5].hidden = true;
    
    Edge e;
    e = add_edge(v1,v2,g).first;
    g[e].index = 1;
    g[e].hidden = true;
    e = add_edge(v2,v3,g).first;
    g[e].index = 2;
    g[e].hidden = false;
    e = add_edge(v3,v4,g).first;
    g[e].index = 3;
    g[e].hidden = false;
    e = add_edge(v4,v5,g).first;
    g[e].index = 4;
    g[e].hidden = true;
    e = add_edge(v5,v1,g).first;
    g[e].index = 5;
    g[e].hidden = true;
    e = add_edge(v5,v1,g).first;
    g[e].index = 6;
    g[e].hidden = true;

    
    std::cout << "coupled:" << std::endl;
    print_graph(g);

    vpred_t vpred(g);
    epred_t epred(g);
    typedef UncoupledAdaptor<Graph, epred_t, vpred_t> Ugraph; 
    Ugraph ug(g, epred, vpred);
    
    std::cout << std::endl << "uncoupled:" << std::endl;
    print_graph(ug);
    
    std::cout << std::endl << "uncoupled (using adjacent_vertices):" 
                           << std::endl;
    print_graph_adjacency(ug);
    std::cout << std::endl;
    
    ParallelFilter<Graph> directed_filter(g);
    filtered_graph<Graph, ParallelFilter<Graph> > gf(g,directed_filter);
    
    std::cout << std::endl<< "filtered coupled:" << std::endl;
    print_graph(gf);
    
    ParallelFilter<Ugraph> undirected_filter(ug);
    filtered_graph<Ugraph, ParallelFilter<Ugraph> > ugf(ug,undirected_filter);
   
    std::cout << std::endl<< "filtered uncoupled:" << std::endl;
    print_graph(ugf);

    typedef graph_traits<filtered_graph<Graph, ParallelFilter<Graph> > >::vertex_descriptor FVertex;
    typedef graph_traits<filtered_graph<Graph, ParallelFilter<Graph> > >::edge_descriptor FEdge;
    
    typedef HiddenFilter<filtered_graph<Graph, ParallelFilter<Graph> >, FVertex> fvpred_t;
    typedef HiddenFilter<filtered_graph<Graph, ParallelFilter<Graph> >, FEdge> fepred_t;


    fvpred_t fvpred(gf);
    fepred_t fepred(gf);
    typedef UncoupledAdaptor<filtered_graph<Graph, ParallelFilter<Graph> >, 
                             fepred_t, fvpred_t> FUgraph; 
    FUgraph fug(gf, fepred, fvpred);
    
    std::cout << std::endl<< "uncoupled filtered:" << std::endl;
    print_graph(fug);
    
    ParallelFilter<FUgraph> filter_undirected_filter(fug);
    filtered_graph<FUgraph, ParallelFilter<FUgraph> > fugf(fug, filter_undirected_filter);
   
    std::cout << std::endl<< "filtered uncoupled filtered:" << std::endl;
    print_graph(fugf);
    

}
