#include <ext/hash_map>
#include <map>
#include <iostream>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/filtered_graph.hpp>
#include "../graph_adaptor.hh"

using namespace boost;


template <class Graph>
void print_graph(Graph &g)
{
    boost::vertex_index_t vertex_index;
    boost::edge_weight_t edge_weight;
    for(typename graph_traits<Graph>::vertex_iterator iter = vertices(g).first; 
        iter != vertices(g).second; ++iter)
    {
        typename graph_traits<Graph>::vertex_descriptor v = *iter;
        std::cout << get(vertex_index, g, v) << ": ";

        for(typename graph_traits<Graph>::out_edge_iterator iter2 = out_edges(v,g).first;
            iter2 != out_edges(v,g).second; ++iter2)
            std::cout << get(vertex_index, g, target(*iter2,g)) 
                      << "(" << get(edge_weight, g, *iter2) << ") ";
        std::cout << std::endl;
    }
}

template <class Graph>
void print_graph_adjacency(Graph &g)
{
    boost::vertex_index_t vertex_index;
    for(typename graph_traits<Graph>::vertex_iterator iter = vertices(g).first; 
        iter != vertices(g).second; ++iter)
    {
        typename graph_traits<Graph>::vertex_descriptor v = *iter;
        std::cout << get(vertex_index, g, v) << ": ";
        for(typename graph_traits<Graph>::adjacency_iterator 
                                       iter2 = adjacent_vertices(*iter,g).first;
            iter2 != adjacent_vertices(*iter,g).second; ++iter2)
            std::cout << get(vertex_index, g, *iter2) <<" ";
        std::cout << std::endl;
    }
}

int main()
{
    typedef  adjacency_list
                < boost::vecS, // boost::vecS, boost::hash_setS
                  boost::vecS, //boost::slistS, 
                  boost::bidirectionalS,
                  boost::property<boost::vertex_index_t, int >,
                  boost::property<boost::edge_weight_t, int > > Graph;

    function_requires<VertexListGraphConcept<UndirectedAdaptor<Graph> > >();
    function_requires<EdgeListGraphConcept<UndirectedAdaptor<Graph> > >();  
    function_requires<VertexMutableGraphConcept<UndirectedAdaptor<Graph> > >();
    function_requires<EdgeMutableGraphConcept<UndirectedAdaptor<Graph> > >();
    function_requires<MutableIncidenceGraphConcept<UndirectedAdaptor<Graph> > >();
    function_requires<AdjacencyGraphConcept<UndirectedAdaptor<Graph> > >();
    function_requires<VertexMutablePropertyGraphConcept<UndirectedAdaptor<Graph> > >();
    function_requires<EdgeMutablePropertyGraphConcept<UndirectedAdaptor<Graph> > >();
    
    Graph g;
    
    typedef graph_traits<Graph>::vertex_descriptor Vertex;
    typedef graph_traits<Graph>::edge_descriptor Edge;
    
    Vertex v1 = add_vertex(1, g);
    Vertex v2 = add_vertex(2, g);
    Vertex v3 = add_vertex(3, g);
    Vertex v4 = add_vertex(4, g);
    Vertex v5 = add_vertex(5, g);
    
    add_edge(v1,v2,1,g);
    add_edge(v2,v1,1,g);
    add_edge(v5,v4,1,g);
    add_edge(v4,v3,1,g);
    add_edge(v3,v1,1,g);
    add_edge(v3,v2,1,g);
    add_edge(v3,v3,1,g);
    add_edge(v3,v4,1,g);
    add_edge(v3,v5,1,g);
    
    std::cout << "Directed:" << std::endl;
    print_graph(g);
    
    typedef UndirectedAdaptor<Graph> Ugraph; 
    Ugraph ug(g);
    
    std::cout << std::endl << "Undirected:" << std::endl;
    print_graph(ug);
    
    std::cout << std::endl << "Undirected (using adjacent_vertices):"  << std::endl;
    print_graph_adjacency(ug);
    
    ParallelFilter<Graph> directed_filter(g);
    filtered_graph<Graph, ParallelFilter<Graph> > gf(g,directed_filter);
    
    std::cout << std::endl<< "Filtered directed:" << std::endl;
    print_graph(gf);
    
    ParallelFilter<Ugraph> undirected_filter(ug);
    filtered_graph<Ugraph, ParallelFilter<Ugraph> > ugf(ug,undirected_filter);
    
    std::cout << std::endl<< "Filtered undirected:" << std::endl;
    print_graph(ugf);
    
    filtered_graph<Ugraph, ParallelFilter<Ugraph> >::edge_iterator e,b;
    for (e = edges(ugf).first; e != edges(ugf).second; e++)
    {
        boost::edge_weight_t tag;
        std::cout << get(tag, ugf, *e) << " " ;
    }
    std::cout << std::endl;
    
    std::cout << "Testing iterators..." << std::endl;
    
    for(graph_traits<Ugraph>::vertex_iterator iter = vertices(ug).first; iter != vertices(ug).second; ++iter)
    {
        graph_traits<Ugraph>::vertex_descriptor v = *iter;
        for(graph_traits<Ugraph>::out_edge_iterator iter2 = out_edges(v,ug).first; iter2 != out_edges(v,ug).second; ++iter2)
        {
            iter2++;
            iter2--;
            std::cout << "(" << get(vertex_index, ug, source(*iter2,ug)) << ", " << get(vertex_index, ug, target(*iter2,ug)) << ") ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl << std::endl;
    for(graph_traits<Ugraph>::vertex_iterator iter = vertices(ug).first; iter != vertices(ug).second; ++iter)
    {
        iter++;
        iter--;
        graph_traits<Ugraph>::vertex_descriptor v = *iter;
        for(graph_traits<Ugraph>::out_edge_iterator iter2 = out_edges(v,ug).second;;)
        {
            if (iter2 == out_edges(v,ug).first)
                break;
            iter2--;
            iter2++;
            iter2--;
            std::cout << "(" << get(vertex_index, ug, source(*iter2,ug)) << ", " << get(vertex_index, ug, target(*iter2,ug)) << ") ";
        }
        std::cout << std::endl;
    }
}
