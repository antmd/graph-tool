#include <ext/hash_map>
#include <map>
#include <iostream>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/connected_components.hpp>


using namespace boost;

int main()
{
    typedef  adjacency_list
                < listS,
                  listS,
                  bidirectionalS,
                  property<vertex_index_t, int, 
                    property<vertex_index2_t, int> > > Graph;

    function_requires<VertexListGraphConcept<Graph> >();
    function_requires<EdgeListGraphConcept<Graph > >();  
    function_requires<VertexMutableGraphConcept<Graph > >();
    function_requires<EdgeMutableGraphConcept<Graph> >();
    function_requires<MutableIncidenceGraphConcept<Graph> >();
    function_requires<AdjacencyGraphConcept<Graph> >();
    function_requires<VertexMutablePropertyGraphConcept<Graph> >();
    function_requires<EdgeMutablePropertyGraphConcept<Graph> >();
       
    Graph g;
    
    typedef graph_traits<Graph>::vertex_descriptor Vertex;
    typedef graph_traits<Graph>::edge_descriptor Edge;
    
    Vertex v[9];
    for (int i = 0; i < 9; ++i)
        v[i] = add_vertex(i, g);
    for (int i = 0; i < 8; ++i)
        add_edge(v[i],v[i+1],g);
      
    add_edge(v[8],v[6],g);
    add_edge(v[6],v[7],g);
    add_edge(v[7],v[1],g);
    add_edge(v[1],v[2],g);
    add_edge(v[2],v[4],g);
    add_edge(v[4],v[5],g);

    vertex_index2_t vertex_index2;    
    strong_components(g, get(vertex_index2, g));
    
    std::cout << "Directed:" << std::endl;
    vertex_index_t vertex_index;
    for(graph_traits<Graph>::vertex_iterator iter = vertices(g).first; 
        iter != vertices(g).second; ++iter)
    {
        graph_traits<Graph>::vertex_descriptor v = *iter;
        std::cout << get(vertex_index, g, v) 
                  << "(" << get(vertex_index2, g, v) << ") : ";

        for(graph_traits<Graph>::out_edge_iterator iter2 = out_edges(v,g).first;
            iter2 != out_edges(v,g).second; ++iter2)
            std::cout << get(vertex_index, g, target(*iter2,g)) << " ";
        std::cout << std::endl;
    }

    
}
