// Copyright (C) 2008  Tiago de Paula Peixoto <tiago@forked.de>
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

// generalized functor to obtain histogram of different types of "degrees"
struct get_vertex_histogram
{
    template <class Graph, class DegreeSelector, class Hist>
    void operator()(const Graph* gp, DegreeSelector deg, Hist& hist) const
    {
        const Graph& g = *gp;
        SharedMap<Hist> s_hist(hist);
        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) \
            firstprivate(s_hist) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            s_hist[deg(v,g)]++;
        }
        s_hist.Gather();
    }
};


// this will return the vertex histogram of degrees or scalar properties
hist_t
GraphInterface::GetVertexHistogram(GraphInterface::deg_t deg) const
{
    hist_t hist;

    try
    {
        run_action<>()(*this, bind<void>(get_vertex_histogram(), _1, _2,
                                         var(hist)),
                       all_selectors())(degree_selector(deg, _properties));
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
    template <class Graph, class Prop, class Hist>
    void operator()(const Graph* gp, Prop eprop, Hist& hist) const
    {
        const Graph& g = *gp;
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
                s_hist[eprop[*e]]++;
        }
        s_hist.Gather();

        typedef typename graph_traits<Graph>::directed_category
            directed_category;
        if(is_convertible<directed_category,undirected_tag>::value)
        {
            for (typeof(s_hist.begin()) iter = s_hist.begin();
                 iter != s_hist.end(); ++iter)
                iter->second /= typename Hist::value_type::second_type(2);
        }
    }
};


// returns the histogram of a given edge property
hist_t GraphInterface::GetEdgeHistogram(string property) const
{
    hist_t hist;
    try
    {
        run_action<>()(*this, bind<void>(get_edge_histogram(), _1, _2,
                                         var(hist)),
                       edge_scalar_properties())
            (prop(property, _edge_index, _properties));
    }
    catch (dynamic_get_failure& e)
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
    void operator()(const Graph* gp, CompMap comp_map) const
    {
        const Graph& g = *gp;
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
    try
    {
        run_action<>()(*this, label_components(), _1, vertex_scalar_selectors())
            (degree_selector(prop, _properties));
    }
    catch (property_not_found)
    {
        typedef vector_property_map<int32_t, vertex_index_map_t> comp_map_t;
        comp_map_t comp_map(_vertex_index);
        _properties.property(prop, comp_map);
        run_action<>()(*this, bind<void>(label_components(), _1, comp_map))();
    }
}

// label parallel edges in the order they are found, starting from 0
struct label_parallel_edges
{
    template <class Graph, class EdgeIndexMap, class ParallelMap>
    void operator()(const Graph* gp, EdgeIndexMap edge_index,
                    ParallelMap parallel) const
    {
        const Graph& g = *gp;
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
        run_action<>()(*this, bind<void>(label_parallel_edges(), _1,
                                         _edge_index, _2),
                       edge_scalar_properties())
            (prop(property, _edge_index, _properties));
    }
    catch (property_not_found)
    {
        typedef vector_property_map<int32_t,edge_index_map_t> parallel_map_t;
        parallel_map_t parallel_map(_edge_index);
        _properties.property(property, parallel_map);
        run_action<>()(*this, bind<void>(label_parallel_edges(), _1,
                                         _edge_index, parallel_map))();
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

