// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2010  Tiago de Paula Peixoto <tiago@forked.de>
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

#include "graph.hh"
#include "graph_filtering.hh"

#include <graph_subgraph_isomorphism.hh>
#include <graph_python_interface.hh>

using namespace graph_tool;
using namespace boost;

template <class Graph1, class Graph2, class Label1, class Label2>
class PropLabelling
{
public:
    PropLabelling(const Graph1& g1, const Graph2& g2,
                  Label1 label1, Label2 label2)
        : _g1(g1), _g2(g2), _label1(label1), _label2(label2) {}

    typedef typename property_traits<Label1>::value_type value_type1;
    typedef typename property_traits<Label2>::value_type value_type2;

    bool operator()(typename graph_traits<Graph1>::vertex_descriptor v1,
                    typename graph_traits<Graph1>::vertex_descriptor v2) const
    {
        if (in_degreeS()(v2, _g2) < in_degreeS()(v1, _g1) ||
            out_degree(v2, _g2) < out_degree(v1, _g1))
            return false;
        return _label1[v1] == _label2[v2];
    }

    bool operator()(typename graph_traits<Graph1>::edge_descriptor e1,
                    typename graph_traits<Graph1>::edge_descriptor e2) const
    {
        return _label1[e1] == _label2[e2];
    }

private:
    const Graph1& _g1;
    const Graph2& _g2;
    Label1 _label1;
    Label2 _label2;
};

struct get_subgraphs
{
    template <class Graph1, class Graph2, class VertexLabel,
              class EdgeLabel>
    void operator()(const Graph1& g1, const Graph2* g2,
                    VertexLabel vertex_label1, boost::any vertex_label2,
                    EdgeLabel edge_label1, boost::any edge_label2,
                    vector<vector<pair<size_t,size_t> > >& F, size_t max_n)
        const
    {
        typedef PropLabelling<Graph1,Graph2,VertexLabel,VertexLabel>
            vlabelling_t;
        typedef PropLabelling<Graph1,Graph2,EdgeLabel,EdgeLabel>
            elabelling_t;

        subgraph_isomorphism
            (g1, *g2, vlabelling_t(g1, *g2, vertex_label1,
                                   any_cast<VertexLabel>(vertex_label2)),
             elabelling_t(g1, *g2, edge_label1,
                          any_cast<EdgeLabel>(edge_label2)), F, max_n);
    }
};

struct get_mapping
{
    template <class Graph1, class Graph2, class EdgeLabel, class VertexMap,
              class EdgeMap, class EdgeIndexMap>
    void operator()(const Graph1& g1, const Graph2* g2, EdgeLabel edge_label1,
                    boost::any edge_label2, vector<pair<size_t, size_t> >& F,
                    VertexMap vmapping, EdgeMap emapping,
                    EdgeIndexMap edge_index2) const
    {
        typedef PropLabelling<Graph1,Graph2,EdgeLabel,EdgeLabel>
            elabelling_t;
        elabelling_t edge_labelling(g1, *g2, edge_label1,
                                    any_cast<EdgeLabel>(edge_label2));
        int i, N = F.size();
        #pragma omp parallel for default(shared) private(i) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            if (vertex(i, g1) == graph_traits<Graph1>::null_vertex())
                continue;
            vmapping[vertex(F[i].first, g1)] = vertex(F[i].second, *g2);
            typename graph_traits<Graph1>::out_edge_iterator e, e_end;
            for (tie(e, e_end) = out_edges(vertex(i, g1), g1); e != e_end; ++e)
            {
                bool found = false;
                typename graph_traits<Graph2>::out_edge_iterator e2, e2_end;
                for (tie(e2, e2_end) = out_edges(vertex(F[i].second, *g2), *g2);
                     e2 != e2_end; ++e2)
                {
                    if (target(*e2, *g2) ==
                        vertex(F[target(*e, g1)].second, *g2) &&
                        edge_labelling(*e, *e2))
                    {
                        emapping[*e] = edge_index2[*e2];
                        found = true;
                    }
                }
                if (!found)
                    throw GraphException("edge not found... "
                                         "can't be isomorphism!!! "
                                         "This is a bug.");
            }
        }
    }
};

struct directed_graph_view_pointers:
    mpl::transform<graph_tool::detail::always_directed,
                   mpl::quote1<add_pointer> >::type {};

struct undirected_graph_view_pointers:
    mpl::transform<graph_tool::detail::never_directed,
                   mpl::quote1<add_pointer> >::type {};


// typedef mpl::push_back<vertex_properties,
//                        ConstantPropertyMap<bool,GraphInterface::vertex_t> >
//     ::type vertex_props_t;

// typedef mpl::push_back<edge_properties,
//                        ConstantPropertyMap<bool,GraphInterface::edge_t> >
//     ::type edge_props_t;

typedef mpl::vector<ConstantPropertyMap<bool,GraphInterface::vertex_t> > vertex_props_t;
typedef mpl::vector<ConstantPropertyMap<bool,GraphInterface::edge_t> > edge_props_t;

void subgraph_isomorphism(GraphInterface& gi1, GraphInterface& gi2,
                          boost::any vertex_label1, boost::any vertex_label2,
                          boost::any edge_label1, boost::any edge_label2,
                          python::list vmapping, python::list emapping,
                          size_t max_n)
{
    if (gi1.GetDirected() != gi2.GetDirected())
        return;

    if (vertex_label1.empty())
    {
        vertex_label1 = vertex_label2 =
            ConstantPropertyMap<bool,GraphInterface::vertex_t>(true);
    }

    if (edge_label1.empty())
    {
        edge_label1 = edge_label2 =
            ConstantPropertyMap<bool,GraphInterface::edge_t>(true);
    }

    vector<vector<pair<size_t,size_t> > > F;

    if (gi1.GetDirected())
    {
        run_action<graph_tool::detail::always_directed>()
            (gi1, bind<void>(get_subgraphs(),
                             _1, _2, _3, vertex_label2, _4, edge_label2,
                             ref(F), max_n),
             directed_graph_view_pointers(), vertex_props_t(),
             edge_props_t())
            (gi2.GetGraphView(), vertex_label1,  edge_label1);
    }
    else
    {
        run_action<graph_tool::detail::never_directed>()
            (gi1, bind<void>(get_subgraphs(),
                             _1, _2, _3, vertex_label2, _4, edge_label2,
                             ref(F), max_n),
             undirected_graph_view_pointers(), vertex_props_t(),
             edge_props_t())
            (gi2.GetGraphView(), vertex_label1, edge_label1);
    }

    for (size_t i = 0; i < F.size(); ++i)
    {
        typedef property_map_type
            ::apply<int64_t, GraphInterface::vertex_index_map_t>::type
            vmap_t;
        typedef property_map_type
            ::apply<int64_t, GraphInterface::edge_index_map_t>::type
            emap_t;
        vmap_t::unchecked_t vm(gi1.GetVertexIndex(), gi1.GetNumberOfVertices());
        emap_t::unchecked_t ep(gi1.GetEdgeIndex(), gi1.GetMaxEdgeIndex()+1);
        if (gi1.GetDirected())
        {
            run_action<graph_tool::detail::always_directed>()
                (gi1, bind<void>(get_mapping(),
                                 _1, _2, _3, edge_label2,
                                 ref(F[i]), ref(vm), ref(ep),
                                 gi2.GetEdgeIndex()),
                 directed_graph_view_pointers(), edge_props_t())
                (gi2.GetGraphView(), edge_label1);
        }
        else
        {
            run_action<graph_tool::detail::never_directed>()
                (gi1, bind<void>(get_mapping(),
                                 _1, _2, _3, edge_label2,
                                 ref(F[i]), ref(vm), ref(ep),
                                 gi2.GetEdgeIndex()),
                 undirected_graph_view_pointers(), edge_props_t())
                (gi2.GetGraphView(), edge_label1);
        }
        vmapping.append(PythonPropertyMap<vmap_t>(vm.get_checked()));
        emapping.append(PythonPropertyMap<emap_t>(ep.get_checked()));
    }
}
