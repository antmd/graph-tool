// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007-2011 Tiago de Paula Peixoto <tiago@skewed.de>
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
                    typename graph_traits<Graph2>::vertex_descriptor v2) const
    {
        if (in_degreeS()(v2, _g2) < in_degreeS()(v1, _g1) ||
            out_degree(v2, _g2) < out_degree(v1, _g1))
            return false;
        return _label1[v1] == _label2[v2];
    }

    bool operator()(typename graph_traits<Graph1>::edge_descriptor e1,
                    typename graph_traits<Graph2>::edge_descriptor e2) const
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
    void operator()(const Graph1& sub, const Graph2* g,
                    VertexLabel vertex_label1, boost::any vertex_label2,
                    EdgeLabel edge_label1, boost::any edge_label2,
                    vector<vector<pair<size_t, size_t> > >& F,
                    vector<size_t>& vlist, pair<size_t,size_t> sn) const
    {
        typedef PropLabelling<Graph1,Graph2,VertexLabel,VertexLabel>
            vlabelling_t;
        typedef PropLabelling<Graph1,Graph2,EdgeLabel,EdgeLabel>
            elabelling_t;
        size_t seed = sn.first;
        size_t max_n = sn.second;

        rng_t rng(static_cast<rng_t::result_type>(seed));
        vlist.resize(num_vertices(*g));
        int i, N = num_vertices(*g);
        for (i = 0; i < N; ++i)
            vlist[i] = i;
        for (i = 0; i < N - 1; ++i)
        {
            tr1::uniform_int<> random(i, N - 1);
            swap(vlist[i], vlist[random(rng)]);
        }

        subgraph_isomorphism
            (sub, *g, vlabelling_t(sub, *g, vertex_label1,
                                   any_cast<VertexLabel>(vertex_label2)),
             elabelling_t(sub, *g, edge_label1,
                          any_cast<EdgeLabel>(edge_label2)), F, vlist, max_n);
    }
};

struct get_mapping
{
    template <class Graph1, class Graph2, class EdgeLabel, class VertexMap,
              class EdgeMap, class EdgeIndexMap>
    void operator()(const Graph1& sub, const Graph2* g, EdgeLabel edge_label1,
                    boost::any edge_label2, vector<pair<size_t, size_t> >& F,
                    VertexMap vmapping, EdgeMap emapping,
                    EdgeIndexMap edge_index2, vector<size_t>& vlist) const
    {
        typedef typename graph_traits<Graph1>::edge_descriptor edge1_t;
        typedef typename graph_traits<Graph2>::edge_descriptor edge2_t;
        typedef PropLabelling<Graph1,Graph2,EdgeLabel,EdgeLabel>
            elabelling_t;
        elabelling_t edge_labelling(sub, *g, edge_label1,
                                    any_cast<EdgeLabel>(edge_label2));
        int i, N = F.size();
        #pragma omp parallel for default(shared) private(i) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            if (vertex(i, sub) == graph_traits<Graph1>::null_vertex())
                continue;
            vmapping[vertex(F[i].first, sub)] = vertex(vlist[F[i].second], *g);
            typename graph_traits<Graph1>::out_edge_iterator e, e_end;
            for (tie(e, e_end) = out_edges(vertex(i, sub), sub); e != e_end;
                 ++e)
            {
                bool found = false;
                typename graph_traits<Graph2>::out_edge_iterator e2, e2_end;
                for (tie(e2, e2_end) =
                         out_edges(vertex(vlist[F[i].second], *g), *g);
                     e2 != e2_end; ++e2)
                {
                    if (target(edge2_t(*e2), *g) ==
                        vertex(vlist[F[target(*e, sub)].second], *g) &&
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
                          size_t max_n, size_t seed)
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
    vector<size_t> vlist;

    if (gi1.GetDirected())
    {
        run_action<graph_tool::detail::always_directed>()
            (gi1, bind<void>(get_subgraphs(),
                             _1, _2, _3, vertex_label2, _4, edge_label2,
                             ref(F), ref(vlist), make_pair(seed, max_n)),
             directed_graph_view_pointers(), vertex_props_t(),
             edge_props_t())
            (gi2.GetGraphView(), vertex_label1,  edge_label1);
    }
    else
    {
        run_action<graph_tool::detail::never_directed>()
            (gi1, bind<void>(get_subgraphs(),
                             _1, _2, _3, vertex_label2, _4, edge_label2,
                             ref(F), ref(vlist), make_pair(seed, max_n)),
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
                                 gi2.GetEdgeIndex(), ref(vlist)),
                 directed_graph_view_pointers(), edge_props_t())
                (gi2.GetGraphView(), edge_label1);
        }
        else
        {
            run_action<graph_tool::detail::never_directed>()
                (gi1, bind<void>(get_mapping(),
                                 _1, _2, _3, edge_label2,
                                 ref(F[i]), ref(vm), ref(ep),
                                 gi2.GetEdgeIndex(), ref(vlist)),
                 undirected_graph_view_pointers(), edge_props_t())
                (gi2.GetGraphView(), edge_label1);
        }
        vmapping.append(PythonPropertyMap<vmap_t>(vm.get_checked()));
        emapping.append(PythonPropertyMap<emap_t>(ep.get_checked()));
    }
}
