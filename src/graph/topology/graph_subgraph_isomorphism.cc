// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2013 Tiago de Paula Peixoto <tiago@skewed.de>
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

#include "random.hh"

#include <graph_subgraph_isomorphism.hh>
#include <graph_python_interface.hh>

#include "tr1_include.hh"
#include TR1_HEADER(unordered_map)

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
                    vector<size_t>& vlist, pair<reference_wrapper<rng_t>,size_t> sn) const
    {
        typedef PropLabelling<Graph1,Graph2,VertexLabel,VertexLabel>
            vlabelling_t;
        typedef PropLabelling<Graph1,Graph2,EdgeLabel,EdgeLabel>
            elabelling_t;
        rng_t& rng = sn.first;
        size_t max_n = sn.second;

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
        #pragma omp parallel for default(shared) private(i) schedule(static) if (N > 100)
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

struct vertex_label_mapping
{
    template <class Graph, class VertexMap, class VertexLabel>
    void operator()(const Graph& g, VertexMap vmap, VertexLabel vlabel,
                    boost::any& vlist) const
    {
        typedef typename property_traits<VertexMap>::value_type value_t;
        tr1::unordered_map<value_t, int32_t> values;
        if (!vlist.empty())
            values = boost::any_cast<tr1::unordered_map<value_t, int32_t> >(vlist);

        typename graph_traits<Graph>::vertex_iterator v, v_end;
        for (tie(v, v_end) = vertices(g); v != v_end; ++v)
        {
            const value_t& val = vmap[*v];
            if (values.find(val) == values.end())
                values[val] = values.size();
            vlabel[*v] = values[val];
        }
        vlist = values;
    }
};

struct edge_label_mapping
{
    template <class Graph, class EdgeMap, class EdgeLabel>
    void operator()(const Graph& g, EdgeMap emap, EdgeLabel elabel,
                    boost::any& vlist) const
    {
        typedef typename property_traits<EdgeMap>::value_type value_t;
        tr1::unordered_map<value_t, int32_t> values;
        if (!vlist.empty())
            values = boost::any_cast<tr1::unordered_map<value_t, int32_t> >(vlist);

        typename graph_traits<Graph>::edge_iterator e, e_end;
        for (tie(e, e_end) = edges(g); e != e_end; ++e)
        {
            const value_t& val = emap[*e];
            if (values.find(val) == values.end())
                values[val] = values.size();
            elabel[*e] = values[val];
        }
        vlist = values;
    }
};

// typedef mpl::push_back<vertex_properties,
//                        ConstantPropertyMap<bool,GraphInterface::vertex_t> >
//     ::type vertex_props_t;

// typedef mpl::push_back<edge_properties,
//                        ConstantPropertyMap<bool,GraphInterface::edge_t> >
//     ::type edge_props_t;

typedef typename property_map_type::apply<int32_t, GraphInterface::vertex_index_map_t>::type vlabel_t;
typedef mpl::vector<vlabel_t, ConstantPropertyMap<bool,GraphInterface::vertex_t> > vertex_props_t;

typedef typename property_map_type::apply<int32_t, GraphInterface::edge_index_map_t>::type elabel_t;
typedef mpl::vector<elabel_t, ConstantPropertyMap<bool,GraphInterface::edge_t> > edge_props_t;

void subgraph_isomorphism(GraphInterface& gi1, GraphInterface& gi2,
                          boost::any vertex_label1, boost::any vertex_label2,
                          boost::any edge_label1, boost::any edge_label2,
                          python::list vmapping, python::list emapping,
                          size_t max_n, rng_t& rng)
{
    if (gi1.GetDirected() != gi2.GetDirected())
        return;

    if (belongs<vertex_props_t>()(vertex_label1))
    {
        vertex_label2 = any_cast<vlabel_t>(vertex_label2).get_unchecked(num_vertices(gi2.GetGraph()));
    }
    else
    {
        if (vertex_label1.empty())
        {
            vertex_label1 = vertex_label2 =
                ConstantPropertyMap<bool,GraphInterface::vertex_t>(true);
        }
        else
        {
            boost::any vlist;
            vlabel_t vlabel1(gi1.GetVertexIndex());
            run_action<>()(gi1, bind<void>(vertex_label_mapping(),
                                           _1, _2, vlabel1, ref(vlist)),
                           vertex_properties())
                (vertex_label1);
            vertex_label1 = vlabel1;

            vlabel_t vlabel2(gi2.GetVertexIndex());
            run_action<>()(gi2, bind<void>(vertex_label_mapping(),
                                           _1, _2, vlabel2, ref(vlist)),
                           vertex_properties())
                (vertex_label2);
            vertex_label2 = vlabel2.get_unchecked(num_vertices(gi2.GetGraph()));
        }
    }

    if (belongs<edge_props_t>()(edge_label1))
    {
        edge_label2 = any_cast<elabel_t>(edge_label2).get_unchecked(gi2.GetGraph().get_last_index());
    }
    else
    {
        if (edge_label1.empty())
        {
            edge_label1 = edge_label2 =
                ConstantPropertyMap<bool,GraphInterface::edge_t>(true);
        }
        else
        {
            boost::any vlist;
            elabel_t elabel1(gi1.GetEdgeIndex());
            run_action<>()(gi1, bind<void>(edge_label_mapping(),
                                           _1, _2, elabel1, ref(vlist)),
                           edge_properties())
                (edge_label1);
            edge_label1 = elabel1;

            elabel_t elabel2(gi2.GetEdgeIndex());
            run_action<>()(gi2, bind<void>(edge_label_mapping(),
                                           _1, _2, elabel2, ref(vlist)),
                           edge_properties())
                (edge_label2);
            edge_label2 = elabel2.get_unchecked(num_edges(gi2.GetGraph()));
        }
    }

    vector<vector<pair<size_t,size_t> > > F;
    vector<size_t> vlist;

    if (gi1.GetDirected())
    {
        run_action<graph_tool::detail::always_directed>()
            (gi1, bind<void>(get_subgraphs(),
                             _1, _2, _3, vertex_label2, _4, edge_label2,
                             ref(F), ref(vlist), make_pair(ref(rng), max_n)),
             directed_graph_view_pointers(), vertex_props_t(),
             edge_props_t())
            (gi2.GetGraphView(), vertex_label1,  edge_label1);
    }
    else
    {
        run_action<graph_tool::detail::never_directed>()
            (gi1, bind<void>(get_subgraphs(),
                             _1, _2, _3, vertex_label2, _4, edge_label2,
                             ref(F), ref(vlist), make_pair(ref(rng), max_n)),
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
