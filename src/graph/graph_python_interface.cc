// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2015 Tiago de Paula Peixoto <tiago@skewed.de>
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

#include "graph_filtering.hh"
#include "graph.hh"
#include "graph_python_interface.hh"

#include <boost/python.hpp>
#include <set>

using namespace std;
using namespace boost;
using namespace graph_tool;

namespace boost
{
size_t hash_value(const boost::python::object& o)
{
    return std::hash<boost::python::object>()(o);
}
}

namespace graph_tool
{

struct get_vertex_iterator
{
    template <class Graph>
    void operator()(Graph& g, python::object& pg,
                    python::object& iter) const
    {
        typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
        iter = python::object(PythonIterator<PythonVertex,
                                             vertex_iterator>(pg, vertices(g)));
    }
};

python::object get_vertices(python::object g)
{
    GraphInterface& gi = python::extract<GraphInterface&>(g().attr("_Graph__graph"));
    python::object iter;
    run_action<>()(gi, std::bind(get_vertex_iterator(), placeholders::_1,
                                 std::ref(g),
                                 std::ref(iter)))();
    return iter;
}

struct get_vertex_soft
{
    template <class Graph>
    void operator()(Graph& g,
                    python::object& pg,
                    size_t i, python::object& v) const
    {
        v = python::object(PythonVertex(pg, vertex(i, g)));
    }
};

struct get_vertex_hard
{
    template <class Graph>
    void operator()(Graph& g, python::object& pg, size_t i,
                    python::object& v) const
    {
        size_t c = 0;
        typename graph_traits<Graph>::vertex_iterator vi, v_end;
        for (tie(vi, v_end) = vertices(g); vi != v_end; ++vi)
        {
            if (c == i)
            {
                v = python::object(PythonVertex(pg, *vi));
                return;
            }
            ++c;
        }
        v = python::object(PythonVertex(pg,
                                        graph_traits<Graph>::null_vertex()));
    }
};

python::object get_vertex(python::object g, size_t i)
{
    GraphInterface& gi = python::extract<GraphInterface&>(g().attr("_Graph__graph"));
    python::object v;
    if (gi.IsVertexFilterActive())
        run_action<>()(gi,
                       std::bind(get_vertex_hard(), placeholders::_1,
                                 std::ref(g), i, std::ref(v)))();
    else
        run_action<>()(gi,
                       std::bind(get_vertex_soft(), placeholders::_1,
                                 std::ref(g), i, std::ref(v)))();
    return v;
}

struct get_edge_iterator
{
    template <class Graph>
    void operator()(Graph& g, const python::object& pg, python::object& iter)
        const
    {
        typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
        iter = python::object(PythonIterator<PythonEdge<Graph>,
                                             edge_iterator>(pg, edges(g)));
    }
};

python::object get_edges(python::object g)
{
    GraphInterface& gi = python::extract<GraphInterface&>(g().attr("_Graph__graph"));
    python::object iter;
    run_action<>()(gi,
                   std::bind(get_edge_iterator(), placeholders::_1,
                             std::ref(g), std::ref(iter)))();
    return iter;
}

python::object add_vertex(python::object g, size_t n)
{
    GraphInterface& gi = python::extract<GraphInterface&>(g().attr("_Graph__graph"));

    if (n > 1)
    {
        for (size_t i = 0; i < n; ++i)
            add_vertex(gi.GetGraph());
        return python::object();
    }
    return python::object(PythonVertex(g, add_vertex(gi.GetGraph())));
}

struct shift_vertex_property
{
    template <class Graph, class PropertyMap>
    void operator()(const Graph& g, size_t vi, PropertyMap pmap)
        const
    {
        size_t N = num_vertices(g);
        if (N <= 1 || vi >= N - 1)
            return;
        for (size_t i = vi; i < N-1; ++i)
        {
            pmap[vertex(i, g)] = pmap[vertex(i+1, g)];
        }
    }
};

void remove_vertex(GraphInterface& gi, const python::object& v, bool fast)
{
    PythonVertex& pv = python::extract<PythonVertex&>(v);
    pv.CheckValid();
    GraphInterface::vertex_t dv = pv.GetDescriptor();
    pv.SetValid(false);

    if (fast)
        remove_vertex_fast(dv, gi.GetGraph());
    else
        remove_vertex(dv, gi.GetGraph());
}

struct add_new_edge
{
    template <class Graph, class EdgeIndexMap>
    void operator()(Graph& g, python::object& pg, GraphInterface&,
                    const PythonVertex& s, const PythonVertex& t,
                    EdgeIndexMap, python::object& new_e) const
    {
        typename graph_traits<Graph>::edge_descriptor e =
            add_edge(s.GetDescriptor(), t.GetDescriptor(), g).first;
        new_e = python::object(PythonEdge<Graph>(pg, e));
    }
};

python::object add_edge(python::object g, const python::object& s,
                        const python::object& t)
{
    PythonVertex& src = python::extract<PythonVertex&>(s);
    PythonVertex& tgt = python::extract<PythonVertex&>(t);
    src.CheckValid();
    tgt.CheckValid();
    GraphInterface& gi = python::extract<GraphInterface&>(g().attr("_Graph__graph"));
    python::object new_e;
    run_action<>()(gi, std::bind(add_new_edge(), placeholders::_1, std::ref(g),
                                 std::ref(gi), src, tgt, gi.GetEdgeIndex(),
                                 std::ref(new_e)))();
    return new_e;
}

struct get_edge_descriptor
{
    template <class Graph>
    void operator()(const Graph& g, const python::object& e,
                    typename GraphInterface::edge_t& edge,
                    bool& found)  const
    {
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;
        PythonEdge<Graph>& pe = python::extract<PythonEdge<Graph>&>(e);
        pe.CheckValid();
        pe.SetValid(false);
        edge = pe.GetDescriptor();
        found = true;
    }
};

void remove_edge(GraphInterface& gi, const python::object& e)
{
    GraphInterface::edge_t de;
    bool found = false;
    run_action<>()(gi, std::bind(get_edge_descriptor(), placeholders::_1,
                                 std::ref(e), std::ref(de), std::ref(found)))();
    remove_edge(de, gi.GetGraph());
    if (!found)
        throw ValueException("invalid edge descriptor");
}

struct get_degree_map
{
    template <class Graph, class DegS, class Weight>
    void operator()(const Graph& g, python::object& odeg_map, DegS deg, Weight weight) const
    {
        typedef typename detail::get_weight_type<Weight>::type weight_t;
        typedef typename mpl::if_<std::is_same<weight_t, size_t>, int32_t, weight_t>::type deg_t;

        typedef typename property_map_type::apply<deg_t,
                                                  GraphInterface::vertex_index_map_t>::type
            map_t;

        map_t cdeg_map(get(vertex_index, g));
        typename map_t::unchecked_t deg_map = cdeg_map.get_unchecked(num_vertices(g));

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) schedule(runtime) if (N > 100)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            deg_map[v] = deg(v, g, weight);
        }

        odeg_map = python::object(PythonPropertyMap<map_t>(cdeg_map));
    }
};

python::object GraphInterface::DegreeMap(string deg, boost::any weight) const
{

    python::object deg_map;

    typedef mpl::push_back<edge_scalar_properties,
                           detail::no_weightS>::type weight_t;
    if (weight.empty())
        weight = detail::no_weightS();

    if (deg == "in")
        run_action<>()(const_cast<GraphInterface&>(*this),
                       std::bind(get_degree_map(), placeholders::_1,
                                 std::ref(deg_map), in_degreeS(), placeholders::_2), weight_t())
            (weight);
    else if (deg == "out")
        run_action<>()(const_cast<GraphInterface&>(*this),
                       std::bind(get_degree_map(), placeholders::_1,
                                 std::ref(deg_map), out_degreeS(), placeholders::_2), weight_t())
            (weight);
    else if (deg == "total")
        run_action<>()(const_cast<GraphInterface&>(*this),
                       std::bind(get_degree_map(), placeholders::_1,
                                 std::ref(deg_map), total_degreeS(), placeholders::_2), weight_t())
            (weight);
    return deg_map;
}

//
// Below are the functions with will properly register all the types to python,
// for every filter, type, etc.
//

// this will register all the Vertex/Edge classes to python
struct export_python_interface
{
    template <class Graph>
    void operator()(const Graph*, set<string>& v_iterators) const
    {
        using namespace boost::python;

        class_<PythonEdge<Graph>, bases<EdgeBase> >
            ("Edge", no_init)
            .def("source", &PythonEdge<Graph>::GetSource,
                 "Return the source vertex.")
            .def("target", &PythonEdge<Graph>::GetTarget,
                 "Return the target vertex.")
            .def("is_valid", &PythonEdge<Graph>::IsValid,
                 "Return whether the edge is valid.")
            .def("get_graph", &PythonEdge<Graph>::GetGraph,
                 "Return the graph to which the edge belongs.")
            .def("__str__", &PythonEdge<Graph>::GetString)
            .def("__hash__", &PythonEdge<Graph>::GetHash);

        typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
        if (v_iterators.find(typeid(vertex_iterator).name()) ==
            v_iterators.end())
        {
            class_<PythonIterator<PythonVertex, vertex_iterator> >
                ("VertexIterator", no_init)
                .def("__iter__", objects::identity_function())
                .def("__next__", &PythonIterator<PythonVertex,
                                                 vertex_iterator>::Next)
                .def("next", &PythonIterator<PythonVertex,
                                             vertex_iterator>::Next);
            v_iterators.insert(typeid(vertex_iterator).name());
        }

        typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
        class_<PythonIterator<PythonEdge<Graph>,
                              edge_iterator> >("EdgeIterator", no_init)
            .def("__iter__", objects::identity_function())
            .def("__next__", &PythonIterator<PythonEdge<Graph>,
                                         edge_iterator>::Next)
            .def("next", &PythonIterator<PythonEdge<Graph>,
                                         edge_iterator>::Next);

        typedef typename graph_traits<Graph>::out_edge_iterator
            out_edge_iterator;
        class_<PythonIterator<PythonEdge<Graph>,
                              out_edge_iterator> >("OutEdgeIterator", no_init)
            .def("__iter__", objects::identity_function())
            .def("__next__", &PythonIterator<PythonEdge<Graph>,
                                             out_edge_iterator>::Next)
            .def("next", &PythonIterator<PythonEdge<Graph>,
                                         out_edge_iterator>::Next);

        typedef typename graph_traits<Graph>::directed_category
            directed_category;
        typedef typename std::is_convertible<directed_category,
                                             boost::directed_tag>::type is_directed;
        if (is_directed::value)
        {
            typedef typename in_edge_iteratorS<Graph>::type in_edge_iterator;
            class_<PythonIterator<PythonEdge<Graph>,
                                  in_edge_iterator> >("InEdgeIterator", no_init)
                .def("__iter__", objects::identity_function())
                .def("__next__", &PythonIterator<PythonEdge<Graph>,
                                                 in_edge_iterator>::Next)
                .def("next", &PythonIterator<PythonEdge<Graph>,
                                             in_edge_iterator>::Next);
        }
    }
};

PythonPropertyMap<GraphInterface::vertex_index_map_t>
get_vertex_index(GraphInterface& g)
{
    return PythonPropertyMap<GraphInterface::vertex_index_map_t>
        (g.GetVertexIndex());
}

PythonPropertyMap<GraphInterface::edge_index_map_t>
do_get_edge_index(GraphInterface& g)
{
    return PythonPropertyMap<GraphInterface::edge_index_map_t>
        (g.GetEdgeIndex());
}


template <class ValueList>
struct add_edge_list
{
    template <class Graph>
    void operator()(Graph& g, python::object aedge_list, bool& found) const
    {
        boost::mpl::for_each<ValueList>(std::bind(dispatch(), std::ref(g),
                                                  std::ref(aedge_list),
                                                  std::ref(found), placeholders::_1));
    }

    struct dispatch
    {
        template <class Graph, class Value>
        void operator()(Graph& g, python::object& aedge_list, bool& found, Value) const
        {
            if (found)
                return;
            try
            {
                boost::multi_array_ref<Value, 2> edge_list = get_array<Value, 2>(aedge_list);

                if (edge_list.shape()[1] < 2)
                    throw GraphException("Second dimension in edge list must be of size (at least) two");

                for (size_t i = 0; i < edge_list.shape()[0]; ++i)
                {
                    size_t s = edge_list[i][0];
                    size_t t = edge_list[i][1];
                    while (s >= num_vertices(g) || t >= num_vertices(g))
                        add_vertex(g);
                    add_edge(vertex(s, g), vertex(t, g), g);
                }
                found = true;
            }
            catch (invalid_numpy_conversion& e) {}
        }
    };
};

void do_add_edge_list(GraphInterface& gi, python::object aedge_list)
{
    typedef mpl::vector<bool, uint8_t, uint32_t, int16_t, int32_t, int64_t, uint64_t,
                        unsigned long int, double, long double> vals_t;
    bool found = false;
    run_action<>()(gi, std::bind(add_edge_list<vals_t>(), placeholders::_1, aedge_list,
                                 std::ref(found)))();
    if (!found)
        throw GraphException("Invalid type for edge list; must be two-dimensional with a scalar type");
}


} // namespace graph_tool

// register everything

void export_python_properties();

void export_python_interface()
{
    using namespace boost::python;

    class_<PythonVertex>
        ("Vertex", no_init)
        .def("__in_degree", &PythonVertex::GetInDegree,
             "Return the in-degree.")
        .def("__weighted_in_degree", &PythonVertex::GetWeightedOutDegree,
             "Return the weighted in-degree.")
        .def("__out_degree", &PythonVertex::GetOutDegree,
             "Return the out-degree.")
        .def("__weighted_out_degree", &PythonVertex::GetWeightedOutDegree,
             "Return the weighted out-degree.")
        .def("in_edges", &PythonVertex::InEdges,
             "Return an iterator over the in-edges.")
        .def("out_edges", &PythonVertex::OutEdges,
             "Return an iterator over the out-edges.")
        .def("is_valid", &PythonVertex::IsValid,
             "Return whether the vertex is valid.")
        .def("get_graph", &PythonVertex::GetGraph,
             "Return the graph to which the vertex belongs.")
        .def("__str__", &PythonVertex::GetString)
        .def("__int__", &PythonVertex::GetIndex)
        .def("__hash__", &PythonVertex::GetHash);
    class_<EdgeBase>("EdgeBase", no_init);

    set<string> v_iterators;
    typedef boost::mpl::transform<graph_tool::detail::all_graph_views,
                                  boost::mpl::quote1<std::add_pointer> >::type graph_views;
    boost::mpl::for_each<graph_views>(std::bind(graph_tool::export_python_interface(),
                                      placeholders::_1, std::ref(v_iterators)));
    export_python_properties();
    def("new_vertex_property",
        &new_property<GraphInterface::vertex_index_map_t>);
    def("new_edge_property", &new_property<GraphInterface::edge_index_map_t>);
    def("new_graph_property",
        &new_property<ConstantPropertyMap<size_t,graph_property_tag> >);

    def("get_vertex", get_vertex);
    def("get_vertices", get_vertices);
    def("get_edges", get_edges);
    def("add_vertex", graph_tool::add_vertex);
    def("add_edge", graph_tool::add_edge);
    def("remove_vertex", graph_tool::remove_vertex);
    def("remove_edge", graph_tool::remove_edge);
    def("add_edge_list", graph_tool::do_add_edge_list);

    def("get_vertex_index", get_vertex_index);
    def("get_edge_index", do_get_edge_index);
}
