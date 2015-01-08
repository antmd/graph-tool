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
#include <cxxabi.h>

using namespace graph_tool;
using namespace graph_tool::detail;
using namespace boost;

bool graph_tool::graph_filtering_enabled()
{
#ifndef NO_GRAPH_FILTERING
    return true;
#else
    return false;
#endif
}

string name_demangle(string name)
{
    int status = 0;
    char *realname = abi::__cxa_demangle(name.c_str(), 0, 0, &status);
    string ret(realname);
    free(realname);
    return ret;
}

// Whenever no implementation is called, the following exception is thrown
graph_tool::ActionNotFound::ActionNotFound(const boost::any& graph_view,
                                           const type_info& action,
                                           const vector<const type_info*>& args)
    : GraphException(""), _graph_view(graph_view),
      _action(action), _args(args) {}

const char * graph_tool::ActionNotFound::what () const throw ()
{
    using python::detail::gcc_demangle;

    string error =
        "No static implementation was found for the desired routine. "
        "This is a graph_tool bug. :-( Please follow bug report "
        "instructions at " PACKAGE_BUGREPORT ". What follows is debug "
        "information.\n\n";

    error += "Graph view: " + name_demangle(_graph_view.type().name()) + "\n\n";
    error += "Action: " + name_demangle(_action.name()) + "\n\n";
    for (size_t i = 0; i < _args.size(); ++i)
    {
        error += "Arg " + lexical_cast<string>(i+1) + ": " +
            name_demangle(_args[i]->name()) + "\n\n";
    }
    return error.c_str();
}

// this function retrieves a graph view stored in graph_views, or stores one if
// non-existent
template <class Graph>
typename std::remove_const<Graph>::type&
retrieve_graph(vector<boost::any>& graph_views, Graph& init)
{
    typedef typename std::remove_const<Graph>::type g_t;
    size_t index = boost::mpl::find<all_graph_views,g_t>::type::pos::value;
    if (index >= graph_views.size())
        graph_views.resize(index+1);
    boost::any gview = graph_views[index];
    std::shared_ptr<g_t>* gptr = any_cast<std::shared_ptr<g_t> >(&gview);
    if (gptr == 0)
    {
        std::shared_ptr<g_t> new_g(new g_t(init));
        gptr = &new_g;
        gview = new_g;
        graph_views[index] = gview;
    }
    return **gptr;
}


// this will check whether a graph is reversed and return the proper view
// encapsulated
template <class Graph>
boost::any check_reverse(const Graph &g, bool reverse,
                         vector<boost::any>& graph_views)
{
    if (reverse)
    {
        typedef typename boost::mpl::if_<std::is_const<Graph>,
                                         const reverse_graph<typename std::remove_const<Graph>::type>,
                                         reverse_graph<Graph> >::type
            reverse_graph_t;

        reverse_graph_t rg(g);
        return &retrieve_graph(graph_views, rg);
    }

    return boost::any(&const_cast<Graph&>(g));
};

// this will check whether a graph is directed and return the proper view
// encapsulated
template <class Graph>
boost::any check_directed(const Graph &g, bool reverse, bool directed,
                          vector<boost::any>& graph_views)
{
    if (directed)
    {
        return check_reverse(g, reverse, graph_views);
    }

    typedef UndirectedAdaptor<Graph> ug_t;
    ug_t ug(g);
    return &retrieve_graph(graph_views, ug);
};

// this will check whether a graph is filtered and return the proper view
// encapsulated
template <class Graph, class EdgeFilter, class VertexFilter>
boost::any
check_filtered(const Graph &g, const EdgeFilter& edge_filter,
               const bool& e_invert, bool e_active, size_t max_eindex,
               const VertexFilter& vertex_filter, const bool& v_invert,
               bool v_active, vector<boost::any>& graph_views, bool reverse,
               bool directed)
{
#ifndef NO_GRAPH_FILTERING
    MaskFilter<EdgeFilter> e_filter(edge_filter, e_invert);
    MaskFilter<VertexFilter> v_filter(vertex_filter, v_invert);

    if (e_active)
    {
        if (!v_active)
            throw GraphException("Edge filter is active but vertex filter is not. This is a bug.");

        if (max_eindex > 0)
            edge_filter.reserve(max_eindex);
        if (num_vertices(g) > 0)
            vertex_filter.reserve(num_vertices(g));
        typedef filtered_graph<Graph, MaskFilter<EdgeFilter>,
                               MaskFilter<VertexFilter> > fg_t;
        fg_t init(g, e_filter, v_filter);
        fg_t& fg = retrieve_graph(graph_views, init);
        fg.m_edge_pred = e_filter;
        fg.m_vertex_pred = v_filter;
        return check_directed(fg, reverse, directed, graph_views);
    }
    else
    {
        if (v_active)
            throw GraphException("Vertex filter is active but edge filter is not. This is a bug.");

        return check_directed(g, reverse, directed, graph_views);
    }
#else
    return check_directed(g, reverse, directed, graph_views);
#endif
}

// gets the correct graph view at run time
boost::any GraphInterface::GetGraphView() const
{
    boost::any graph =
        check_filtered(*_mg, _edge_filter_map, _edge_filter_invert,
                       _edge_filter_active, _mg->get_last_index(),
                       _vertex_filter_map, _vertex_filter_invert,
                       _vertex_filter_active,
                       const_cast<vector<boost::any>&>(_graph_views), _reversed,
                       _directed);
    return graph;
}

// these test whether or not the vertex and edge filters are active
bool GraphInterface::IsVertexFilterActive() const
{ return _vertex_filter_active; }

bool GraphInterface::IsEdgeFilterActive() const
{ return _edge_filter_active; }


// this function will reindex all the edges, in the order in which they are
// found
void GraphInterface::ReIndexEdges()
{
    _mg->reindex_edges();
}

// this will definitively remove all the edges from the graph, which are being
// currently filtered out. This will also disable the edge filter
void GraphInterface::PurgeEdges()
{
    if (!IsEdgeFilterActive())
        return;

    MaskFilter<edge_filter_t> filter(_edge_filter_map, _edge_filter_invert);
    graph_traits<multigraph_t>::vertex_iterator v, v_end;
    graph_traits<multigraph_t>::out_edge_iterator e, e_end;
    vector<graph_traits<multigraph_t>::edge_descriptor> deleted_edges;
    for (tie(v, v_end) = vertices(*_mg); v != v_end; ++v)
    {
        for (tie(e, e_end) = out_edges(*v, *_mg); e != e_end; ++e)
            if (!filter(*e))
                deleted_edges.push_back(*e);
        for (typeof(deleted_edges.begin()) iter = deleted_edges.begin();
             iter != deleted_edges.end(); ++iter)
            remove_edge(*iter, *_mg);
        deleted_edges.clear();
    }
}


// this will definitively remove all the vertices from the graph, which are
// being currently filtered out. This will also disable the vertex filter
void GraphInterface::PurgeVertices(boost::any aold_index)
{
    if (!IsVertexFilterActive())
        return;

    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        index_prop_t;
    index_prop_t old_index = any_cast<index_prop_t>(aold_index);

    MaskFilter<vertex_filter_t> filter(_vertex_filter_map,
                                       _vertex_filter_invert);
    size_t N = num_vertices(*_mg);
    vector<bool> deleted(N, false);
    for (size_t i = 0; i < N; ++i)
        deleted[i] = !filter(vertex(i, *_mg));
    vector<int> old_indexes;

    vector<graph_traits<multigraph_t>::edge_descriptor> edges;

    //remove vertices
    for (int i = N-1; i >= 0; --i)
    {
        if (deleted[i])
        {
            graph_traits<multigraph_t>::vertex_descriptor v =
                vertex(i, *_mg);
            remove_vertex(v, *_mg);
        }
        else
        {
            old_indexes.push_back(i);
        }
    }

    N = old_indexes.size();
    for (int i = N-1; i >= 0; --i)
    {
        old_index[vertex((N - 1) - i, *_mg)] = old_indexes[i];
    }
}

void GraphInterface::SetVertexFilterProperty(boost::any property, bool invert)
{
#ifdef NO_GRAPH_FILTERING
    throw GraphException("graph filtering was not enabled at compile time");
#endif

    try
    {
        _vertex_filter_map =
            any_cast<vertex_filter_t::checked_t>(property).get_unchecked();
        _vertex_filter_invert = invert;
        _vertex_filter_active = true;
    }
    catch(bad_any_cast&)
    {
        if (!property.empty())
            throw GraphException("Invalid vertex filter property!");
        _vertex_filter_active = false;
    }
}

void GraphInterface::SetEdgeFilterProperty(boost::any property, bool invert)
{
#ifdef NO_GRAPH_FILTERING
    throw GraphException("graph filtering was not enabled at compile time");
#endif

    try
    {
        _edge_filter_map =
            any_cast<edge_filter_t::checked_t>(property).get_unchecked();
        _edge_filter_invert = invert;
        _edge_filter_active = true;
    }
    catch(bad_any_cast&)
    {
        if (!property.empty())
            throw GraphException("Invalid edge filter property!");
        _edge_filter_active = false;
    }
}
