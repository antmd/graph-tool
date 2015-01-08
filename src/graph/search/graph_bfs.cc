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
#include "graph_python_interface.hh"

#include <boost/python.hpp>
#include <boost/graph/breadth_first_search.hpp>

#include "graph.hh"
#include "graph_selectors.hh"
#include "graph_util.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

class BFSVisitorWrapper
{
public:
    BFSVisitorWrapper(python::object gi, python::object vis)
        : _gi(gi), _vis(vis) {}

    template <class Vertex, class Graph>
    void initialize_vertex(Vertex u, const Graph& g)
    {
        _vis.attr("initialize_vertex")(PythonVertex(_gi, u));
    }

    template <class Vertex, class Graph>
    void discover_vertex(Vertex u, const Graph& g)
    {
        _vis.attr("discover_vertex")(PythonVertex(_gi, u));
    }

    template <class Vertex, class Graph>
    void examine_vertex(Vertex u, const Graph& g)
    {
        _vis.attr("examine_vertex")(PythonVertex(_gi, u));
    }

    template <class Edge, class Graph>
    void examine_edge(Edge e, const Graph& g)
    {
        _vis.attr("examine_edge")
            (PythonEdge<Graph>(_gi, e));
    }

    template <class Edge, class Graph>
    void tree_edge(Edge e, const Graph& g)
    {
        _vis.attr("tree_edge")
            (PythonEdge<Graph>(_gi, e));
    }

    template <class Edge, class Graph>
    void non_tree_edge(Edge e, const Graph& g)
    {
        _vis.attr("non_tree_edge")
            (PythonEdge<Graph>(_gi, e));
    }

    template <class Edge, class Graph>
    void gray_target(Edge e, const Graph& g)
    {
        _vis.attr("gray_target")
            (PythonEdge<Graph>(_gi, e));
    }

    template <class Edge, class Graph>
    void black_target(Edge e, const Graph& g)
    {
        _vis.attr("black_target")
            (PythonEdge<Graph>(_gi, e));
    }

    template <class Vertex, class Graph>
    void finish_vertex(Vertex u, const Graph& g)
    {
        _vis.attr("finish_vertex")(PythonVertex(_gi, u));
    }

private:
    python::object _gi, _vis;
};

struct do_bfs
{
    template <class Graph>
    void operator()(Graph& g, size_t s, BFSVisitorWrapper vis) const
    {
        breadth_first_search(g, vertex(s, g), visitor(vis));
    }
};

void bfs_search(GraphInterface& g, python::object gi, size_t s,
                python::object vis)
{
    run_action<graph_tool::detail::all_graph_views,mpl::true_>()
        (g, std::bind(do_bfs(), placeholders::_1, s,
                      BFSVisitorWrapper(gi, vis)))();
}

void export_bfs()
{
    using namespace boost::python;
    def("bfs_search", &bfs_search);
}
