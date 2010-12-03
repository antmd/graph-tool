// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "graph_filtering.hh"
#include "graph_python_interface.hh"

#include <boost/python.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/undirected_dfs.hpp>

#include "graph.hh"
#include "graph_selectors.hh"
#include "graph_util.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;


class DFSVisitorWrapper
{
public:
    DFSVisitorWrapper(python::object& gi, python::object vis)
        : _gi(gi), _vis(vis) {}


    template <class Vertex, class Graph>
    void initialize_vertex(Vertex u, const Graph& g)
    {
        _vis.attr("initialize_vertex")(PythonVertex(_gi, u));
    }
    template <class Vertex, class Graph>
    void start_vertex(Vertex u, const Graph& g)
    {
        _vis.attr("start_vertex")(PythonVertex(_gi, u));
    }
    template <class Vertex, class Graph>
    void discover_vertex(Vertex u, const Graph& g)
    {
        _vis.attr("discover_vertex")(PythonVertex(_gi, u));
    }

    template <class Edge, class Graph>
    void examine_edge(Edge e, const Graph& g)
    {
        _vis.attr("examine_edge")
            (PythonEdge<typename Graph::orig_graph_t>(_gi, e));
    }

    template <class Edge, class Graph>
    void tree_edge(Edge e, const Graph& g)
    {
        _vis.attr("tree_edge")
            (PythonEdge<typename Graph::orig_graph_t>(_gi, e));
    }

    template <class Edge, class Graph>
    void back_edge(Edge e, const Graph& g)
    {
        _vis.attr("back_edge")
            (PythonEdge<typename Graph::orig_graph_t>(_gi, e));
    }

    template <class Edge, class Graph>
    void forward_or_cross_edge(Edge e, const Graph& g)
    {
        _vis.attr("forward_or_cross_edge")
            (PythonEdge<typename Graph::orig_graph_t>(_gi, e));
    }

    template <class Vertex, class Graph>
    void finish_vertex(Vertex u, const Graph& g)
    {
        _vis.attr("finish_vertex")(PythonVertex(_gi, u));
    }

private:
    python::object _gi, _vis;
};

struct do_dfs
{
    template <class Graph, class EdgeIndexMap>
    void operator()(Graph& g, EdgeIndexMap edge_index, size_t num_e, size_t s,
                    DFSVisitorWrapper vis) const
    {
        dfs_dispatch(g, edge_index, num_e, s, vis,
                     typename is_directed::apply<Graph>::type());
    }

    template <class Graph, class EdgeIndexMap>
    void dfs_dispatch(Graph& g, EdgeIndexMap edge_index, size_t num_e, size_t s,
                      DFSVisitorWrapper vis, mpl::true_ is_directed) const
    {
        depth_first_search(g, visitor(vis).root_vertex(vertex(s, g)));
    }

    template <class Graph, class EdgeIndexMap>
    void dfs_dispatch(Graph& g, EdgeIndexMap edge_index, size_t num_e, size_t s,
                      DFSVisitorWrapper vis, mpl::false_ is_directed) const
    {
        typename property_map_type::apply<default_color_type,
                                          EdgeIndexMap>::type
            ecolor(edge_index);
        undirected_dfs(g, visitor(vis).edge_color_map(ecolor).
                       root_vertex(vertex(s, g)));
    }
};


void dfs_search(GraphInterface& g, python::object gi, size_t s,
                python::object vis)
{
    run_action<graph_tool::detail::all_graph_views,mpl::true_>()
        (g, bind<void>(do_dfs(), _1, g.GetEdgeIndex(),
                       g.GetMaxEdgeIndex() + 1, s,
                       DFSVisitorWrapper(gi, vis)))();
}

void export_dfs()
{
    using namespace boost::python;
    def("dfs_search", &dfs_search);
}
