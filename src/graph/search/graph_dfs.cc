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
    void initialize_vertex(Vertex u, const Graph&)
    {
        _vis.attr("initialize_vertex")(PythonVertex(_gi, u));
    }
    template <class Vertex, class Graph>
    void start_vertex(Vertex u, const Graph&)
    {
        _vis.attr("start_vertex")(PythonVertex(_gi, u));
    }
    template <class Vertex, class Graph>
    void discover_vertex(Vertex u, const Graph&)
    {
        _vis.attr("discover_vertex")(PythonVertex(_gi, u));
    }

    template <class Edge, class Graph>
    void examine_edge(Edge e, const Graph&)
    {
        _vis.attr("examine_edge")
            (PythonEdge<Graph>(_gi, e));
    }

    template <class Edge, class Graph>
    void tree_edge(Edge e, const Graph&)
    {
        _vis.attr("tree_edge")
            (PythonEdge<Graph>(_gi, e));
    }

    template <class Edge, class Graph>
    void back_edge(Edge e, const Graph&)
    {
        _vis.attr("back_edge")
            (PythonEdge<Graph>(_gi, e));
    }

    template <class Edge, class Graph>
    void forward_or_cross_edge(Edge e, const Graph&)
    {
        _vis.attr("forward_or_cross_edge")
            (PythonEdge<Graph>(_gi, e));
    }

    template <class Vertex, class Graph>
    void finish_vertex(Vertex u, const Graph&)
    {
        _vis.attr("finish_vertex")(PythonVertex(_gi, u));
    }

private:
    python::object _gi, _vis;
};

struct do_dfs
{
    template <class Graph, class VertexIndexMap>
    void operator()(Graph& g, VertexIndexMap vertex_index, size_t s,
                    DFSVisitorWrapper vis) const
    {
        typename property_map_type::apply<default_color_type,
                                          VertexIndexMap>::type
            color(vertex_index);
        depth_first_visit(g, vertex(s, g), vis, color);
    }
};


void dfs_search(GraphInterface& g, python::object gi, size_t s,
                python::object vis)
{
    run_action<graph_tool::detail::all_graph_views,mpl::true_>()
        (g, std::bind(do_dfs(), placeholders::_1, g.GetVertexIndex(),
                      s, DFSVisitorWrapper(gi, vis)))();
}

void export_dfs()
{
    using namespace boost::python;
    def("dfs_search", &dfs_search);
}
