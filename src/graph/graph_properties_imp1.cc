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

#include "graph_python_interface.hh"
#include "graph.hh"
#include "graph_properties.hh"
#include "graph_filtering.hh"
#include "graph_selectors.hh"
#include "graph_util.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

template <bool src>
struct do_edge_endpoint
{
    template <class Graph, class EdgeIndexMap, class VertexPropertyMap>
    void operator()(Graph& g, EdgeIndexMap, VertexPropertyMap prop,
                    boost::any aeprop) const
    {
        typedef typename property_traits<VertexPropertyMap>::value_type vval_t;
        typedef typename boost::mpl::if_<std::is_same<vval_t, size_t>, int64_t, vval_t>::type
            val_t;
        typedef typename property_map_type::apply<val_t, EdgeIndexMap>::type
            eprop_t;
        eprop_t eprop = any_cast<eprop_t>(aeprop);
        eprop.reserve(num_edges(g));

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            typename graph_traits<Graph>::out_edge_iterator e, e_end;
            for (auto e : out_edges_range(v, g))
            {
                auto s = v;
                auto t = target(e, g);
                if (!is_directed::apply<Graph>::type::value && s > t)
                    continue;
                if (src)
                    eprop[e] = prop[s];
                else
                    eprop[e] = prop[t];
            }
        }
    }
};

void edge_endpoint(GraphInterface& gi, boost::any prop,
                   boost::any eprop, std::string endpoint)
{
    if (endpoint == "source")
        run_action<>()(gi, bind<void>(do_edge_endpoint<true>(), _1,
                                      gi.GetEdgeIndex(), _2, eprop),
                       vertex_properties())(prop);
    else
        run_action<>()(gi, bind<void>(do_edge_endpoint<false>(), _1,
                                      gi.GetEdgeIndex(), _2, eprop),
                       vertex_properties())(prop);
}
