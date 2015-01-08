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
#include "graph_properties.hh"

#include "random.hh"

#include <boost/graph/random_spanning_tree.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

struct get_random_span_tree
{
    template <class Graph, class IndexMap, class WeightMap, class TreeMap,
              class RNG>
    void operator()(const Graph& g, size_t root, IndexMap vertex_index,
                    WeightMap weights, TreeMap tree_map, RNG& rng) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;

        unchecked_vector_property_map<vertex_t,IndexMap>
            pred_map(vertex_index, num_vertices(g));
        random_spanning_tree(g, rng, predecessor_map(pred_map).
                             root_vertex(vertex(root, g)).
                             weight_map(weights).
                             vertex_index_map(vertex_index));

        // convert the predecessor map to a tree map, and avoid trouble with
        // parallel edges
        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) schedule(runtime) if (N > 100)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            vector<edge_t> edges;
            vector<typename property_traits<WeightMap>::value_type> ws;
            typename graph_traits<Graph>::out_edge_iterator e, e_end;
            for (tie(e,e_end) = out_edges(v, g); e != e_end; ++e)
            {
                if (target(*e,g) == pred_map[v])
                {
                    edges.push_back(*e);
                    ws.push_back(weights[*e]);
                }
            }
            if (!edges.empty())
            {
                edge_t e = *(edges.begin() +
                             size_t(min_element(ws.begin(),
                                                ws.end())-ws.begin()));
                tree_map[e] = 1;
            }
        }
    }
};

typedef property_map_types::apply<mpl::vector<uint8_t>,
                                  GraphInterface::edge_index_map_t,
                                  mpl::bool_<false> >::type
    tree_properties;

void get_random_spanning_tree(GraphInterface& gi, size_t root,
                              boost::any weight_map, boost::any tree_map,
                              rng_t& rng)
{
    typedef ConstantPropertyMap<size_t,GraphInterface::edge_t> cweight_t;

    if (weight_map.empty())
        weight_map = cweight_t(1);

    typedef mpl::push_back<writable_edge_scalar_properties, cweight_t>::type
        weight_maps;

    run_action<>()
        (gi, std::bind(get_random_span_tree(), placeholders::_1, root, gi.GetVertexIndex(),
            placeholders::_2, placeholders::_3, std::ref(rng)),
         weight_maps(), tree_properties())(weight_map, tree_map);
}
