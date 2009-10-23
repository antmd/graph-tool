// Copyright (C) 2008  Tiago de Paula Peixoto <tiago@forked.de>
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
#include "graph_properties.hh"
#include "graph_selectors.hh"

#include <boost/python.hpp>

#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

struct do_all_pairs_search
{
    template <class Graph, class VertexIndexMap, class DistMap, class WeightMap>
    void operator()(const Graph& g, VertexIndexMap vertex_index,
                    DistMap dist_map, WeightMap weight, bool dense) const
    {
        typedef typename property_traits<DistMap>::value_type::value_type
            dist_t;

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            dist_map[v].clear();
            dist_map[v].resize(num_vertices(g), 0);
        }

        if (dense)
        {
            floyd_warshall_all_pairs_shortest_paths
                (g, dist_map,
                 weight_map(ConvertedPropertyMap<WeightMap,dist_t>(weight)).
                 vertex_index_map(vertex_index));
        }
        else
        {
            johnson_all_pairs_shortest_paths
                (g, dist_map,
                 weight_map(ConvertedPropertyMap<WeightMap,dist_t>(weight)).
                 vertex_index_map(vertex_index));
        }
    }
};

void get_all_dists(GraphInterface& gi, boost::any dist_map, boost::any weight,
                   bool dense)
{
    typedef ConstantPropertyMap<size_t,GraphInterface::edge_t> cweight_map_t;

    if (weight.empty())
        weight = boost::any(cweight_map_t(1));

    run_action<>()
        (gi,
         bind<void>(do_all_pairs_search(), _1, gi.GetVertexIndex(),
                    _2, _3, dense),
         vertex_scalar_vector_properties(),
         mpl::push_back<edge_scalar_properties,cweight_map_t>::type())
        (dist_map, weight);
}

void export_all_dists()
{
    python::def("get_all_dists", &get_all_dists);
};

