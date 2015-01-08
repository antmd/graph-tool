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
        #pragma omp parallel for default(shared) private(i) schedule(runtime) if (N > 100)
        for (i = 0; i < N; ++i)
        {
            dist_map[i].clear();
            dist_map[i].resize(num_vertices(g), 0);
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
        (gi, std::bind(do_all_pairs_search(), placeholders::_1,
                       gi.GetVertexIndex(), placeholders::_2, placeholders::_3,
                       dense),
         vertex_scalar_vector_properties(),
         mpl::push_back<edge_scalar_properties,cweight_map_t>::type())
        (dist_map, weight);
}

void export_all_dists()
{
    python::def("get_all_dists", &get_all_dists);
};
