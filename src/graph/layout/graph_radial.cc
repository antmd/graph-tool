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

#include "graph_selectors.hh"
#include "graph_properties.hh"

#include <cmath>

using namespace std;
using namespace boost;
using namespace graph_tool;


struct do_get_radial
{
    template <class Graph, class PosProp, class LevelMap, class OrderMap>
    void operator()(Graph& g, PosProp tpos, LevelMap level, OrderMap order,
                    size_t root, bool weighted, double r) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        typedef property_map_type::apply<int, GraphInterface::vertex_index_map_t>::type
            vcount_t;
        vcount_t::unchecked_t count(get(vertex_index, g), num_vertices(g));

        if (!weighted)
        {
            typename graph_traits<Graph>::vertex_iterator v, v_end;
            for(tie(v, v_end) = vertices(g); v != v_end; ++v)
                count[*v] = 1;
        }
        else
        {
            deque<vertex_t> q;
            typename graph_traits<Graph>::vertex_iterator v, v_end;
            for(tie(v, v_end) = vertices(g); v != v_end; ++v)
                if (out_degree(*v, g) == 0)
                {
                    q.push_back(*v);
                    count[*v] = 1;
                }

            typedef property_map_type::apply<uint8_t, GraphInterface::vertex_index_map_t>::type
                vmark_t;
            vmark_t::unchecked_t mark(get(vertex_index, g), num_vertices(g));

            while (!q.empty())
            {
                vertex_t v = q.front();
                q.pop_front();
                typename graph_traits<Graph>::in_edge_iterator e, e_end;
                for(tie(e, e_end) = in_edges(v, g); e != e_end; ++e)
                {
                    vertex_t w = source(*e, g);
                    count[w] += count[v];
                    if (!mark[w])
                    {
                        q.push_back(w);
                        mark[w] = true;
                    }
                }
            }
        }

        vector<vector<vertex_t> > layers(1);
        layers[0].push_back(root);

        bool last = false;
        while (!last)
        {
            layers.resize(layers.size() + 1);
            vector<vertex_t>& new_layer = layers[layers.size() - 1];
            vector<vertex_t>& last_layer = layers[layers.size() - 2];

            last = true;
            for (size_t i = 0; i < last_layer.size(); ++i)
            {
                vertex_t v = last_layer[i];
                typename graph_traits<Graph>::out_edge_iterator e, e_end;
                for(tie(e, e_end) = out_edges(v, g); e != e_end; ++e)
                {
                    vertex_t w = target(*e, g);
                    new_layer.push_back(w);

                    if (int(layers.size()) - 1 == int(level[w]))
                        last = false;
                }

                std::sort(new_layer.end() - out_degree(v, g), new_layer.end(),
                          [&] (vertex_t u, vertex_t v) -> bool { return order[u] < order[v]; });

                if (out_degree(v, g) == 0)
                    new_layer.push_back(v);
            }

            if (last)
                layers.pop_back();
        }


        typedef property_map_type::apply<double, GraphInterface::vertex_index_map_t>::type
            vangle_t;
        vangle_t::unchecked_t angle(get(vertex_index, g), num_vertices(g));

        double d_sum = 0;
        vector<vertex_t>& outer_layer = layers.back();
        for (size_t i = 0; i < outer_layer.size(); ++i)
            d_sum += count[outer_layer[i]];
        for (size_t i = 0; i < outer_layer.size(); ++i)
            angle[outer_layer[i]] = (i * 2 * M_PI * count[outer_layer[i]]) / d_sum;

        for (size_t i = 0; i < layers.size(); ++i)
        {
            vector<vertex_t>& vs = layers[layers.size() - 1 - i];
            for (size_t j = 0; j < vs.size(); ++j)
            {
                vertex_t v = vs[j];
                d_sum = 0;
                typename graph_traits<Graph>::out_edge_iterator e, e_end;
                for(tie(e, e_end) = out_edges(v, g); e != e_end; ++e)
                {
                    vertex_t w = target(*e, g);
                    d_sum += count[w];
                }
                for(tie(e, e_end) = out_edges(v, g); e != e_end; ++e)
                {
                    vertex_t w = target(*e, g);
                    angle[v] += angle[w] * count[w] / d_sum;
                }
                double d = level[v] * r;
                tpos[v].resize(2);
                tpos[v][0] = d * cos(angle[v]);
                tpos[v][1] = d * sin(angle[v]);
            }
        }
    }
};

void get_radial(GraphInterface& gi, boost::any otpos, boost::any olevels,
                boost::any oorder, size_t root, bool weighted, double r)
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        vmap_t;

    vmap_t levels = boost::any_cast<vmap_t>(olevels);

    run_action<graph_tool::detail::always_directed>()
        (gi, std::bind(do_get_radial(), placeholders::_1, placeholders::_2,
                       levels, placeholders::_3, root, weighted,  r),
         vertex_scalar_vector_properties(),
         vertex_properties())(otpos, oorder);
}

#include <boost/python.hpp>

void export_radial()
{
    python::def("get_radial", &get_radial);
}
