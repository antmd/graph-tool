// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007  Tiago de Paula Peixoto <tiago@forked.de>
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

#include <boost/lambda/bind.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

struct get_kruskal_min_span_tree
{
    template <class Graph, class IndexMap, class WeightMap, class TreePropMap>
    void operator()(const Graph* gp, IndexMap vertex_index, WeightMap weights,
                    TreePropMap tree_map) const
    {
        const Graph& g = *gp;
        typedef vector<typename graph_traits<Graph>::edge_descriptor>
            tree_edges_t;
        tree_edges_t tree_edges;
        back_insert_iterator<tree_edges_t> tree_inserter(tree_edges);

        kruskal_minimum_spanning_tree(g, tree_inserter, weight_map(weights));

        typename graph_traits<Graph>::edge_iterator e, e_end;
        for(tie(e, e_end) = edges(g); e != e_end; ++e)
            tree_map[*e] = 0;
        for(typeof(tree_edges.begin()) te = tree_edges.begin();
            te != tree_edges.end(); ++te)
            tree_map[*te] = 1;
    }
};

void GraphInterface::GetMinimumSpanningTree(string weight, string property)
{
    boost::any weight_map, tree_map;

    try
    {
        tree_map = prop(property, _edge_index, _properties);
    }
    catch (property_not_found)
    {
        typedef vector_property_map<bool, edge_index_map_t> tree_map_t;
        tree_map_t new_tree_map(_edge_index);
        _properties.property(property, new_tree_map);
        tree_map = new_tree_map;
    }

    if(weight != "")
    {
        try
        {
            weight_map = prop(weight, _edge_index, _properties);
        }
        catch (property_not_found)
        {
            throw GraphException("weight edge property " + weight +
                                 " not found");
        }
    }
    else
    {
        weight_map = ConstantPropertyMap<size_t,edge_t>(1);
    }

    bool directed = _directed;
    _directed = false;
    typedef mpl::push_back<edge_scalar_properties,
                           ConstantPropertyMap<size_t,edge_t> >::type
        weight_maps;
    run_action<detail::never_directed>()
        (*this, bind<void>(get_kruskal_min_span_tree(), _1, _vertex_index,
                           _2, _3),
         weight_maps(), edge_scalar_properties())(weight_map, tree_map);

    _directed = directed;
}
