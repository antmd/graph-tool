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
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>

#include "graph.hh"
#include "histogram.hh"
#include "graph_filtering.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

struct get_kruskal_min_span_tree
{
    template <class Graph, class IndexMap, class WeightMap, class TreePropMap>
    void operator()(Graph& g, IndexMap vertex_index, WeightMap weights, TreePropMap tree_map) const
    {
        typedef vector<typename graph_traits<Graph>::edge_descriptor> tree_edges_t;
        tree_edges_t tree_edges;
        back_insert_iterator<tree_edges_t> tree_inserter(tree_edges);

        HashedDescriptorMap<IndexMap, size_t> ranks(vertex_index);
        HashedDescriptorMap<IndexMap, typename graph_traits<Graph>::vertex_descriptor> preds(vertex_index);

        kruskal_minimum_spanning_tree(g, tree_inserter, weight_map(weights).rank_map(ranks).predecessor_map(preds));

        typename graph_traits<Graph>::edge_iterator e, e_end;
        for(tie(e, e_end) = edges(g); e != e_end; ++e)
            tree_map[*e] = 0;
        for(typeof(tree_edges.begin()) te = tree_edges.begin(); te != tree_edges.end(); ++te)
            tree_map[*te] = 1;
    }
};

void GraphInterface::GetMinimumSpanningTree(string weight, string property)
{
    typedef vector_property_map<size_t, edge_index_map_t> tree_map_t;
    tree_map_t tree_map(_edge_index);

    bool directed = _directed;
    _directed = false;

    if(weight != "")
    {
        try 
        {
            dynamic_property_map& weight_prop = find_property_map(_properties, weight, typeid(graph_traits<multigraph_t>::edge_descriptor));
            if (get_static_property_map<vector_property_map<double,edge_index_map_t> >(&weight_prop))
            {
                vector_property_map<double,edge_index_map_t> weight_map = 
                    get_static_property_map<vector_property_map<double,edge_index_map_t> >(weight_prop);
                check_filter(*this, bind<void>(get_kruskal_min_span_tree(), _1, var(_vertex_index),var(weight_map),var(tree_map)), 
                             reverse_check(), always_undirected());
            }
            else
            {
                DynamicPropertyMapWrap<double,graph_traits<multigraph_t>::edge_descriptor> weight_map(weight_prop);
                check_filter(*this, bind<void>(get_kruskal_min_span_tree(), _1, var(_vertex_index),var(weight_map),var(tree_map)), 
                             reverse_check(), always_undirected());
            }
        }
        catch (property_not_found& e)
        {
            throw GraphException("error getting scalar property: " + string(e.what()));
        }
    }
    else
    {
        ConstantPropertyMap<double,graph_traits<multigraph_t>::edge_descriptor> weight_map(1.0);
        check_filter(*this, bind<void>(get_kruskal_min_span_tree(), _1, var(_vertex_index),var(weight_map),var(tree_map)), 
                     reverse_check(), always_undirected());
    }

    _directed = directed;

    try
    {
        find_property_map(_properties, property, typeid(graph_traits<multigraph_t>::edge_descriptor));
        RemoveVertexProperty(property);
    }
    catch (property_not_found) {}

    _properties.property(property, tree_map);
}
