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

#ifndef GRAPH_COMPONENTS_HH
#define GRAPH_COMPONENTS_HH

#include <boost/graph/connected_components.hpp>
#include <boost/graph/strong_components.hpp>

namespace graph_tool
{
using namespace std;
using namespace boost;

// this will label the components of a graph to a given vertex property, from
// [0, number of components - 1]. If the graph is directed the "strong
// components" are used.
struct label_components
{
    template <class Graph, class CompMap>
    void operator()(const Graph* gp, CompMap comp_map) const
    {
        const Graph& g = *gp;
        typedef typename graph_traits<Graph>::directed_category
            directed_category;
        get_components(g, comp_map,
                       typename is_convertible<directed_category,
                                               directed_tag>::type());
    }

    template <class Graph, class CompMap>
    void get_components(const Graph& g, CompMap comp_map,
                        boost::true_type is_directed) const
    {
        strong_components(g, comp_map);
    }

    template <class Graph, class CompMap>
    void get_components(const Graph& g, CompMap comp_map,
                        boost::false_type is_directed) const
    {
        connected_components(g, comp_map);
    }
};

} // graph_tool namespace

#endif // GRAPH_COMPONENTS_HH
