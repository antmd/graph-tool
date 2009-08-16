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
#include <boost/graph/biconnected_components.hpp>

namespace graph_tool
{
using namespace std;
using namespace boost;

// this will label the components of a graph to a given vertex property, from
// [0, number of components - 1]. If the graph is directed the strong
// components are used.
struct label_components
{
    template <class Graph, class CompMap>
    void operator()(const Graph& g, CompMap comp_map) const
    {
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

struct label_biconnected_components
{
    template <class ArtMap>
    class vertex_inserter
    {
    public:
        vertex_inserter(ArtMap art_map): _art_map(art_map) {}

        vertex_inserter& operator++() { return *this; }
        vertex_inserter& operator++(int) { return *this; }
        vertex_inserter& operator*() { return *this; }

        vertex_inserter& operator=
        (const typename property_traits<ArtMap>::key_type& e)
        {
            put(_art_map, e, 1);
            return *this;
        }

    private:
        ArtMap _art_map;
    };

    template <class Graph, class CompMap, class ArtMap>
    void operator()(const Graph& g, CompMap comp_map, ArtMap art_map,
                    size_t& nc) const
    {
        nc = biconnected_components(g, comp_map,
                                    vertex_inserter<ArtMap>(art_map)).first;
    }
};


} // graph_tool namespace

#endif // GRAPH_COMPONENTS_HH
