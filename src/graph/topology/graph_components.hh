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

#ifndef GRAPH_COMPONENTS_HH
#define GRAPH_COMPONENTS_HH

#include <boost/graph/connected_components.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/biconnected_components.hpp>

namespace graph_tool
{
template <class PropertyMap>
class HistogramPropertyMap;
}

namespace boost
{
template <class PropertyMap>
struct property_traits<graph_tool::HistogramPropertyMap<PropertyMap> >:
        public property_traits<PropertyMap> {};
}

namespace graph_tool
{
using namespace std;
using namespace boost;

// this wraps an existing property map, and makes a simple histogram of the
// values (assumed to be an integer type)

template <class PropertyMap>
class HistogramPropertyMap
{
public:
    typedef typename property_traits<PropertyMap>::value_type value_type;
    typedef typename property_traits<PropertyMap>::key_type key_type;
    typedef typename property_traits<PropertyMap>::category category;

    HistogramPropertyMap(PropertyMap base_map, size_t max, vector<size_t>& hist)
        : _base_map(base_map), _max(max), _hist(hist) {}
    HistogramPropertyMap(){}

    value_type get(const key_type& k) const
    {
        return boost::get(_base_map, k);
    }

    void put(const key_type& k, const value_type& v)
    {
        boost::put(_base_map, k, v);

        vector<size_t>& h = _hist;
        size_t bin = v;
        if (bin > _max)
            return;
        if (bin >= h.size())
            h.resize(bin + 1);
        ++h[bin];
    }

private:
    PropertyMap _base_map;
    size_t _max;
    boost::reference_wrapper<vector<size_t> > _hist;
};

template <class PropertyMap>
typename property_traits<PropertyMap>::value_type
get(const HistogramPropertyMap<PropertyMap>& pmap,
    const typename property_traits<PropertyMap>::key_type& k)
{
    return pmap.get(k);
}

template <class PropertyMap>
void put(HistogramPropertyMap<PropertyMap> pmap,
         const typename property_traits<PropertyMap>::key_type& k,
         const typename property_traits<PropertyMap>::value_type& val)
{
    pmap.put(k, val);
}


// this will label the components of a graph to a given vertex property, from
// [0, number of components - 1], and keep an histogram. If the graph is
// directed the strong components are used.
struct label_components
{
    template <class Graph, class CompMap>
    void operator()(Graph& g, CompMap comp_map, vector<size_t>& hist)
        const
    {
        typedef typename graph_traits<Graph>::directed_category
            directed_category;
        HistogramPropertyMap<CompMap> cm(comp_map, num_vertices(g), hist);
        get_components(g, cm,
                       typename std::is_convertible<directed_category,
                                                    directed_tag>::type());
    }

    template <class Graph, class CompMap>
    void get_components(Graph& g, CompMap comp_map,
                        std::true_type) const
    {
        boost::strong_components(g, comp_map);
    }

    template <class Graph, class CompMap>
    void get_components(Graph& g, CompMap comp_map,
                        std::false_type) const
    {
        boost::connected_components(g, comp_map);
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
    void operator()(Graph& g, CompMap comp_map, ArtMap art_map,
                    vector<size_t>& hist) const
    {
        HistogramPropertyMap<CompMap> cm(comp_map, num_edges(g), hist);
        biconnected_components(g, cm,
                               vertex_inserter<ArtMap>(art_map));
    }
};


struct label_out_component
{
    template <class CompMap>
    class marker_visitor: public bfs_visitor<>
    {
    public:
        marker_visitor() { }
        marker_visitor(CompMap comp) : _comp(comp) { }

        template <class Vertex, class Graph>
        void discover_vertex(Vertex u, const Graph&)
        {
            _comp[u] = true;
        }
    private:
        CompMap _comp;
    };

    template <class Graph, class CompMap>
    void operator()(Graph& g, CompMap comp_map, size_t root) const
    {
        marker_visitor<CompMap> marker(comp_map);
        breadth_first_search(g, vertex(root, g), visitor(marker));
    }
};

struct label_attractors
{
    template <class Graph, class CompMap, class AttrVec>
    void operator()(Graph& g, CompMap comp_map, AttrVec attr)
        const
    {
        typedef typename property_traits<CompMap>::value_type c_type;
        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) schedule(runtime) if (N > 100)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v =
                vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            c_type c = get(comp_map, v);
            if (attr[size_t(c)] == false)
                continue;

            typename graph_traits<Graph>::adjacency_iterator a, a_end;
            for (tie(a, a_end) = adjacent_vertices(v, g); a != a_end; ++a)
            {
                if (get(comp_map, *a) != c)
                {
                    attr[size_t(c)] = false;
                    break;
                }
            }
        }
    }
};

} // graph_tool namespace

#endif // GRAPH_COMPONENTS_HH
