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

#ifndef GRAPH_SIMILARITY_HH
#define GRAPH_SIMILARITY_HH

#include <unordered_set>


namespace graph_tool
{
using namespace std;
using namespace boost;


template <class Keys, class Set1, class Set2>
size_t intersection_size(Keys& ks, Set1& s1, Set2& s2)
{
    size_t s = 0;
    for (auto k : ks)
    {
        int c1 = s1.count(k);
        int c2 = s2.count(k);
        s += max(c1, c2) - abs(c1 - c2);
    }
    return s;
}


struct get_similarity
{
    template <class Graph1, class Graph2, class LabelMap>
    void operator()(const Graph1& g1, const Graph2* g2p, LabelMap l1,
                    boost::any al2, size_t& s) const
    {
        const Graph2& g2 = *g2p;

        LabelMap l2 = boost::any_cast<typename LabelMap::checked_t>(al2).get_unchecked(num_vertices(g2));

        typedef typename property_traits<LabelMap>::value_type label_t;

        std::unordered_map<label_t, typename graph_traits<Graph1>::vertex_descriptor,
                           boost::hash<label_t>>
            lmap1;
        std::unordered_map<label_t, typename graph_traits<Graph2>::vertex_descriptor,
                           boost::hash<label_t>>
            lmap2;

        for (auto v : vertices_range(g1))
            lmap1[get(l1, v)] = v;
        for (auto v : vertices_range(g2))
            lmap2[get(l2, v)] = v;

        s = 0;
        for (auto& lv1 : lmap1)
        {
            auto v1 = lv1.second;

            auto li2 = lmap2.find(lv1.first);
            if (li2 == lmap2.end())
                continue;
            auto v2 = li2->second;

            std::unordered_set<label_t, boost::hash<label_t>> keys;
            std::unordered_multiset<label_t, boost::hash<label_t>> adj1;
            std::unordered_multiset<label_t, boost::hash<label_t>> adj2;

            for (auto a1 : adjacent_vertices_range(v1, g1))
            {
                adj1.insert(get(l1, a1));
                keys.insert(get(l1, a1));
            }

            for (auto a2 : adjacent_vertices_range(v2, g2))
            {
                adj2.insert(get(l2, a2));
                keys.insert(get(l2, a2));
            }

            s += intersection_size(keys, adj1, adj2);
        }
    }
};

struct get_similarity_fast
{
    template <class Graph1, class Graph2, class LabelMap>
    void operator()(const Graph1& g1, const Graph2* g2p, LabelMap l1,
                    boost::any al2, size_t& s) const
    {
        const Graph2& g2 = *g2p;

        LabelMap l2 = boost::any_cast<LabelMap>(al2);

        typedef typename property_traits<LabelMap>::value_type label_t;

        vector<typename graph_traits<Graph1>::vertex_descriptor> lmap1;
        vector<typename graph_traits<Graph1>::vertex_descriptor> lmap2;

        for (auto v : vertices_range(g1))
        {
            size_t i = get(l1, v);
            if (lmap1.size() <= i)
                lmap1.resize(i + 1);
            lmap1[i] = v;
        }

        for (auto v : vertices_range(g2))
        {
            size_t i = get(l2, v);
            if (lmap2.size() <= i)
                lmap2.resize(i + 1);
            lmap2[i] = v;
        }

        size_t ss = 0;

        int i, N = lmap1.size();
        #pragma omp parallel for default(shared) private(i) schedule(runtime) \
            reduction(+:ss) if (N > 100)
        for (i = 0; i < N; ++i)
        {
            auto v1 = lmap1[i];
            auto v2 = lmap2[i];

            std::unordered_set<label_t> keys;
            std::unordered_multiset<label_t> adj1;
            std::unordered_multiset<label_t> adj2;

            for (auto a1 : adjacent_vertices_range(v1, g1))
            {
                adj1.insert(get(l1, a1));
                keys.insert(get(l1, a1));
            }

            for (auto a2 : adjacent_vertices_range(v2, g2))
            {
                adj2.insert(get(l2, a2));
                keys.insert(get(l2, a2));
            }

            ss += intersection_size(keys, adj1, adj2);
        }

        s = ss;
    }
};

} // graph_tool namespace

#endif // GRAPH_SIMILARITY_HH
