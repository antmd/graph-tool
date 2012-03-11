// Copyright (C) 2008  Tiago de Paula Peixoto <tiago@skewed.de>
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

#if (GCC_VERSION >= 40400)
#   include <tr1/unordered_set>
#else
#   include <boost/tr1/unordered_set.hpp>
#endif

namespace graph_tool
{
using namespace std;
using namespace boost;


template <class Keys, class Set>
size_t intersection_size(Keys& ks, Set& s1, Set& s2)
{
    size_t s = 0;
    for (typeof(ks.begin()) k = ks.begin(); k != ks.end(); ++k)
    {
        int c1 = s1.count(*k);
        int c2 = s2.count(*k);
        s += max(c1, c2) - abs(c1 - c2);
    }
    return s;
}


struct get_similarity
{
    template <class Graph1, class Graph2, class LabelMap>
    void operator()(const Graph1& g1, const Graph2* g2p, LabelMap l1,
                    any l2a, size_t& s) const
    {
        LabelMap l2 = any_cast<LabelMap>(l2a);
        const Graph2& g2 = *g2p;

        typedef typename property_traits<LabelMap>::value_type label_t;

        tr1::unordered_map<label_t, typename graph_traits<Graph1>::vertex_descriptor>
            lmap1;
        tr1::unordered_map<label_t, typename graph_traits<Graph2>::vertex_descriptor>
            lmap2;

        typename graph_traits<Graph1>::vertex_iterator v1, v1_end;
        for (tie(v1, v1_end) = vertices(g1); v1 != v1_end; ++v1)
            lmap1[get(l1, *v1)] = *v1;
        typename graph_traits<Graph2>::vertex_iterator v2, v2_end;
        for (tie(v2, v2_end) = vertices(g2); v2 != v2_end; ++v2)
            lmap2[get(l2, *v2)] = *v2;

        s = 0;
        for (typeof(lmap1.begin()) li = lmap1.begin(); li != lmap1.end(); ++li)
        {
            typename graph_traits<Graph1>::vertex_descriptor v1 = li->second;

            typeof(lmap2.begin()) li2 = lmap2.find(li->first);
            if (li2 == lmap2.end())
                continue;
            typename graph_traits<Graph2>::vertex_descriptor v2 = li2->second;

            tr1::unordered_set<label_t> keys;
            tr1::unordered_multiset<label_t> adj1;
            tr1::unordered_multiset<label_t> adj2;

            typename graph_traits<Graph1>::adjacency_iterator a1, a1_end;
            for(tie(a1, a1_end) = adjacent_vertices(v1, g1); a1 != a1_end; ++a1)
            {
                adj1.insert(get(l1, *a1));
                keys.insert(get(l1, *a1));
            }

            typename graph_traits<Graph2>::adjacency_iterator a2, a2_end;
            for(tie(a2, a2_end) = adjacent_vertices(v2, g2); a2 != a2_end; ++a2)
            {
                adj2.insert(get(l2, *a2));
                keys.insert(get(l2, *a2));
            }

            s += intersection_size(keys, adj1, adj2);
        }
    }
};

} // graph_tool namespace

#endif // GRAPH_SIMILARITY_HH
