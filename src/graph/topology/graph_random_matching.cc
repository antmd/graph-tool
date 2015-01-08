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
#include "graph_util.hh"

#include "random.hh"

#include <boost/python.hpp>


using namespace std;
using namespace boost;
using namespace graph_tool;

struct do_random_matching
{
    template <class Graph, class VertexIndex, class WeightMap, class MatchMap,
              class RNG>
    void operator()(const Graph& g, VertexIndex vertex_index, WeightMap weight,
                    MatchMap match, bool minimize, RNG& rng) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;
        typedef typename property_traits<WeightMap>::value_type wval_t;

        vector<vertex_t> vlist;
        typename graph_traits<Graph>::vertex_iterator v, v_end;
        for (tie(v, v_end) = vertices(g); v != v_end; ++v)
            vlist.push_back(*v);

        unchecked_vector_property_map<uint8_t, VertexIndex>
            matched(vertex_index, num_vertices(g));

        typedef random_permutation_iterator<typename vector<vertex_t>::iterator,
                                            RNG>
            random_vertex_iter;
        random_vertex_iter vr(vlist.begin(), vlist.end(), rng),
            vr_end(vlist.end(), vlist.end(), rng);
        for (; vr != vr_end; ++vr)
        {
            vertex_t v = *vr;
            if (matched[v])
                continue;
            wval_t min_w = minimize ? numeric_limits<wval_t>::max() :
                numeric_limits<wval_t>::min() ;
            vector<edge_t> candidates;
            typename graph_traits<Graph>::out_edge_iterator e, e_end;
            for(tie(e, e_end) = out_edges(v, g); e != e_end; ++e)
            {
                vertex_t w = target(*e, g);
                if (matched[w])
                    continue;
                if ((minimize && (weight[*e] < min_w)) ||
                    (!minimize && (weight[*e] > min_w)))
                {
                    min_w = weight[*e];
                    candidates.clear();
                }
                if (weight[*e] == min_w)
                    candidates.push_back(*e);
            }

            if (!candidates.empty())
            {
                uniform_int_distribution<> sample(0, candidates.size() - 1);
                size_t j = sample(rng);
                match[candidates[j]] = true;
                matched[v] = true;
                matched[target(candidates[j], g)] = true;
            }
        }
    }
};

void random_matching(GraphInterface& gi, boost::any weight, boost::any match,
                     bool minimize, rng_t& rng)
{
    typedef ConstantPropertyMap<int32_t,GraphInterface::edge_t> weight_map_t;
    typedef mpl::push_back<edge_scalar_properties, weight_map_t>::type
        edge_props_t;

    if(weight.empty())
        weight = weight_map_t(1);

    run_action<>()
        (gi, std::bind(do_random_matching(), placeholders::_1, gi.GetVertexIndex(),
                       placeholders::_2, placeholders::_3, minimize, std::ref(rng)),
         edge_props_t(), writable_edge_scalar_properties())(weight, match);
}

void export_random_matching()
{
    python::def("random_matching", &random_matching);
}
