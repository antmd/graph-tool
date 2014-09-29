// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2016 Tiago de Paula Peixoto <tiago@skewed.de>
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

#ifndef GRAPH_BLOCKMODEL_COVARIATE_HH
#define GRAPH_BLOCKMODEL_COVARIATE_HH

#include <cmath>
#include <iostream>
#include <queue>
#include <boost/math/special_functions/zeta.hpp>
#include <boost/functional/hash.hpp>

#include "config.h"
#include <unordered_set>
#include <unordered_map>
#include <tuple>

#ifdef HAVE_SPARSEHASH
#include SPARSEHASH_INCLUDE(dense_hash_set)
#include SPARSEHASH_INCLUDE(dense_hash_map)
#endif

namespace graph_tool
{

#ifdef HAVE_SPARSEHASH
using google::dense_hash_set;
using google::dense_hash_map;
#endif

using namespace boost;

// this will label each edge covariate in a consecutive range [0, C-1]

struct ec_hist
{
    template <class Graph, class EVMap, class EMap>
    void operator()(Graph&g, const EVMap& ev, EMap& ec) const
    {
        typedef typename property_traits<EVMap>::value_type val_t;
        unordered_map<val_t, size_t> ehist;
        for (auto e : edges_range(g))
        {
            auto iter = ehist.find(ev[e]);
            size_t x;
            if (iter == ehist.end())
            {
                x = ehist.size();
                ehist[ev[e]] = x;
            }
            else
            {
                x = iter->second;
            }
            ec[e] = x;
        }
    }
};

// this will split an edge-valued graph into layers

struct split_graph
{
    template <class Graph, class ECMap, class EMap, class VMap, class VVMap, class BMap>
    void operator()(Graph& g, ECMap& ec, VMap& b, EMap& eweight, VMap& vweight,
                    VVMap& vc, VVMap& vmap,
                    std::vector<std::reference_wrapper<GraphInterface>>& us,
                    std::vector<std::reference_wrapper<VMap>>& ub,
                    std::vector<std::reference_wrapper<VMap>>& uvweight,
                    std::vector<std::reference_wrapper<EMap>>& ueweight,
                    std::vector<BMap>& block_map,
                    std::vector<std::reference_wrapper<VMap>>& block_rmap,
                    std::vector<std::reference_wrapper<VMap>>& uvmap) const
    {
        std::vector<unordered_map<size_t, size_t>> vhmap(num_vertices(g));

        auto get_v = [&] (size_t v, size_t l) -> size_t
        {
            auto iter = vhmap[v].find(l);
            if (iter == vhmap[v].end())
            {
                size_t u = add_vertex(us[l].get().GetGraph());
                vhmap[v][l] = u;
                size_t pos = lower_bound(vc[v].begin(), vc[v].end(), l) - vc[v].begin();
                vc[v].insert(vc[v].begin() + pos, l);
                vmap[v].insert(vmap[v].begin() + pos, u);
                uvmap[l].get()[u] = v;
                uvweight[l].get()[u] = vweight[v];
                size_t r =  b[v];
                size_t u_r;

                if (block_map.size() <= l + 1)
                {
                    size_t n = block_map.size();
                    block_map.resize(l + 2);
#ifdef HAVE_SPARSEHASH
                    for (size_t i = n; i < block_map.size(); ++i)
                    {
                        block_map[i].set_empty_key(numeric_limits<size_t>::max());
                        block_map[i].set_deleted_key(numeric_limits<size_t>::max() - 1);
                    }
#endif
                }
                auto& bmap = block_map[l + 1];
                auto riter = bmap.find(r);
                if (riter == bmap.end())
                {
                    u_r = bmap.size();
                    bmap[r] = u_r;
                    block_rmap[l].get()[u_r] = r;
                }
                else
                {
                    u_r = riter->second;
                }
                ub[l].get()[u] = u_r;
                return u;
            }
            else
            {
                return iter->second;
            }
        };

        for (auto e : edges_range(g))
        {
            auto s = source(e, g);
            auto t = target(e, g);
            size_t l = ec[e];

            auto u_s = get_v(s, l);
            auto u_t = get_v(t, l);
            auto ne = add_edge(u_s, u_t, us[l].get().GetGraph()).first;
            ueweight[l].get()[ne] = eweight[e];
        }
    }
};

#ifdef HAVE_SPARSEHASH
    typedef vector<dense_hash_map<size_t, size_t, std::hash<size_t>>> bmap_t;
#else
    typedef vector<unordered_map<size_t, size_t>> bmap_t;
#endif


template <class BlockState>
size_t get_block_map(BlockState& state, typename bmap_t::value_type& bmap,
                     size_t r)
{
    size_t r_u;
    #pragma omp critical (covariate_block_map)
    {
        auto iter = bmap.find(r);
        if (iter == bmap.end())
        {
            if (state.free_blocks.empty())
            {
                r_u = bmap.size();
            }
            else
            {
                r_u = state.free_blocks.back();
                state.free_blocks.pop_back();
            }
            bmap[r] = r_u;
            state.block_rmap[r_u] = r;
        }
        else
        {
            r_u = iter->second;
        }
        assert(r_u < num_vertices(state.bg));
    }
    // assert(state.block_rmap[r_u] == r);
    // assert(r_u < num_vertices(state.bg));
    return r_u;
}

template <class BlockState>
void remove_block_map(BlockState& state, typename bmap_t::value_type& bmap,
                      size_t r)
{
    #pragma omp critical (covariate_block_map)
    {
        auto iter = bmap.find(r);
        if (iter != bmap.end()) // another thread may have removed it already
        {
            state.free_blocks.push_back(iter->second);
            bmap.erase(iter);
        }
        // assert(bmap.find(r) == bmap.end());
    }
}


}; // graph_tool namespace

#endif // GRAPH_BLOCKMODEL_COVARIATE_HH
