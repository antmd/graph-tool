// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007-2010 Tiago de Paula Peixoto <tiago@forked.de>
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

#ifndef GRAPH_PRICE_HH
#define GRAPH_PRICE_HH

#include <iostream>
#include <boost/functional/hash.hpp>
#include "graph_util.hh"

#if (GCC_VERSION >= 40400)
#   include <tr1/random>
#   include <tr1/unordered_set>
#else
#   include <boost/tr1/random.hpp>
#   include <boost/tr1/unordered_set>
#endif

#include <map>
#include <iostream>

namespace graph_tool
{
using namespace std;
using namespace boost;

typedef tr1::mt19937 rng_t;

struct get_price
{
    template <class Graph>
    void operator()(Graph& g, size_t N, double gamma, double c, size_t m,
                    size_t seed) const
    {
        rng_t rng(static_cast<rng_t::result_type>(seed));

        typedef typename mpl::if_<typename is_directed::apply<Graph>::type,
                                  in_degreeS, out_degreeS>::type DegSelector;

        map<double, typename graph_traits<Graph>::vertex_descriptor> probs;

        double p = 0;
        typename graph_traits<Graph>::vertex_iterator vi, vi_end;
        for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
        {
            p += pow(DegSelector()(*vi, g) + c, gamma);
            probs.insert(make_pair(p, *vi));
        }

        if (probs.rbegin()->first <= 0)
            throw GraphException("Cannot connect edges: probabilities are <= 0!");

        for (size_t i = 0; i < N; ++i)
        {
            unordered_set<typename graph_traits<Graph>::vertex_descriptor>
                visited;
            typename graph_traits<Graph>::vertex_descriptor v = add_vertex(g);
            for (size_t j = 0; j < m; ++j)
            {
                tr1::variate_generator<rng_t&, tr1::uniform_real<> >
                    sample(rng, tr1::uniform_real<>(0, probs.rbegin()->first));
                double r = sample();
                typeof(probs.begin()) iter = probs.lower_bound(r);
                typename graph_traits<Graph>::vertex_descriptor w =
                    iter->second;

                if (visited.find(w) != visited.end())
                {
                    --j;
                    continue;
                }
                visited.insert(w);
                add_edge(v, w, g);

                p = abs(pow(DegSelector()(w, g) + c, gamma)
                        - pow(DegSelector()(w, g) + c - 1, gamma));
                probs.insert(make_pair(probs.rbegin()->first + p, w));
            }
            p = pow(DegSelector()(v, g) + c, gamma);
            probs.insert(make_pair(probs.rbegin()->first + p, v));
        }
    }
};

} // namespace graph_tool

#endif // GRAPH_PRICE_HH
