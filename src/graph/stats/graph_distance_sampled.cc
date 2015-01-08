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

#include "graph_filtering.hh"
#include "graph.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"

#include "graph_distance_sampled.hh"

#include "random.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

typedef Histogram<size_t, size_t, 1> hist_t;

python::object sampled_distance_histogram(GraphInterface& gi, boost::any weight,
                                          const vector<long double>& bins,
                                          size_t n_samples, rng_t& rng)
{
    python::object ret;

    if (weight.empty())
    {
        run_action<>()(gi,
                       std::bind(get_sampled_distance_histogram(),
                                 placeholders::_1,
                                 gi.GetVertexIndex(), no_weightS(),
                                 n_samples, std::ref(bins),
                                 std::ref(ret), std::ref(rng)))();
    }
    else
    {
        run_action<>()(gi,
                       std::bind(get_sampled_distance_histogram(),
                                 placeholders::_1,
                                 gi.GetVertexIndex(), placeholders::_2,
                                 n_samples,
                                 std::ref(bins), std::ref(ret),
                                 std::ref(rng)),
                           edge_scalar_properties())(weight);
    }
    return ret;
}

void export_sampled_distance()
{
    python::def("sampled_distance_histogram", &sampled_distance_histogram);
}
