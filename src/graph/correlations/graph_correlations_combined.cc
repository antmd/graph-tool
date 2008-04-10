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
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "graph_filtering.hh"

#include <boost/lambda/bind.hpp>
#include <boost/python.hpp>

#include "graph.hh"
#include "histogram.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"

#include "graph_correlations.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

typedef ConstantPropertyMap<int,GraphInterface::edge_t> dummy_weight;

python::object
get_vertex_combined_correlation_histogram(const GraphInterface& gi,
                                          GraphInterface::deg_t deg1,
                                          GraphInterface::deg_t deg2,
                                          const vector<long double>& xbin,
                                          const vector<long double>& ybin)
{
    python::object hist;
    python::object ret_bins;

    array<vector<long double>,2> bins;
    bins[0] = xbin;
    bins[1] = ybin;

    run_action<>()(gi, lambda::bind<void>
                   (get_correlation_histogram<GetCombinedPair>(hist, bins,
                                                               ret_bins),
                    lambda::_1, lambda::_2, lambda::_3, dummy_weight(0)),
                   all_selectors(), all_selectors())
        (degree_selector(deg1, gi), degree_selector(deg2, gi));

    return python::make_tuple(hist, ret_bins);
}

using namespace boost::python;

void export_combined_vertex_correlations()
{
    def("vertex_combined_correlation_histogram",
        &get_vertex_combined_correlation_histogram);
}
