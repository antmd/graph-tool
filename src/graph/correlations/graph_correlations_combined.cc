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

#include <boost/python.hpp>

#include "graph.hh"
#include "histogram.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"

#include "graph_corr_hist.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

typedef ConstantPropertyMap<int,GraphInterface::edge_t> dummy_weight;

boost::python::object
get_vertex_combined_correlation_histogram(GraphInterface& gi,
                                          GraphInterface::deg_t deg1,
                                          GraphInterface::deg_t deg2,
                                          const vector<long double>& xbin,
                                          const vector<long double>& ybin)
{
    boost::python::object hist;
    boost::python::object ret_bins;

    std::array<vector<long double>,2> bins;
    bins[0] = xbin;
    bins[1] = ybin;

    run_action<>()(gi, get_correlation_histogram<GetCombinedPair>(hist, bins,
                                                                  ret_bins),
                   scalar_selectors(), scalar_selectors(),
                   boost::mpl::vector<dummy_weight>())
        (any(degree_selector(deg1)), any(degree_selector(deg2)),
         boost::any(dummy_weight()));

    return boost::python::make_tuple(hist, ret_bins);
}

using namespace boost::python;

void export_combined_vertex_correlations()
{
    def("vertex_combined_correlation_histogram",
        &get_vertex_combined_correlation_histogram);
}
