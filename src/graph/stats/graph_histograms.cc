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

#include "graph_histograms.hh"

#include <boost/python.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

struct deg_check : public boost::static_visitor<>
{
    void operator()(boost::any a) const
    {
        if (!belongs<vertex_scalar_properties>()(a))
            throw ValueException("Vertex property must be of scalar type.");
    }

    template<class T>
    void operator()(const T&) const {}
};


// this will return the vertex histogram of degrees or scalar properties
python::object
get_vertex_histogram(GraphInterface& gi, GraphInterface::deg_t deg,
                     const vector<long double>& bins)
{
    boost::apply_visitor(deg_check(), deg);

    python::object hist;
    python::object ret_bins;

    run_action<>()(gi, get_histogram<VertexHistogramFiller>(hist, bins,
                                                            ret_bins),
         scalar_selectors())(degree_selector(deg));
    return python::make_tuple(hist, ret_bins);
}

// this will return the vertex histogram of degrees or scalar properties
python::object
get_edge_histogram(GraphInterface& gi, boost::any prop,
                   const vector<long double>& bins)
{
    if (!belongs<edge_scalar_properties>()(prop))
        throw ValueException("Edge property must be of scalar type.");

    python::object hist;
    python::object ret_bins;

    bool directed = gi.GetDirected();
    gi.SetDirected(true);
    run_action<graph_tool::detail::always_directed>()
        (gi, get_histogram<EdgeHistogramFiller>(hist, bins, ret_bins),
         edge_scalar_properties())(prop);
    gi.SetDirected(directed);

    return python::make_tuple(hist, ret_bins);
}

using namespace boost::python;

void export_histograms()
{
    def("get_vertex_histogram", &get_vertex_histogram);
    def("get_edge_histogram", &get_edge_histogram);
}
