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

#ifndef GRAPH_HISTOGRAMS_HH
#define GRAPH_HISTOGRAMS_HH

#include <algorithm>
#include <boost/numeric/conversion/bounds.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/python/object.hpp>
#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>
#include "numpy_bind.hh"
#include "histogram.hh"
#include "shared_map.hh"

namespace graph_tool
{
using namespace std;
using namespace boost;

class VertexHistogramFiller
{
public:
    template <class Graph, class DegreeSelector, class Hist>
    void operator()(const Graph& g,
                    typename graph_traits<Graph>::vertex_descriptor v,
                    DegreeSelector& deg, Hist& hist)
    {
        typename Hist::point_t p;
        p[0] = deg(v, g);
        hist.PutValue(p);
    }
};

class EdgeHistogramFiller
{
public:
    template <class Graph, class EdgeProperty, class Hist>
    void operator()(const Graph& g,
                    typename graph_traits<Graph>::vertex_descriptor v,
                    EdgeProperty& eprop, Hist& hist)
    {
        typename Hist::point_t p;
        typename graph_traits<Graph>::out_edge_iterator e, e_begin, e_end;
        tie(e_begin,e_end) = out_edges(v,g);
        for(e = e_begin; e != e_end; ++e)
        {
            p[0] = eprop[*e];
            hist.PutValue(p);
        }
    }
};

// generalized functor to obtain histogram of different types of "degrees"
template <class HistogramFiller>
struct get_histogram
{
    get_histogram(python::object& hist, const vector<long double>& bins,
                  python::object& ret_bins)
        : _hist(hist), _bins(bins), _ret_bins(ret_bins) {}

    template <class Graph, class DegreeSelector>
    void operator()(const Graph& g, DegreeSelector deg) const
    {
        typedef typename DegreeSelector::value_type value_type;
        typedef Histogram<value_type, size_t, 1> hist_t;

        HistogramFiller filler;

        vector<value_type> bins(_bins.size());
        for (size_t i = 0; i < bins.size(); ++i)
        {
            // we'll attempt to recover from out of bounds conditions
            try
            {
                bins[i] = numeric_cast<value_type,long double>(_bins[i]);
            }
            catch (boost::numeric::negative_overflow&)
            {
                bins[i] = boost::numeric::bounds<value_type>::lowest();
            }
            catch (boost::numeric::positive_overflow&)
            {
                bins[i] = boost::numeric::bounds<value_type>::highest();
            }
        }

        // sort the bins
        sort(bins.begin(), bins.end());

        // clean bins of zero size
        vector<value_type> temp_bin(1);
        temp_bin[0] = bins[0];
        for (size_t j = 1; j < bins.size(); ++j)
        {
            if (bins[j] > bins[j-1])
                temp_bin.push_back(bins[j]);
        }
        bins = temp_bin;

        std::array<vector<value_type>, 1> bin_list;
        bin_list[0] = bins;

        hist_t hist(bin_list);
        SharedHistogram<hist_t> s_hist(hist);

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) \
            firstprivate(s_hist) schedule(runtime) if (N > 100)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            filler(g, v, deg, s_hist);
        }
        s_hist.Gather();

        bin_list = hist.GetBins();
        python::object ret_bins = wrap_vector_owned(bin_list[0]);
        _ret_bins = ret_bins;
        _hist = wrap_multi_array_owned<size_t,1>(hist.GetArray());
    }
    python::object& _hist;
    const vector<long double>& _bins;
    python::object& _ret_bins;
};

} // graph_tool namespace

#endif // GRAPH_HISTOGRAMS_HH
