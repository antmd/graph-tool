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

#ifndef GRAPH_CORR_HIST_HH
#define GRAPH_CORR_HIST_HH

#include "graph_correlations.hh"

namespace graph_tool
{
using namespace std;
using namespace boost;

// retrieves the generalized vertex-vertex correlation histogram
template <class GetDegreePair>
struct get_correlation_histogram
{
    get_correlation_histogram(python::object& hist,
                              const std::array<vector<long double>,2>& bins,
                              python::object& ret_bins)
        : _hist(hist), _bins(bins), _ret_bins(ret_bins) {}

    template <class Graph, class DegreeSelector1, class DegreeSelector2,
              class WeightMap>
    void operator()(Graph& g, DegreeSelector1 deg1, DegreeSelector2 deg2,
                    WeightMap weight) const
    {
        GetDegreePair put_point;

        typedef typename DegreeSelector1::value_type type1;
        typedef typename DegreeSelector2::value_type type2;

        typedef typename graph_tool::detail::
            select_float_and_larger::apply<type1,type2>::type
            val_type;
        typedef typename property_traits<WeightMap>::value_type count_type;
        typedef Histogram<val_type, count_type, 2> hist_t;

        std::array<vector<val_type>,2> bins;
        for (size_t i = 0; i < bins.size(); ++i)
            clean_bins(_bins[i], bins[i]);

        hist_t hist(bins);
        {
            SharedHistogram<hist_t> s_hist(hist);

            int i, N = num_vertices(g);
            #pragma omp parallel for default(shared) private(i) \
                firstprivate(s_hist) schedule(runtime) if (N > 100)
            for (i = 0; i < N; ++i)
            {
                typename graph_traits<Graph>::vertex_descriptor v =
                    vertex(i, g);
                if (v == graph_traits<Graph>::null_vertex())
                    continue;
                put_point(v, deg1, deg2, g, weight, s_hist);
            }
        }

        bins = hist.GetBins();
        python::list ret_bins;
        ret_bins.append(wrap_vector_owned(bins[0]));
        ret_bins.append(wrap_vector_owned(bins[1]));
        _ret_bins = ret_bins;
        _hist = wrap_multi_array_owned<count_type,2>(hist.GetArray());
    }
    python::object& _hist;
    const std::array<vector<long double>,2>& _bins;
    python::object& _ret_bins;
};

} // graph_tool namespace

#endif
