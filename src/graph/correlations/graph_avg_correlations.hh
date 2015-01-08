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

#ifndef GRAPH_AVG_CORR_HH
#define GRAPH_AVG_CORR_HH

#include "graph_correlations.hh"

namespace graph_tool
{
using namespace std;
using namespace boost;

// retrieves the generalized correlation
template <class GetDegreePair>
struct get_avg_correlation
{
    get_avg_correlation(python::object& avg, python::object& dev,
                        const vector<long double>& bins,
                        python::object& ret_bins)
        : _avg(avg), _dev(dev), _bins(bins), _ret_bins(ret_bins) {}

    template <class Graph, class DegreeSelector1, class DegreeSelector2,
              class WeightMap>
    void operator()(Graph& g, DegreeSelector1 deg1, DegreeSelector2 deg2,
                    WeightMap weight) const
    {
        GetDegreePair put_point;

        typedef typename DegreeSelector1::value_type type1;
        typedef typename DegreeSelector2::value_type type2;

        typedef typename graph_tool::detail::
            select_float_and_larger::apply<type2,double>::type
            avg_type;
        typedef type1 val_type;

        typedef typename property_traits<WeightMap>::value_type count_type;
        typedef Histogram<type1,count_type,1> count_t;
        typedef Histogram<val_type,avg_type,1> sum_t;

        std::array<vector<val_type>,1> bins;
        bins[0].resize(_bins.size());
        clean_bins(_bins, bins[0]);

        sum_t sum(bins);
        sum_t sum2(bins);
        count_t count(bins);

        SharedHistogram<sum_t> s_sum(sum);
        SharedHistogram<sum_t> s_sum2(sum2);
        SharedHistogram<count_t> s_count(count);

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) \
            firstprivate(s_sum, s_sum2, s_count) schedule(runtime) if (N > 100)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            put_point(v, deg1, deg2, g, weight, s_sum, s_sum2, s_count);
        }
        s_sum.Gather();
        s_sum2.Gather();
        s_count.Gather();

        for (i = 0; i < int(sum.GetArray().size()); ++i)
        {
            sum.GetArray()[i] /= count.GetArray()[i];
            sum2.GetArray()[i] =
                sqrt(abs(sum2.GetArray()[i]/count.GetArray()[i] -
                         sum.GetArray()[i] * sum.GetArray()[i])) /
                sqrt(count.GetArray()[i]);
        }

        bins = sum.GetBins();
        python::list ret_bins;
        ret_bins.append(wrap_vector_owned(bins[0]));
        _ret_bins = ret_bins;
        _avg = wrap_multi_array_owned<avg_type,1>(sum.GetArray());
        _dev = wrap_multi_array_owned<avg_type,1>(sum2.GetArray());
    }
    python::object& _avg;
    python::object& _dev;
    const vector<long double>& _bins;
    python::object& _ret_bins;
};

} // graph_tool namespace

#endif
