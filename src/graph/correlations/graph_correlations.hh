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

#ifndef GRAPH_CORRELATIONS_HH
#define GRAPH_CORRELATIONS_HH

#include <algorithm>
#include <boost/numeric/conversion/bounds.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/python/object.hpp>
#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>
#include "histogram.hh"
#include "numpy_bind.hh"
#include "shared_map.hh"

namespace graph_tool
{
using namespace std;
using namespace boost;

// get degrees pairs from source and of neighbours
class GetNeighboursPairs
{
public:

    template <class Graph, class Deg1, class Deg2, class Hist, class WeightMap>
    void operator()(typename graph_traits<Graph>::vertex_descriptor v,
                    Deg1& deg1, Deg2& deg2, Graph& g, WeightMap& weight,
                    Hist& hist)
    {
        typename Hist::point_t k;
        k[0] = deg1(v, g);
        typename graph_traits<Graph>::out_edge_iterator e, e_end;
        for (tie(e,e_end) = out_edges(v, g); e != e_end; ++e)
        {
            k[1] = deg2(target(*e,g),g);
            hist.PutValue(k, get(weight, *e));
        }
    }

    template <class Graph, class Deg1, class Deg2, class Sum, class Count,
              class WeightMap>
    void operator()(typename graph_traits<Graph>::vertex_descriptor v,
                    Deg1& deg1, Deg2& deg2, Graph& g, WeightMap& weight,
                    Sum& sum, Sum& sum2, Count& count)
    {
        typename Sum::point_t k1;
        k1[0] = deg1(v, g);
        typename Sum::count_type k2;
        typename graph_traits<Graph>::out_edge_iterator e, e_end;
        for (tie(e,e_end) = out_edges(v, g); e != e_end; ++e)
        {
            k2 = deg2(target(*e,g),g)*get(weight, *e);
            sum.PutValue(k1, k2);
            sum2.PutValue(k1, k2*k2);
            count.PutValue(k1, get(weight, *e));
        }
    }
};

// get degrees pairs from one single vertex
class GetCombinedPair
{
public:

    template <class Graph, class Deg1, class Deg2, class Hist, class Dummy>
    void operator()(typename graph_traits<Graph>::vertex_descriptor v,
                    Deg1& deg1, Deg2& deg2, Graph& g, const Dummy&,
                    Hist& hist)
    {
        typename Hist::point_t k;
        k[0] = deg1(v, g);
        k[1] = deg2(v, g);
        hist.PutValue(k);
    }

    template <class Graph, class Deg1, class Deg2, class Sum, class Count,
              class WeightMap>
    void operator()(typename graph_traits<Graph>::vertex_descriptor v,
                    Deg1& deg1, Deg2& deg2, Graph& g, WeightMap& weight,
                    Sum& sum, Sum& sum2, Count& count)
    {
        typename Sum::point_t k1;
        k1[0] = deg1(v, g);
        typename Sum::count_type k2;
        k2 = deg2(v, g);
        sum.PutValue(k1, k2);
        sum2.PutValue(k1, k2*k2);
        count.PutValue(k1, 1);
    }
};


namespace detail
{
struct select_larger_type
{
    template <class Type1, class Type2>
    struct apply
    {
        typedef typename mpl::if_<
            typename mpl::greater<typename mpl::sizeof_<Type1>::type,
                                  typename mpl::sizeof_<Type2>::type>::type,
            Type1,
            Type2>::type type;
    };
};

struct select_float_and_larger
{
    template <class Type1, class Type2>
    struct apply
    {
        typedef typename mpl::if_<
            typename mpl::and_<is_floating_point<Type1>,
                               is_floating_point<Type2> >::type,
            typename select_larger_type::apply<Type1,Type2>::type,
            typename mpl::if_<is_floating_point<Type1>,
                              Type1,
                              Type2>::type>::type type;
    };
};

}

template <class Value>
void clean_bins(const vector<long double>& obins, vector<Value>& rbins)
{
    typedef Value val_type;
    rbins.resize(obins.size());
    for (size_t j = 0; j < rbins.size(); ++j)
    {
        // we'll attempt to recover from out of bounds conditions
        try
        {
            rbins[j] = numeric_cast<val_type,long double>(obins[j]);
        }
        catch (boost::numeric::negative_overflow&)
        {
            rbins[j] = boost::numeric::bounds<val_type>::lowest();
        }
        catch (boost::numeric::positive_overflow&)
        {
            rbins[j] = boost::numeric::bounds<val_type>::highest();
        }
    }
    // sort the bins
    sort(rbins.begin(), rbins.end());
    // clean bins of zero size
    vector<val_type> temp_bin(1);
    temp_bin[0] = rbins[0];
    for (size_t j = 1; j < rbins.size(); ++j)
    {
        if (rbins[j] > rbins[j-1])
                    temp_bin.push_back(rbins[j]);
    }
    rbins = temp_bin;
}


// retrieves the generalized vertex-vertex correlation histogram
template <class GetDegreePair>
struct get_correlation_histogram
{
    get_correlation_histogram(python::object& hist,
                              const array<vector<long double>,2>& bins,
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

        array<vector<val_type>,2> bins;
        for (size_t i = 0; i < bins.size(); ++i)
            clean_bins(_bins[i], bins[i]);

        hist_t hist(bins);
        SharedHistogram<hist_t> s_hist(hist);

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) \
            firstprivate(s_hist) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            put_point(v, deg1, deg2, g, weight, s_hist);
        }
        s_hist.Gather();

        bins = hist.GetBins();
        python::list ret_bins;
        ret_bins.append(wrap_vector_owned(bins[0]));
        ret_bins.append(wrap_vector_owned(bins[1]));
        _ret_bins = ret_bins;
        _hist = wrap_multi_array_owned<count_type,2>(hist.GetArray());
    }
    python::object& _hist;
    const array<vector<long double>,2>& _bins;
    python::object& _ret_bins;
};

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

        array<vector<val_type>,1> bins;
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
            firstprivate(s_sum, s_sum2, s_count) schedule(dynamic)
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

        for (size_t i = 0; i < sum.GetArray().size(); ++i)
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
