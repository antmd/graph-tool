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
        typedef typename Hist::point_t::value_type val_t;
        k[0] = val_t(deg1(v, g));
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
                    Deg1& deg1, Deg2& deg2, Graph& g, WeightMap&,
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
            typename mpl::and_<std::is_floating_point<Type1>,
                               std::is_floating_point<Type2> >::type,
            typename select_larger_type::apply<Type1,Type2>::type,
            typename mpl::if_<std::is_floating_point<Type1>,
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

} // graph_tool namespace

#endif
