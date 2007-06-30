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
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#include <algorithm>
#include <tr1/unordered_set>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/random.hpp>

#include "graph.hh"
#include "histogram.hh"
#include "graph_filtering.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"
#include "shared_map.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

//==============================================================================
// GetCombinedVertexHistogram()
// retrieves the distribution of combined (deg1,deg2) degrees
//==============================================================================

template <class DegreeSelector1, class DegreeSelector2>
struct get_combined_degree_histogram
{
    get_combined_degree_histogram(DegreeSelector1& deg1, DegreeSelector2& deg2): _deg1(deg1), _deg2(deg2) {}

    template <class Graph, class Hist>
    void operator()(Graph &g, Hist &hist) const
    {
        SharedMap<Hist> s_hist(hist);

        typedef typename Hist::key_type::first_type first_type;
        typedef typename Hist::key_type::second_type second_type;

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) firstprivate(s_hist) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            s_hist[make_pair(first_type(_deg1(v,g)), second_type(_deg2(v,g)))]++;
        }

        s_hist.Gather();
    }
    DegreeSelector1& _deg1;
    DegreeSelector2& _deg2;
};

template <class SecondDegreeSelectors>
struct choose_combined_degree_histogram
{
    choose_combined_degree_histogram(const GraphInterface &g, GraphInterface::deg_t deg1, 
                                     GraphInterface::deg_t deg2, GraphInterface::hist2d_t &hist)
        : _g(g), _hist(hist) 
    {
        tie(_deg1, _deg_name1) = get_degree_type(deg1);
        tie(_deg2, _deg_name2) = get_degree_type(deg2);
    }

    template <class DegreeSelector1>
    struct check_second_degree
    {
        check_second_degree(choose_combined_degree_histogram<SecondDegreeSelectors> &parent):_parent(parent) {}
        template <class DegreeSelector2>
        void operator()(DegreeSelector2)
        {
            if (mpl::at<degree_selector_index,DegreeSelector2>::type::value == _parent._deg2)
            {
                DegreeSelector1 deg1(_parent._deg_name1, _parent._g);
                DegreeSelector2 deg2(_parent._deg_name2, _parent._g);
                check_filter(_parent._g, bind<void>(get_combined_degree_histogram<DegreeSelector1,DegreeSelector2>(deg1,deg2), _1, var(_parent._hist)), 
                             reverse_check(),directed_check());
            }
        }
        choose_combined_degree_histogram<SecondDegreeSelectors>& _parent;
    };

    template <class DegreeSelector>
    void operator()(DegreeSelector)
    {
        if (mpl::at<degree_selector_index,DegreeSelector>::type::value == _deg1)
            mpl::for_each<SecondDegreeSelectors>(check_second_degree<DegreeSelector>(*this));
    }
    const GraphInterface &_g;
    GraphInterface::hist2d_t &_hist;
    GraphInterface::degree_t _deg1;
    string _deg_name1;
    GraphInterface::degree_t _deg2;
    string _deg_name2;
    
};

GraphInterface::hist2d_t
GraphInterface::GetCombinedVertexHistogram(deg_t deg1, deg_t deg2) const
{
    hist2d_t hist;
    typedef mpl::vector<in_degreeS, out_degreeS, total_degreeS, scalarS> degree_selectors;
    try 
    {
        mpl::for_each<degree_selectors>(choose_combined_degree_histogram<degree_selectors>(*this, deg1, deg2, hist));
    }
    catch (dynamic_get_failure &e)
    {
        throw GraphException("error getting scalar property: " + string(e.what()));
    }
    return hist;
}


//==============================================================================
// GetAverageCombinedVertexCorrelation()
// retrieves the average of deg2 in function of deg1
//==============================================================================

template <class DegreeSelector1, class DegreeSelector2>
struct get_average_combined_degree_correlation
{
    get_average_combined_degree_correlation(DegreeSelector1& deg1, DegreeSelector2& deg2): _deg1(deg1), _deg2(deg2) {}

    template <class Graph, class AvgDeg>
    void operator()(Graph &g, AvgDeg &avg_deg) const
    {
        tr1::unordered_map<double,size_t> count;

        typename graph_traits<Graph>::vertex_iterator v, v_begin, v_end;
        tie(v_begin, v_end) = vertices(g);
        for(v = v_begin; v != v_end; ++v)
        {
            typename AvgDeg::key_type d1 = _deg1(*v,g);
            typename AvgDeg::value_type::second_type::first_type d2 = _deg2(*v,g);
            avg_deg[d1].first += d2;
            avg_deg[d1].second += d2*d2;
            count[d1]++;
        }

        for (typeof(avg_deg.begin()) iter = avg_deg.begin(); iter != avg_deg.end(); ++iter)
        {
            size_t N = count[iter->first];
            iter->second.first /= N;
            if (N > 1)
            {
                double err = (iter->second.second - N*iter->second.first*iter->second.first)/(N*(N-1));
                iter->second.second = (err<0.0)?0.0:sqrt(err);
            }
            else
            {
                iter->second.second = 0.0;
            }
        }
    }
    DegreeSelector1& _deg1;
    DegreeSelector2& _deg2;
};

template <class SecondDegreeSelectors>
struct choose_average_combined_degree_correlation
{
    choose_average_combined_degree_correlation(const GraphInterface &g, GraphInterface::deg_t deg1, 
                                               GraphInterface::deg_t deg2, GraphInterface::avg_corr_t &avg_corr)
        : _g(g), _avg_corr(avg_corr) 
    {
        tie(_deg1, _deg_name1) = get_degree_type(deg1);
        tie(_deg2, _deg_name2) = get_degree_type(deg2);
    }

    template <class DegreeSelector1>
    struct check_second_degree
    {
        check_second_degree(choose_average_combined_degree_correlation<SecondDegreeSelectors> &parent):_parent(parent) {}
        template <class DegreeSelector2>
        void operator()(DegreeSelector2)
        {
            if (mpl::at<degree_selector_index,DegreeSelector2>::type::value == _parent._deg2)
            {
                DegreeSelector1 deg1(_parent._deg_name1, _parent._g);
                DegreeSelector2 deg2(_parent._deg_name2, _parent._g);
                check_filter(_parent._g, bind<void>(get_average_combined_degree_correlation<DegreeSelector1,DegreeSelector2>(deg1,deg2), 
                                                    _1, var(_parent._avg_corr)), 
                             reverse_check(),directed_check());
            }
        }
        choose_average_combined_degree_correlation<SecondDegreeSelectors>& _parent;
    };

    template <class DegreeSelector>
    void operator()(DegreeSelector)
    {
        if (mpl::at<degree_selector_index,DegreeSelector>::type::value == _deg1)
            mpl::for_each<SecondDegreeSelectors>(check_second_degree<DegreeSelector>(*this));
    }
    const GraphInterface &_g;
    GraphInterface::avg_corr_t &_avg_corr;
    GraphInterface::degree_t _deg1;
    string _deg_name1;
    GraphInterface::degree_t _deg2;
    string _deg_name2;
    
};

GraphInterface::avg_corr_t
GraphInterface::GetAverageCombinedVertexCorrelation(deg_t deg1, deg_t deg2) const
{
    avg_corr_t avg_corr;
    typedef mpl::vector<in_degreeS, out_degreeS, total_degreeS, scalarS> degree_selectors;
    try 
    {
        mpl::for_each<degree_selectors>(choose_average_combined_degree_correlation<degree_selectors>(*this, deg1, deg2, avg_corr));
    }
    catch (dynamic_get_failure &e)
    {
        throw GraphException("error getting scalar property: " + string(e.what()));
    }
    return avg_corr;
}
