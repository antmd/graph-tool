// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006  Tiago de Paula Peixoto <tiago@forked.de>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
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

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

//==============================================================================
// GetAssortativityCoefficient(type)
//==============================================================================

template <class DegreeSelector>
struct get_assortativity_coefficient
{
    get_assortativity_coefficient(DegreeSelector& deg): _deg(deg) {}

    template <class Graph>
    void operator()(const Graph& g, double& r, double& r_err) const
    {
        size_t n_edges = 0;
        int e_kk = 0;
        tr1::unordered_map<double,int> a, b;
        
        typename graph_traits<Graph>::edge_iterator e, e_begin, e_end;
        tie(e_begin,e_end) = edges(g);
        for (e = e_begin; e != e_end; ++e)
        {
            double k1, k2;
            k1 = _deg(source(*e,g),g);
            k2 = _deg(target(*e,g),g);
            a[k1]++;
            b[k2]++;
            if (k1 == k2)
        	e_kk++;
            n_edges++;
            if(is_convertible<typename graph_traits<Graph>::directed_category, undirected_tag>::value)
            {
        	a[k2]++;
        	a[k1]++;
        	if (k1 == k2)
        	    e_kk++;
        	n_edges++;
            }
        }
        
        double t1=double(e_kk)/n_edges, t2=0.0;
        
        for (typeof(a.begin()) iter = a.begin(); iter != a.end(); ++iter)
            if (b.find(iter->second) != b.end())
        	t2 += double(iter->second * b[iter->first]);
        t2 /= n_edges*n_edges;
        
        r = (t1 - t2)/(1.0 - t2);

        // "jackknife" variance
        r_err = 0.0;
        for (e = e_begin; e != e_end; ++e)
        {
            double k1, k2;
            k1 = _deg(source(*e,g),g);
            k2 = _deg(target(*e,g),g);
            int one = 1;
            if(is_convertible<typename graph_traits<Graph>::directed_category, undirected_tag>::value)
        	one = 2; // nice! :-)
            double tl2 = (t2*(n_edges*n_edges) - b[k1] - a[k2])/((n_edges-one)*(n_edges-one));
            double tl1 = t1*n_edges;
            if (k1==k2)
        	tl1 -= one;
            tl1 /= n_edges - one;
            double rl = (tl1 - tl2)/(1.0 - tl2);
            r_err += (r-rl)*(r-rl);
        }
        r_err = sqrt(r_err);
    }
    
    DegreeSelector& _deg;
};

struct choose_assortativity_coefficient
{
    choose_assortativity_coefficient(const GraphInterface& g, GraphInterface::deg_t deg, double& a, double& a_err)
        : _g(g), _a(a), _a_err(a_err) 
    {
        tie(_deg, _deg_name) = get_degree_type(deg);
    }
    
    
    template <class DegreeSelector>
    void operator()(DegreeSelector)
    {        
        if (mpl::at<degree_selector_index, DegreeSelector>::type::value == _deg)
        {
            DegreeSelector deg(_deg_name, _g);
            check_filter(_g, bind<void>(get_assortativity_coefficient<DegreeSelector>(deg), _1, var(_a), var(_a_err)), 
        		 reverse_check(),directed_check());
        }
    };


    const GraphInterface& _g;
    double& _a;
    double& _a_err;
    GraphInterface::degree_t _deg;
    string _deg_name;
};


pair<double,double>
GraphInterface::GetAssortativityCoefficient(GraphInterface::deg_t deg) const
{
    double a, a_err;
    try
    {
        typedef mpl::vector<in_degreeS, out_degreeS, total_degreeS, scalarS> degrees;
        mpl::for_each<degrees>(choose_assortativity_coefficient(*this, deg, a, a_err));
    }
    catch (dynamic_get_failure &e)
    {
        throw GraphException("error getting scalar property: " + string(e.what()));
    }

    return make_pair(a, a_err);
}


//==============================================================================
// GetScalarAssortativityCoefficient(type)
//==============================================================================

template <class DegreeSelector>
struct get_scalar_assortativity_coefficient
{
    get_scalar_assortativity_coefficient(DegreeSelector& deg): _deg(deg) {}

    template <class Graph>
    void operator()(const Graph& g, double& r, double& r_err) const
    {
        size_t n_edges = 0;
        double e_xy = 0.0;
        tr1::unordered_map<double,int> a, b;
        
        typename graph_traits<Graph>::edge_iterator e, e_begin, e_end;
        tie(e_begin,e_end) = edges(g);
        for (e = e_begin; e != e_end; ++e)
        {
            double k1, k2;
            k1 = _deg(source(*e,g),g);
            k2 = _deg(target(*e,g),g);
            a[k1]++;
            b[k2]++;
            e_xy += k1*k2;
            n_edges++;
            if(is_convertible<typename graph_traits<Graph>::directed_category, undirected_tag>::value)
            {
        	a[k2]++;
        	b[k1]++;
        	e_xy += k1*k2;
        	n_edges++;
            }
        }
        
        double t1 = e_xy/n_edges;
        double avg_a = GetHistogramMean(a), avg_b = GetHistogramMean(b);
        double sa = GetHistogramDeviation(a,avg_a), sb = GetHistogramDeviation(b,avg_b);
        
        if (sa*sb > 0)
            r = (t1 - avg_a*avg_b)/(sa*sb);
        else
            r = (t1 - avg_a*avg_b);

        // "jackknife" variance
        r_err = 0.0;
        double diff_a = 0.0, diff_b = 0.0;
        for (typeof(a.begin()) iter = a.begin(); iter != a.end(); ++iter)
            diff_a += (iter->first - avg_a)*iter->second;
        for (typeof(b.begin()) iter = b.begin(); iter != b.end(); ++iter)
            diff_b += (iter->first - avg_b)*iter->second;
        for (e = e_begin; e != e_end; ++e)
        {
            double k1, k2;
            k1 = _deg(source(*e,g),g);
            k2 = _deg(target(*e,g),g);
            int one = 1;
            if(is_convertible<typename graph_traits<Graph>::directed_category, undirected_tag>::value)
        	one = 2;
            double t1l = (e_xy - k1*k2)/(n_edges-one);
            double avg_al = (avg_a*n_edges - k1)/(n_edges-one), avg_bl = (avg_b*n_edges - k2)/(n_edges-one);
            double sal = sa - 2*diff_a*(avg_al-avg_a) + (avg_al-avg_a)*(avg_al-avg_a);
            double sbl = sb - 2*diff_b*(avg_bl-avg_b) + (avg_bl-avg_b)*(avg_bl-avg_b);
            double rl;
            if (sal*sbl > 0)
        	rl = (t1l - avg_al*avg_bl)/(sal*sbl);
            else
        	rl = (t1l - avg_al*avg_bl);
            r_err += (r-rl)*(r-rl);
        }
        r_err = sqrt(r_err);
    }
    
    DegreeSelector& _deg;
};

struct choose_scalar_assortativity_coefficient
{
    choose_scalar_assortativity_coefficient(const GraphInterface& g, GraphInterface::deg_t deg, double& a, double& a_err)
        : _g(g), _a(a), _a_err(a_err) 
    {
        tie(_deg, _deg_name) = get_degree_type(deg);
    }
    
    
    template <class DegreeSelector>
    void operator()(DegreeSelector)
    {        
        if (mpl::at<degree_selector_index, DegreeSelector>::type::value == _deg)
        {
            DegreeSelector deg(_deg_name, _g);
            check_filter(_g, bind<void>(get_scalar_assortativity_coefficient<DegreeSelector>(deg), _1, var(_a), var(_a_err)), 
        		 reverse_check(),directed_check());
        }
    };


    const GraphInterface& _g;
    double& _a;
    double& _a_err;
    GraphInterface::degree_t _deg;
    string _deg_name;
};


pair<double,double>
GraphInterface::GetScalarAssortativityCoefficient(GraphInterface::deg_t deg) const
{
    double a, a_err;
    try
    {
        typedef mpl::vector<in_degreeS, out_degreeS, total_degreeS, scalarS> degrees;
        mpl::for_each<degrees>(choose_scalar_assortativity_coefficient(*this, deg, a, a_err));
    }
    catch (dynamic_get_failure &e)
    {
        throw GraphException("error getting scalar property: " + string(e.what()));
    }

    return make_pair(a, a_err);
}
