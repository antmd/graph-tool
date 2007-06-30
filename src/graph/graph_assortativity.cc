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
        SharedMap<tr1::unordered_map<double,int> > sa(a), sb(b);

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) firstprivate(sa,sb) schedule(dynamic) reduction(+:e_kk, n_edges)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            double k1 = _deg(v, g);
            typename graph_traits<Graph>::out_edge_iterator e, e_end;
            for (tie(e,e_end) = out_edges(v, g); e != e_end; ++e)
            {
                double k2 = _deg(target(*e, g), g);
                sa[k1]++;
                sb[k2]++;
                if (k1 == k2)
                    e_kk++;
                n_edges++;
            }
        }

        sa.Gather();
        sb.Gather();

        double t1=double(e_kk)/n_edges, t2=0.0;
        
        for (typeof(a.begin()) iter = a.begin(); iter != a.end(); ++iter)
            if (b.find(iter->second) != b.end())
                t2 += double(iter->second * b[iter->first]);
        t2 /= n_edges*n_edges;
        
        r = (t1 - t2)/(1.0 - t2);

        // "jackknife" variance
        double err = 0.0;
        #pragma omp parallel for default(shared) private(i) schedule(dynamic) reduction(+:err)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            double k1 = _deg(v, g);
            typename graph_traits<Graph>::out_edge_iterator e, e_end;
            for (tie(e,e_end) = out_edges(v, g); e != e_end; ++e)
            {
                double k2 = _deg(target(*e,g),g);                
                double tl2 = (t2*(n_edges*n_edges) - b[k1] - a[k2])/((n_edges-1)*(n_edges-1));
                double tl1 = t1*n_edges;
                if (k1 == k2)
                    tl1 -= 1;
                tl1 /= n_edges - 1;
                double rl = (tl1 - tl2)/(1.0 - tl2);
                err += (r-rl)*(r-rl);
            }
        }
        r_err = sqrt(err);
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
        SharedMap<tr1::unordered_map<double,int> > sa(a), sb(b);

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) firstprivate(sa,sb) schedule(dynamic) reduction(+:e_xy, n_edges)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            double k1 = _deg(v, g);
            typename graph_traits<Graph>::out_edge_iterator e, e_end;
            for (tie(e,e_end) = out_edges(v, g); e != e_end; ++e)
            {
                double k2 = _deg(target(*e,g),g);
                sa[k1]++;
                sb[k2]++;
                e_xy += k1*k2;
                n_edges++;
            }
        }

        sa.Gather();
        sb.Gather();

         
        double t1 = e_xy/n_edges;
        double avg_a = GetHistogramMean(a), avg_b = GetHistogramMean(b);
        double da = GetHistogramDeviation(a,avg_a), db = GetHistogramDeviation(b,avg_b);
        
        if (da*db > 0)
            r = (t1 - avg_a*avg_b)/(da*db);
        else
            r = (t1 - avg_a*avg_b);

        // "jackknife" variance
        r_err = 0.0;
        double diff_a = 0.0, diff_b = 0.0;
        for (typeof(a.begin()) iter = a.begin(); iter != a.end(); ++iter)
            diff_a += (iter->first - avg_a)*iter->second;
        for (typeof(b.begin()) iter = b.begin(); iter != b.end(); ++iter)
            diff_b += (iter->first - avg_b)*iter->second;
        
        double err = 0.0;
        #pragma omp parallel for default(shared) private(i) schedule(dynamic) reduction(+:err)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            double k1 = _deg(v, g);
            typename graph_traits<Graph>::out_edge_iterator e, e_end;
            for (tie(e,e_end) = out_edges(v, g); e != e_end; ++e)
            {
                double k2 = _deg(target(*e, g), g);                
                double t1l = (e_xy - k1*k2)/(n_edges-1);
                double avg_al = (avg_a*n_edges - k1)/(n_edges-1), avg_bl = (avg_b*n_edges - k2)/(n_edges-1);
                double dal = da - 2*diff_a*(avg_al-avg_a) + (avg_al-avg_a)*(avg_al-avg_a);
                double dbl = db - 2*diff_b*(avg_bl-avg_b) + (avg_bl-avg_b)*(avg_bl-avg_b);
                double rl;
                if (dal*dbl > 0)
                    rl = (t1l - avg_al*avg_bl)/(dal*dbl);
                else
                    rl = (t1l - avg_al*avg_bl);
                err += (r-rl)*(r-rl);
            }
        }
        r_err = sqrt(err);
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
