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
    void operator()(const Graph &g, double &a) const
    {
	tr1::unordered_map<double,int> e_kk;
	size_t n_edges = 0, n_vertices = 0;
	tr1::unordered_map<double,int> n_k;
	
	typename graph_traits<Graph>::vertex_iterator v, v_begin, v_end;
	tie(v_begin,v_end) = vertices(g);
	for(v = v_begin; v != v_end; ++v)
	{
	    double k1, k2;
	    k1 = _deg(*v,g);
	    n_k[k1]++;
	    n_vertices++;
	    
	    typename graph_traits<Graph>::out_edge_iterator e, e_begin, e_end;
	    tie(e_begin, e_end) = out_edges(*v,g);
	    for (e = e_begin; e != e_end; ++e)
	    {
		k2 = _deg(target(*e,g),g);
		if (k1 == k2)
		    e_kk[k1]++;
		n_edges++;
	    }
	}
	
	double t1=0.0, t2=0.0;
	
	for (typeof(e_kk.begin()) iter = e_kk.begin(); iter != e_kk.end(); ++iter)
	    t1 += double(iter->second)/n_edges;
	
	double avg_k = GetHistogramMean(n_k);
	
	for (typeof(n_k.begin()) iter = n_k.begin(); iter != n_k.end(); ++iter)
	    t2 += double(iter->second * iter->second * iter->first * iter->first)/(avg_k*avg_k*n_vertices*n_vertices);
	
	a = (t1 == 1.0)?1.0:(t1 - t2)/(1.0 - t2);
    }
    
    DegreeSelector& _deg;
};

struct choose_assortativity_coefficient
{
    choose_assortativity_coefficient(const GraphInterface &g, GraphInterface::deg_t deg, double &a)
	: _g(g), _a(a) 
    {
	tie(_deg, _deg_name) = get_degree_type(deg);
    }
    
    
    template <class DegreeSelector>
    void operator()(DegreeSelector)
    {	
	if (mpl::at<degree_selector_index, DegreeSelector>::type::value == _deg)
	{
	    DegreeSelector deg(_deg_name, _g);
	    check_filter(_g, bind<void>(get_assortativity_coefficient<DegreeSelector>(deg), _1, var(_a)), 
			 reverse_check(),directed_check());
	}
    };


    const GraphInterface &_g;
    double &_a;
    GraphInterface::degree_t _deg;
    string _deg_name;
};


double 
GraphInterface::GetAssortativityCoefficient(GraphInterface::deg_t deg) const
{
    double a;
    try
    {
	typedef mpl::vector<in_degreeS, out_degreeS, total_degreeS, scalarS> degrees;
	mpl::for_each<degrees>(choose_assortativity_coefficient(*this, deg, a));
    }
    catch (dynamic_get_failure &e)
    {
	throw GraphException("error getting scalar property: " + string(e.what()));
    }

    return a;
}
