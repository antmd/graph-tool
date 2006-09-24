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

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/random.hpp>
#include <tr1/unordered_set>

#include "graph.hh"
#include "histogram.hh"
#include "graph_filtering.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

using std::tr1::unordered_map;
using std::tr1::unordered_set;

typedef boost::mt19937 rng_t;

template <class Graph>
size_t out_degree_no_loops(typename graph_traits<Graph>::vertex_descriptor v, const Graph &g)
{
    size_t k = 0;
    typename graph_traits<Graph>::adjacency_iterator a,a_end;
    for (tie(a,a_end) = adjacent_vertices(v,g); a != a_end; ++a)
	if (*a != v)
	    k++;
    return k;
}

//==============================================================================
// GetCommunityStructure()
// computes the community structure through a spin glass system with 
// simulated annealing
//==============================================================================

template <template <class G> class NNKS>
struct get_communities
{
    template <class Graph, class WeightMap, class CommunityMap>
    void operator()(Graph& g, WeightMap weights, CommunityMap s, double gamma, size_t n_iter, size_t seed) const
    {
	rng_t rng(seed);
	boost::uniform_real<double> uniform_p(0.0,1.0);

	unordered_map<size_t, size_t> Nk; // degree histogram
	unordered_map<size_t, size_t> Ns; // spin histogram
	unordered_map<size_t, unordered_map<size_t, size_t> > kNs; // spin histogram per degree
	unordered_map<size_t, size_t> Ks; // sum of degrees per spin
	unordered_map<size_t, map<double, unordered_set<size_t> > > global_term; // global energy term

	NNKS<Graph> Nnnks(g); // this will retrieve the expected number of neighbours with given spin, in funcion of degree

	// init spins from [0,N-1] and global info	
	size_t index = 0;
	typename graph_traits<Graph>::vertex_iterator v,v_end;
	for (tie(v,v_end) = vertices(g); v != v_end; ++v)
	{
	    s[*v] = index++; 
	    Ns[s[*v]]++;
	    size_t k = out_degree_no_loops(*v,g);
	    Nk[k]++;
	    kNs[k][s[*v]]++;
	    Ks[s[*v]] += k;
	}

	for (tie(v,v_end) = vertices(g); v != v_end; ++v)
	{
	    size_t k = out_degree_no_loops(*v,g);
	    global_term[k][gamma*Nnnks(k,s[*v],Ns,kNs,Ks)].insert(s[*v]);
	}

	// define cooling rate so that temperature starts at "infinity" (numeric_limits::max()) at temp_count == 0
	// and reaches "zero" (numeric_limits::epsilon()) at temp_count == n_iter - 1
	double cooling_rate = -(log(numeric_limits<double>::epsilon())-log(numeric_limits<double>::max()))/(n_iter-1);

	for (size_t temp_count = 0; temp_count < n_iter; ++temp_count)
	{
	    double T = numeric_limits<double>::max()*exp(-cooling_rate*temp_count);
	    
            // sample a new spin for every vertex
	    for (tie(v,v_end) = vertices(g); v != v_end; ++v)
	    {
		double E = -T*log(1-uniform_p(rng)); // sampled target jump energy
		double local_term = 0;
		
		unordered_map<size_t, double> ns; // number of neighbours with spin 's' (weighted)
		
		// neighbourhood spins info
		typename graph_traits<Graph>::out_edge_iterator e,e_end;
		for (tie(e,e_end) = out_edges(*v,g); e != e_end; ++e)
		{
		    typename graph_traits<Graph>::vertex_descriptor t = target(*e,g);
		    if (t != *v)
			ns[s[t]] += get(weights, *e);
		}
		    
		size_t k = out_degree_no_loops(*v,g);
		
		// local energy term
		local_term = ns[s[*v]] + gamma*Nnnks(k,s[*v],Ns,kNs,Ks);
		
		// update global info with local info
		for (typeof(ns.begin()) iter = ns.begin(); iter != ns.end(); ++iter)
		{
		    double old_E = gamma*Nnnks(k,iter->first,Ns,kNs,Ks);
		    global_term[k][old_E].erase(iter->first);
		    if (global_term[k][old_E].empty())
			global_term[k].erase(old_E);
		    global_term[k][old_E - ns[iter->first]].insert(iter->first);
		}
		    
		// fetch spins with closest energy
		E -= local_term;
		
		typeof(global_term[k].begin()) lower = global_term[k].lower_bound(E);
		if (lower == global_term[k].end())
		    lower = global_term[k].begin();
		typeof(global_term[k].begin()) upper = global_term[k].upper_bound(E);
		if (upper == global_term[k].end())
		    --upper;
		    
		typeof(global_term[k].begin()) closest = (abs(E - lower->first) < abs(E-upper->first))?lower:upper;
		    
		//new spin (randomly chosen amongst those with equal energy)
		uniform_int<size_t> sample_spin(0,closest->second.size()-1);
		typeof(closest->second.begin()) iter = closest->second.begin();
		advance(iter, sample_spin(rng));
		size_t a = *iter;
		
		//cleanup global info
		for (typeof(ns.begin()) iter = ns.begin(); iter != ns.end(); ++iter)
		{
		    double old_E = gamma*Nnnks(k,iter->first,Ns,kNs,Ks) - ns[iter->first];		    
		    global_term[k][old_E].erase(iter->first);
		    if (global_term[k][old_E].empty())
			global_term[k].erase(old_E);
		    global_term[k][old_E + ns[iter->first]].insert(iter->first);
		}

		//update global info
		double old_E = gamma*Nnnks(k,s[*v],Ns,kNs,Ks);
		global_term[k][old_E].erase(s[*v]);
		if (global_term[k][old_E].empty())
		    global_term[k].erase(old_E);
		Ns[s[*v]]--;
		kNs[k][s[*v]]--;
		Ks[s[*v]] -= k;
		Ns[a]++;
		kNs[k][a]++;
		Ks[a] += k;
		global_term[k][gamma*Nnnks(k,s[*v],Ns,kNs,Ks)].insert(s[*v]);
		global_term[k][gamma*Nnnks(k,a,Ns,kNs,Ks)].insert(a);
		
		// update spin
		s[*v] = a;
	    }		
	}	
    }
};

template <class Graph>
class NNKSErdosReyni
{
public:
    NNKSErdosReyni(Graph &g)
    {
	size_t N = 0;
	double _avg_k = 0.0;
	typename graph_traits<Graph>::vertex_iterator v,v_end;
	for (tie(v,v_end) = vertices(g); v != v_end; ++v)
	{
	    size_t k = out_degree_no_loops(*v,g); 
	    _avg_k += k;
	    N++;
	}
	_p = _avg_k/(N*N);
    }

    double operator()(size_t k, size_t s, unordered_map<size_t,size_t>& Ns, unordered_map<size_t,unordered_map<size_t,size_t> >& kNs,  
		      unordered_map<size_t, size_t>& Ks) const 
    {
	return _p*Ns[s];
    }

private:
    double _p;
};

template <class Graph>
class NNKSUncorr
{
public:
    NNKSUncorr(Graph &g)
	:_K(0)
    {
	typename graph_traits<Graph>::vertex_iterator v,v_end;
	for (tie(v,v_end) = vertices(g); v != v_end; ++v)
	    _K += out_degree_no_loops(*v,g); 
    }

    double operator()(size_t k, size_t s, unordered_map<size_t,size_t>& Ns, unordered_map<size_t,unordered_map<size_t,size_t> >& kNs,  
		      unordered_map<size_t, size_t>& Ks) const 
    {
	return k*Ks[s]/double(_K);
    }
private:
    size_t _K;
};

template <class Graph>
class NNKSCorr
{
public:
    NNKSCorr(Graph &g)
    {
	typename graph_traits<Graph>::vertex_iterator v,v_end;
	for (tie(v,v_end) = vertices(g); v != v_end; ++v)
	{
	    size_t k = out_degree_no_loops(*v,g);
	    _Nk[k]++;
	    typename graph_traits<Graph>::adjacency_iterator a,a_end;
	    for (tie(a,a_end) = adjacent_vertices(*v,g); a != a_end; ++a)
		if (*a != *v)
		    _Nkk[k][out_degree_no_loops(*a,g)]++;
	}
    }

    double operator()(size_t k, size_t s, unordered_map<size_t,size_t>& Ns, unordered_map<size_t,unordered_map<size_t,size_t> >& kNs,  
		      unordered_map<size_t, size_t>& Ks)
    {
	double N = 0;
	unordered_map<size_t, size_t>& nkk = _Nkk[k];
	for (typeof(nkk.begin()) iter = nkk.begin(); iter != nkk.end(); ++iter)
	    N += iter->second * kNs[iter->first][s]/double(_Nk[iter->first]);
	return N;
    }
private:
    unordered_map<size_t, size_t> _Nk;
    unordered_map<size_t, unordered_map<size_t, size_t> > _Nkk;
};


struct get_communities_selector
{
    get_communities_selector(GraphInterface::comm_corr_t corr):_corr(corr) {}
    GraphInterface::comm_corr_t _corr;

    template <class Graph, class WeightMap, class CommunityMap>
    void operator()(Graph& g, WeightMap weights, CommunityMap s, double gamma, size_t n_iter, size_t seed) const
    {
	switch (_corr)
	{
	case GraphInterface::ERDOS_REYNI:
	    get_communities<NNKSErdosReyni>()(g, weights, s, gamma, n_iter, seed);
	    break;
	case GraphInterface::UNCORRELATED:
	    get_communities<NNKSUncorr>()(g, weights, s, gamma, n_iter, seed);
	    break;
	case GraphInterface::CORRELATED:
	    get_communities<NNKSCorr>()(g, weights, s, gamma, n_iter, seed);
	    break;
	}
    }
};

void GraphInterface::GetCommunityStructure(double gamma, comm_corr_t corr, size_t n_iter, size_t seed, string weight, string property)
{
    typedef HashedDescriptorMap<vertex_index_map_t,size_t> comm_map_t;
    comm_map_t comm_map(_vertex_index);

    bool directed = _directed;
    _directed = false;

    if(weight != "")
    {
	try 
	{
	    dynamic_property_map& weight_prop = find_property_map(_properties, weight, typeid(graph_traits<multigraph_t>::edge_descriptor));
	    if (get_static_property_map<vector_property_map<double,edge_index_map_t> >(&weight_prop))
	    {
		vector_property_map<double,edge_index_map_t> weight_map = 
		    get_static_property_map<vector_property_map<double,edge_index_map_t> >(weight_prop);
		check_filter(*this, bind<void>(get_communities_selector(corr), _1, var(weight_map), var(comm_map), gamma, n_iter, seed),
			     reverse_check(), always_undirected());
	    }
	    else
	    {
		DynamicPropertyMapWrap<double,graph_traits<multigraph_t>::edge_descriptor> weight_map(weight_prop);
		check_filter(*this, bind<void>(get_communities_selector(corr), _1, var(weight_map), var(comm_map), gamma, n_iter, seed),
			     reverse_check(), always_undirected());
	    }
	}
	catch (property_not_found& e)
	{
	    throw GraphException("error getting scalar property: " + string(e.what()));
	}
    }
    else
    {
	ConstantPropertyMap<double,graph_traits<multigraph_t>::edge_descriptor> weight_map(1.0);
	check_filter(*this, bind<void>(get_communities_selector(corr), _1, var(weight_map), var(comm_map), gamma, n_iter, seed),
		     reverse_check(), always_undirected());
    }
    _directed = directed;

    try
    {
	find_property_map(_properties, property, typeid(graph_traits<multigraph_t>::edge_descriptor));
	RemoveVertexProperty(property);
    }
    catch (property_not_found) {}

    _properties.property(property, comm_map);
}

//==============================================================================
// GetModularity()
// get Newman's modularity of a given community partition
//==============================================================================

struct get_modularity
{
    template <class Graph, class WeightMap, class CommunityMap>
    void operator()(Graph& g, WeightMap weights, CommunityMap s, double& modularity) const
    {
	modularity = 0.0;
	
	size_t E = 0;
	double W = 0;

	typename graph_traits<Graph>::edge_iterator e,e_end;
	for (tie(e,e_end) = edges(g); e != e_end; ++e)
	    if (target(*e,g) != source(*e,g))
	    {
		W += get(weights, *e);
		E++;
	    }

	typename graph_traits<Graph>::vertex_iterator v,v_end;
	for (tie(v,v_end) = vertices(g); v != v_end; ++v)
	{
	    size_t k = out_degree_no_loops(*v,g);
	    typename graph_traits<Graph>::out_edge_iterator e,e_end;
	    for(tie(e, e_end) = out_edges(*v,g); e != e_end; ++e)
	    {
		typename graph_traits<Graph>::vertex_descriptor t = target(*e,g);
		if (t != *v)
		{
		    if (get(s, t) == get(s, *v))
			modularity += get(weights, *e) - k*out_degree_no_loops(t,g)/double(2*E);
		}
	    }
	}
	modularity /= 2*W;
    }
};

double GraphInterface::GetModularity(string weight, string property)
{
    double modularity = 0;

    bool directed = _directed;
    _directed = false;
    try
    {
	dynamic_property_map& comm_prop = find_property_map(_properties, property, typeid(graph_traits<multigraph_t>::vertex_descriptor));
	DynamicPropertyMapWrap<size_t,graph_traits<multigraph_t>::vertex_descriptor> comm_map(comm_prop);
	if(weight != "")
	{
	    dynamic_property_map& weight_prop = find_property_map(_properties, weight, typeid(graph_traits<multigraph_t>::edge_descriptor));
	    DynamicPropertyMapWrap<double,graph_traits<multigraph_t>::edge_descriptor> weight_map(weight_prop);
	    check_filter(*this, bind<void>(get_modularity(), _1, var(weight_map), var(comm_map), var(modularity)),
			 reverse_check(), always_undirected());	    
	}
	else
	{
	    ConstantPropertyMap<double,graph_traits<multigraph_t>::edge_descriptor> weight_map(1.0);
	    check_filter(*this, bind<void>(get_modularity(), _1, var(weight_map), var(comm_map), var(modularity)),
			 reverse_check(), always_undirected());
	}
    }
    catch (property_not_found& e)
    {
	throw GraphException("error getting scalar property: " + string(e.what()));
    }
    _directed = directed;

    return modularity;
}
