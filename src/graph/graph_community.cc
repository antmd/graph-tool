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
// ManagedUnorderedMap
//==============================================================================
template <class Key, class Value>
class ManagedUnorderedMap: public tr1::unordered_map<Key,Value>
{
    typedef tr1::unordered_map<Key,Value> base_t;
public:
    void erase(typename base_t::iterator pos)
    {
	static_cast<base_t*>(this)->erase(pos);
	manage();
    }

    size_t erase(const Key& k)
    {
	size_t n = static_cast<base_t*>(this)->erase(k);
	manage();
	return n;
    }

    void manage()
    {
	if (this->bucket_count() > 2*this->size())
	{
	    base_t* new_map = new base_t;
	    for (typeof(this->begin()) iter = this->begin(); iter != this->end(); ++iter)
		(*new_map)[iter->first] = iter->second;
	    *static_cast<base_t*>(this) = *new_map;
	    delete new_map;
	}
    }

};


//==============================================================================
// GetCommunityStructure()
// computes the community structure through a spin glass system with 
// simulated annealing
//==============================================================================

template <template <class G, class CommunityMap> class NNKS>
struct get_communities
{
    template <class Graph, class WeightMap, class CommunityMap>
    void operator()(Graph& g, WeightMap weights, CommunityMap s, double gamma, size_t n_iter, double Tmin, double Tmax, size_t seed, bool verbose) const
    {
	rng_t rng(static_cast<rng_t::result_type>(seed));
	boost::uniform_real<double> uniform_p(0.0,1.0);

	ManagedUnorderedMap<size_t, size_t> Nk; // degree histogram
	ManagedUnorderedMap<size_t, size_t> Ns; // spin histogram
	ManagedUnorderedMap<size_t, map<double, unordered_set<size_t> > > global_term; // global energy term

	// init spins from [0,N-1] and global info	
	size_t index = 0;
	typename graph_traits<Graph>::vertex_iterator v,v_end;
	for (tie(v,v_end) = vertices(g); v != v_end; ++v)
	{
	    s[*v] = index++; 
	    Ns[s[*v]]++;
	    Nk[out_degree_no_loops(*v,g)]++;
	}

	NNKS<Graph,CommunityMap> Nnnks(g, s); // this will retrieve the expected number of neighbours with given spin, in funcion of degree

	for (typeof(Nk.begin()) iter = Nk.begin(); iter != Nk.end(); ++iter)
	    for (tie(v,v_end) = vertices(g); v != v_end; ++v)
		global_term[iter->first][gamma*Nnnks(iter->first,s[*v])].insert(s[*v]);

	// define cooling rate so that temperature starts at Tmax at temp_count == 0
	// and reaches Tmin at temp_count == n_iter - 1
	if (Tmin < numeric_limits<double>::epsilon())
	    Tmin = numeric_limits<double>::epsilon();
	double cooling_rate = -(log(Tmin)-log(Tmax))/(n_iter-1);

	for (size_t temp_count = 0; temp_count < n_iter; ++temp_count)
	{
	    double T = Tmax*exp(-cooling_rate*temp_count);
	    
            // sample a new spin for every vertex
	    for (tie(v,v_end) = vertices(g); v != v_end; ++v)
	    {
		ManagedUnorderedMap<size_t, double> ns; // number of neighbours with spin 's' (weighted)
		
		// neighbourhood spins info
		typename graph_traits<Graph>::out_edge_iterator e,e_end;
		for (tie(e,e_end) = out_edges(*v,g); e != e_end; ++e)
		{
		    typename graph_traits<Graph>::vertex_descriptor t = target(*e,g);
		    if (t != *v)
			ns[s[t]] += get(weights, *e);
		}
		    
		size_t k = out_degree_no_loops(*v,g);
		
		double E = -T*log(1-uniform_p(rng)) - (gamma+1)*k; // sampled target jump energy

		// local energy term
		double local_term = ns[s[*v]] - gamma*Nnnks(k,s[*v]);
		
		// update global info with local info
		for (typeof(ns.begin()) iter = ns.begin(); iter != ns.end(); ++iter)
		{
		    double old_E = gamma*Nnnks(k,iter->first);
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
		    double old_E = gamma*Nnnks(k,iter->first) - ns[iter->first];
		    global_term[k][old_E].erase(iter->first);
		    if (global_term[k][old_E].empty())
			global_term[k].erase(old_E);
		    global_term[k][old_E + ns[iter->first]].insert(iter->first);
		}
		
		//update global info
		if (s[*v] != a)
		{
		    double old_E = gamma*Nnnks(k,s[*v]);
		    global_term[k][old_E].erase(s[*v]);
		    if (global_term[k][old_E].empty())
			global_term[k].erase(old_E);
		    old_E = gamma*Nnnks(k,a);
		    global_term[k][old_E].erase(a);
		    if (global_term[k][old_E].empty())
			global_term[k].erase(old_E);
		    Nnnks.Update(*v,s[*v],a);
		    Ns[s[*v]]--;
		    if (Ns[s[*v]] == 0)
			Ns.erase(s[*v]);
		    Ns[a]++;
		    global_term[k][gamma*Nnnks(k,s[*v])].insert(s[*v]);
		    global_term[k][gamma*Nnnks(k,a)].insert(a);
		    
		    // update spin
		    s[*v] = a;
		}
	    }

	    if (verbose)
	    {
		static stringstream str;
		for (size_t j = 0; j < str.str().length(); ++j)
		    cout << "\b";
		str.str("");
		str << temp_count << " of " << n_iter << " (" << (temp_count+1)*100/n_iter << "%) " << "temperature: " << T << " spins: " << Ns.size() << " energy levels: ";
		size_t n_energy = 0;
		for (typeof(global_term.begin()) iter = global_term.begin(); iter != global_term.end(); ++iter)
		    n_energy += iter->second.size();
		str << n_energy;
		cout << str.str() << flush;
	    }
	}
	if (verbose)
	    cout << endl;
    }
};

template <class Graph, class CommunityMap>
class NNKSErdosReyni
{
public:
    NNKSErdosReyni(Graph &g, CommunityMap s)
    {
	size_t N = 0;
	double _avg_k = 0.0;
	typename graph_traits<Graph>::vertex_iterator v,v_end;
	for (tie(v,v_end) = vertices(g); v != v_end; ++v)
	{
	    size_t k = out_degree_no_loops(*v,g); 
	    _avg_k += k;
	    N++;
	    _Ns[s[*v]]++;
	}
	_p = _avg_k/(N*N);
    }

    void Update(typename graph_traits<Graph>::vertex_descriptor v, size_t old_s, size_t s)
    {
	_Ns[old_s]--;
	if (_Ns[old_s] == 0)
	    _Ns.erase(old_s);
	_Ns[s]++;
    }

    double operator()(size_t k, size_t s) const 
    {
	size_t ns = 0;
	typeof(_Ns.begin()) iter = _Ns.find(s);
	if (iter != _Ns.end())
	    ns = iter->second;
	return _p*ns;
    }

private:
    double _p;
    ManagedUnorderedMap<size_t,size_t> _Ns;
};

template <class Graph, class CommunityMap>
class NNKSUncorr
{
public:
    NNKSUncorr(Graph &g, CommunityMap s): _g(g), _K(0)
    {
	typename graph_traits<Graph>::vertex_iterator v,v_end;
	for (tie(v,v_end) = vertices(_g); v != v_end; ++v)
	{
	    size_t k = out_degree_no_loops(*v, _g);
	    _K += k;
	    _Ks[s[*v]] += k; 
	}
    }

    void Update(typename graph_traits<Graph>::vertex_descriptor v, size_t old_s, size_t s)
    {
	size_t k = out_degree_no_loops(v, _g);
	_Ks[old_s] -= k;
	if (_Ks[old_s] == 0)
	    _Ks.erase(old_s);
	_Ks[s] += k; 
    }


    double operator()(size_t k, size_t s) const 
    {
	size_t ks = 0;
	typeof(_Ks.begin()) iter = _Ks.find(s);
	if (iter != _Ks.end())
	    ks = iter->second;
	return k*ks/double(_K);
    }
private:
    Graph& _g;
    size_t _K;
    ManagedUnorderedMap<size_t, size_t> _Ks;
};

template <class Graph, class CommunityMap>
class NNKSCorr
{
public:
    NNKSCorr(Graph &g, CommunityMap s): _g(g)
    {
	typename graph_traits<Graph>::vertex_iterator v,v_end;
	for (tie(v,v_end) = vertices(_g); v != v_end; ++v)
	{
	    size_t k = out_degree_no_loops(*v, _g);
	    _Nk[k]++;
	    _kNs[k][s[*v]]++;
	    unordered_set<typename graph_traits<Graph>::vertex_descriptor> neighbours;
	    neighbours.insert(*v);
	    typename graph_traits<Graph>::adjacency_iterator a,a_end;
	    for (tie(a,a_end) = adjacent_vertices(*v, _g); a != a_end; ++a)
		if (neighbours.find(*a) == neighbours.end())
		{
		    _Nkk[k][out_degree_no_loops(*a, _g)]++;
		    neighbours.insert(*a);
		}
	}

	for (typeof(_Nk.begin()) k_iter = _Nk.begin(); k_iter != _Nk.end(); ++k_iter)
	    for (tie(v,v_end) = vertices(_g); v != v_end; ++v)
	    {
		size_t k = out_degree_no_loops(*v, _g);
		_kNNs[k_iter->first][s[*v]] += _Nkk[k_iter->first][k] * _kNs[k][s[*v]]/double(_Nk[k_iter->first]*_Nk[k]);
	    }
    }

    void Update(typename graph_traits<Graph>::vertex_descriptor v, size_t old_s, size_t s)
    {
	size_t k = out_degree_no_loops(v, _g);
	for (typeof(_Nk.begin()) k_iter = _Nk.begin(); k_iter != _Nk.end(); ++k_iter)
	{
	    _kNNs[k_iter->first][old_s] -= _Nkk[k_iter->first][k] * _kNs[k][old_s]/double(_Nk[k_iter->first]*_Nk[k]);
	    if (_kNNs[k_iter->first][old_s] == 0.0)
		_kNNs[k_iter->first].erase(old_s);
	    _kNNs[k_iter->first][s] += _Nkk[k_iter->first][k] * _kNs[k][s]/double(_Nk[k_iter->first]*_Nk[k]);
	}
	
    }

    double operator()(size_t k, size_t s)
    {
	double knns = 0.0;
	typeof(_kNNs[k].begin()) iter = _kNNs[k].find(s);
	if (iter != _kNNs[k].end())
	    knns = iter->second;
	return knns;
    }

private:
    Graph& _g;
    ManagedUnorderedMap<size_t, size_t> _Nk;
    ManagedUnorderedMap<size_t, ManagedUnorderedMap<size_t, size_t> > _Nkk;
    ManagedUnorderedMap<size_t, ManagedUnorderedMap<size_t,size_t> > _kNs;
    ManagedUnorderedMap<size_t, ManagedUnorderedMap<size_t, double> > _kNNs;
};


struct get_communities_selector
{
    get_communities_selector(GraphInterface::comm_corr_t corr):_corr(corr) {}
    GraphInterface::comm_corr_t _corr;

    template <class Graph, class WeightMap, class CommunityMap>
    void operator()(Graph& g, WeightMap weights, CommunityMap s, double gamma, size_t n_iter, double Tmin, double Tmax, size_t seed, bool verbose) const
    {
	switch (_corr)
	{
	case GraphInterface::ERDOS_REYNI:
	    get_communities<NNKSErdosReyni>()(g, weights, s, gamma, n_iter, Tmin, Tmax, seed, verbose);
	    break;
	case GraphInterface::UNCORRELATED:
	    get_communities<NNKSUncorr>()(g, weights, s, gamma, n_iter, Tmin, Tmax, seed, verbose);
	    break;
	case GraphInterface::CORRELATED:
	    get_communities<NNKSCorr>()(g, weights, s, gamma, n_iter, Tmin, Tmax, seed, verbose);
	    break;
	}
    }
};

void GraphInterface::GetCommunityStructure(double gamma, comm_corr_t corr, size_t n_iter, double Tmin, double Tmax, size_t seed, bool verbose, string weight, string property)
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
		check_filter(*this, bind<void>(get_communities_selector(corr), _1, var(weight_map), var(comm_map), gamma, n_iter, Tmin, Tmax, seed, verbose),
			     reverse_check(), always_undirected());
	    }
	    else
	    {
		DynamicPropertyMapWrap<double,graph_traits<multigraph_t>::edge_descriptor> weight_map(weight_prop);
		check_filter(*this, bind<void>(get_communities_selector(corr), _1, var(weight_map), var(comm_map), gamma, n_iter, Tmin, Tmax, seed, verbose),
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
	check_filter(*this, bind<void>(get_communities_selector(corr), _1, var(weight_map), var(comm_map), gamma, n_iter, Tmin, Tmax, seed, verbose),
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
