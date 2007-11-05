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
#include <boost/functional/hash.hpp>

#include "graph.hh"
#include "graph_filtering.hh"
#include "shared_map.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

typedef boost::mt19937 rng_t;

struct get_in_degree
{
    template <class Graph>
    size_t operator()(typename graph_traits<Graph>::vertex_descriptor v, const Graph& g)
    {
        return get_degree(v, g, typename is_convertible<typename graph_traits<Graph>::directed_category,directed_tag>::type());
    }

    template <class Graph>
    size_t get_degree(typename graph_traits<Graph>::vertex_descriptor v, const Graph& g, true_type)
    {
        return in_degree(v,g);
    }

    template <class Graph>
    size_t get_degree(typename graph_traits<Graph>::vertex_descriptor v, const Graph& g, false_type)
    {
        return out_degree(v,g);
    }
};

struct source_in
{
    template <class Graph>
    typename graph_traits<Graph>::vertex_descriptor operator()(typename graph_traits<Graph>::edge_descriptor e, const Graph& g)
    {
        return get_source(e, g, typename is_convertible<typename graph_traits<Graph>::directed_category,directed_tag>::type());
    }

    template <class Graph>
    typename graph_traits<Graph>::vertex_descriptor get_source(typename graph_traits<Graph>::edge_descriptor e, const Graph& g, true_type)
    {
        return source(e, g);
    }

    template <class Graph>
    typename graph_traits<Graph>::vertex_descriptor get_source(typename graph_traits<Graph>::edge_descriptor e, const Graph& g, false_type)
    {
        return target(e, g);
    }
};

struct target_in
{
    template <class Graph>
    typename graph_traits<Graph>::vertex_descriptor operator()(typename graph_traits<Graph>::edge_descriptor e, const Graph& g)
    {
        return get_target(e, g, typename is_convertible<typename graph_traits<Graph>::directed_category,directed_tag>::type());
    }

    template <class Graph>
    typename graph_traits<Graph>::vertex_descriptor get_target(typename graph_traits<Graph>::edge_descriptor e, const Graph& g, true_type)
    {
        return target(e, g);
    }

    template <class Graph>
    typename graph_traits<Graph>::vertex_descriptor get_target(typename graph_traits<Graph>::edge_descriptor e, const Graph& g, false_type)
    {
        return source(e, g);
    }
};

template <class Graph>
struct get_in_edges
{
    typedef typename mpl::if_<typename is_convertible<typename graph_traits<Graph>::directed_category,directed_tag>::type,
                              typename graph_traits<Graph>::in_edge_iterator,
                              typename graph_traits<Graph>::out_edge_iterator>::type iterator;

    pair<iterator,iterator> operator()(typename graph_traits<Graph>::vertex_descriptor v, const Graph& g)
    {
        return get_edges(v, g, typename is_convertible<typename graph_traits<Graph>::directed_category,directed_tag>::type());
    }

    pair<iterator,iterator> get_edges(typename graph_traits<Graph>::vertex_descriptor v, const Graph& g, true_type)
    {
        return in_edges(v,g);
    }

    pair<iterator,iterator> get_edges(typename graph_traits<Graph>::vertex_descriptor v, const Graph& g, false_type)
    {
        return out_edges(v,g);
    }
};

template <class Graph>
bool is_adjacent(typename graph_traits<Graph>::vertex_descriptor u, typename graph_traits<Graph>::vertex_descriptor v, Graph& g )
{
    typename graph_traits<Graph>::out_edge_iterator e, e_end;
    for (tie(e, e_end) = out_edges(u, g); e != e_end; ++e)
    {
        if (target(*e,g) == v)
            return true;
    }
    return false;
}

template <class Graph, class EdgeIsNew>
bool is_adjacent_in_new(typename graph_traits<Graph>::vertex_descriptor u, typename graph_traits<Graph>::vertex_descriptor v, Graph& g , EdgeIsNew& edge_is_new,
                        const typename graph_traits<Graph>::edge_descriptor& e1, const typename graph_traits<Graph>::edge_descriptor& e2)
{
    typename graph_traits<Graph>::out_edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = out_edges(u, g); ei != ei_end; ++ei)
    {
        if ( (target(*ei, g) == v) && edge_is_new[*ei] && (*ei != e1) && (*ei != e2) )
            return true;
    }
    return false;
}

template <class EdgeIterator>
typename iterator_traits<EdgeIterator>::value_type sample_edge(const EdgeIterator& e_begin, const EdgeIterator& e_end, rng_t& rng)
{
    vector<typename iterator_traits<EdgeIterator>::value_type> candidates;
    EdgeIterator e;
    for (e = e_begin; e != e_end; ++e)
        candidates.push_back(*e);

    uniform_int<size_t> edge_sample(0, candidates.size() - 1);
    size_t edge = edge_sample(rng);
    return candidates[edge];
}

template <class Graph, class EdgeIndexMap, class EdgesType>
void swap_edge_triad( typename graph_traits<Graph>::edge_descriptor& e, typename graph_traits<Graph>::edge_descriptor& se,
                      typename graph_traits<Graph>::edge_descriptor& te, EdgesType& edges, EdgeIndexMap edge_index, Graph& g )
{
            typename graph_traits<Graph>::edge_descriptor ne, nse, nte;
//cout << "antes : e " << source(e,g) << ',' << target(e,g) << " se " << source(se,g) << ',' << target(se,g) << " te " << source_in()(te,g) << ',' << target_in()(te,g) << endl;
//cout << "edge index:" << " e " << edge_index[e] << " se " << edge_index[se] << " te " << edge_index[te] << endl;
            // split cases where different combinations of the three edges are the same
            if( se != te )
            {
                ne = add_edge(source(se, g), target_in()(te, g), g).first;
                if(e != se)
                {
                    nse = add_edge(source(e, g), target(se, g), g).first;
                    edge_index[nse] = edge_index[se];
                    remove_edge(se, g);
                    edges[edge_index[nse]] = nse;
se = nse;
                }
else { se = ne; }
                if(e != te)
                {
                    nte = add_edge(source_in()(te, g), target(e, g), g).first;
                    edge_index[nte] = edge_index[te];
                    remove_edge(te, g);
                    edges[edge_index[nte]] = nte;
te = nte;
                }
else { te = ne; }
                edge_index[ne] = edge_index[e];
                remove_edge(e, g);
                edges[edge_index[ne]] = ne;
e = ne;
            }
            else
            {
                if(e != se)
                {
                    ne = se; nse = e;
                    tie(edge_index[ne], edge_index[nse]) = make_pair(edge_index[e], edge_index[se]);
                    edges[edge_index[ne]] = ne; edges[edge_index[nse]] = nse;
e = ne; se = nse; te = nse;
                }
            }
//cout << "depois : e " << source(e,g) << ',' << target(e,g) << " se " << source(se,g) << ',' << target(se,g) << " te " << source(te,g) << ',' << target(te,g) << endl;
//cout << "edge index:" << " e " << edge_index[e] << " se " << edge_index[se] << " te " << edge_index[te] << endl;
}


template <template <class Graph, class EdgeIndexMap> class RewireStrategy>
struct graph_rewire
{
    template <class Graph, class EdgeIndexMap>
    void operator()(const Graph& og, EdgeIndexMap edge_index, size_t seed, bool self_loops, bool parallel_edges) const
    {
        Graph& g = const_cast<Graph&>(og);

        rng_t rng(static_cast<rng_t::result_type>(seed));

        if (!self_loops)
        {
            // check the existence of self-loops
            bool has_self_loops = false;
            int i, N = num_vertices(g);
            #pragma omp parallel for default(shared) private(i) schedule(dynamic)
            for (i = 0; i < N; ++i)
            {
                typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
                if (v == graph_traits<Graph>::null_vertex())
                    continue;
                if (is_adjacent(v, v, g))
                    has_self_loops = true;
            }
            if (has_self_loops)
                throw GraphException("Self-loop detected. Can't rewire graph without self-loops if it already contains self-loops!");
        }

        if (!parallel_edges)
        {
            // check the existence of parallel edges
            bool has_parallel_edges = false;
            int i, N = num_vertices(g);
            #pragma omp parallel for default(shared) private(i) schedule(dynamic)
            for (i = 0; i < N; ++i)
            {
                typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
                if (v == graph_traits<Graph>::null_vertex())
                    continue;

                tr1::unordered_set<typename graph_traits<Graph>::vertex_descriptor> targets;
                typename graph_traits<Graph>::out_edge_iterator e, e_end;
                for (tie(e, e_end) = out_edges(v, g); e != e_end; ++e)
                {
                    if (targets.find(target(*e, g)) != targets.end())
                        has_parallel_edges = true;
                    else
                        targets.insert(target(*e, g));
                }
            }

            if (has_parallel_edges)
                throw GraphException("Parallel edge detected. Can't rewire graph without parallel edges if it already contains parallel edges!");
        }

        RewireStrategy<Graph, EdgeIndexMap> rewire(g, rng);

        vector<typename graph_traits<Graph>::edge_descriptor> edges(num_edges(g));
        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            typename graph_traits<Graph>::out_edge_iterator e, e_end;
            for (tie(e, e_end) = out_edges(v, g); e != e_end; ++e)
                edges[edge_index[*e]] = *e;
        }

        // for each edge simultaneously rewire its source and target
        for (i = 0; i < int(edges.size()); ++i)
        {
            typename graph_traits<Graph>::edge_descriptor e = edges[i];
            typename graph_traits<Graph>::edge_descriptor se, te;
            tie(se, te) = rewire(e, self_loops, parallel_edges);
            swap_edge_triad(e, se, te, edges, edge_index, g);
        }
    }
};

template <class Graph>
class RandomRewireStrategy
{
public:
    RandomRewireStrategy (const Graph& g, rng_t& rng): _g(g), _rng(rng)
    {
        int i, N = num_vertices(_g);
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, _g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            size_t k = get_in_degree()(v, _g);
            if (k == 0)
                continue;
            if (!_vertices.empty())
                _vertices[_vertices.rbegin()->first + k] = v;
            else
                _vertices[k] = v;
        }
    }

    pair<typename graph_traits<Graph>::vertex_descriptor, typename graph_traits<Graph>::vertex_descriptor>
    operator()(const typename graph_traits<Graph>::edge_descriptor& e, bool self_loops)
    {
        uniform_int<size_t> vertex_sample(0, _vertices.rbegin()->first-1);
        size_t v1,v2;
        typename graph_traits<Graph>::vertex_descriptor s,t;
        do
        {
            v1 = vertex_sample(_rng);
            v2 = vertex_sample(_rng);
            s = _vertices.upper_bound(v1)->second;
            t = _vertices.upper_bound(v2)->second;
        }
        while (!self_loops && (s == t));
        return make_pair(s, t);
    }

private:
    const Graph& _g;
    rng_t& _rng;
    typedef map<size_t, typename graph_traits<Graph>::vertex_descriptor> vertices_t;
    vertices_t _vertices;
};

template <class Graph, class EdgeIndexMap>
class CorrelatedRewireStrategy
{
public:
    CorrelatedRewireStrategy (const Graph& g, rng_t& rng): _g(g), _rng(rng)
    {
        int i, N = num_vertices(_g);
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, _g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            size_t j = get_in_degree()(v, _g);
            size_t k = out_degree(v, _g);

            if (j > 0)
                _deg_target_vertices[make_pair(j,k)].push_back(v);
            
            if (k > 0)
                _deg_source_vertices[make_pair(j,k)].push_back(v);
        }
    }

    pair<typename graph_traits<Graph>::edge_descriptor, typename graph_traits<Graph>::edge_descriptor>
    operator()(const typename graph_traits<Graph>::edge_descriptor& e, bool self_loops, bool parallel_edges)
    {
        typename deg_vertices_t::value_type::second_type& source_vertices = _deg_source_vertices[make_pair(get_in_degree()(source(e, _g), _g),
                                                                                                 out_degree(source(e, _g), _g))];
        typename deg_vertices_t::value_type::second_type& target_vertices = _deg_target_vertices[make_pair(get_in_degree()(target(e, _g), _g),
                                                                                                 out_degree(target(e, _g), _g))];

        //shuffle the sources
        random_number_generator<rng_t, size_t> generator(_rng);
        random_shuffle(source_vertices.begin(), source_vertices.end(), generator);
        random_shuffle(target_vertices.begin(), target_vertices.end(), generator);

        //try new targets until one passes all the consistency checks
        typename graph_traits<Graph>::edge_descriptor se, te;
        typename deg_vertices_t::value_type::second_type::iterator vs, vs_begin = source_vertices.begin(), vs_end = source_vertices.end();
        typename deg_vertices_t::value_type::second_type::iterator vt, vt_begin = target_vertices.begin(), vt_end = target_vertices.end();
        bool pass=false;
        for(vs = vs_begin ; vs != vs_end; ++vs)
        {
            for(vt = vt_begin; vt != vt_end; ++vt)
            {
                if( (!self_loops && (*vs == *vt)) || (!parallel_edges && is_adjacent_in_new(*vs, *vt, _g, _edge_is_new, e, e)) )
                    continue;

                // set up the shuffled edges
                vector<typename graph_traits<Graph>::edge_descriptor> in_edges_vt, out_edges_vs;
                typename get_in_edges<Graph>::iterator ie, ie_end;
                typename graph_traits<Graph>::out_edge_iterator oe, oe_end;
                for( tie(ie, ie_end) = get_in_edges<Graph>()(*vt, _g); ie != ie_end ; ++ie)
                    in_edges_vt.push_back(*ie);
                for( tie(oe, oe_end) = out_edges(*vs, _g); oe != oe_end ; ++oe)
                    out_edges_vs.push_back(*oe);
                random_shuffle(in_edges_vt.begin(), in_edges_vt.end(), generator);
                random_shuffle(out_edges_vs.begin(), out_edges_vs.end(), generator);

                // for combinations of in_vt and out_vs...
                typename vector<typename graph_traits<Graph>::edge_descriptor>::iterator out_vs_i, out_vs_i_begin=out_edges_vs.begin(), out_vs_i_end=out_edges_vs.end();
                typename vector<typename graph_traits<Graph>::edge_descriptor>::iterator in_vt_i, in_vt_i_begin=in_edges_vt.begin(), in_vt_i_end=in_edges_vt.end();
                for(out_vs_i = out_vs_i_begin; out_vs_i != out_vs_i_end; ++out_vs_i)
                {
                    for(in_vt_i = in_vt_i_begin; in_vt_i != in_vt_i_end; ++in_vt_i)
                    {
                        se = *out_vs_i;
                        te = *in_vt_i;
                            if( ( !self_loops && (source_in()(*in_vt_i, _g) == target(e, _g))) ||
                                ( !parallel_edges &&
                                  is_adjacent_in_new(source_in()(*in_vt_i, _g), target(e, _g), _g, _edge_is_new, *out_vs_i, *in_vt_i) ) )
                                continue;
                            if( ( !self_loops && (source(e, _g) == target(*out_vs_i, _g))) ||
                                ( !parallel_edges &&
                                  is_adjacent_in_new(source(e, _g), target(*out_vs_i, _g), _g, _edge_is_new, *out_vs_i, *in_vt_i) ) )
                                break;
                            if( !parallel_edges &&
                                _edge_is_new[*in_vt_i] && ( (source_in()(*in_vt_i, _g)==*vs) && (target(e, _g)==*vt) ) )
                                continue;
                            if( !parallel_edges &&
                                _edge_is_new[*out_vs_i] && ( (source(e, _g)==*vs) && (target(*out_vs_i, _g)==*vt) ) )
                                break;
                            if( !parallel_edges && (_edge_is_new[*in_vt_i] && _edge_is_new[*out_vs_i]) &&
                                (source(e, _g)==source_in()(*in_vt_i, _g)) && (target(e, _g)==target(*out_vs_i, _g)) &&
                                (*in_vt_i != *out_vs_i) )
                                continue;
                        pass = true;
                        break;
                    }
                    if(pass) break;
                }
                if(pass) break;
            }
            if(pass) break;
        }
//cout << "pass " << pass << endl;
        _edge_is_new[e]=true;
        if(pass==false)
            throw GraphException("Bad things happen when you're not one with your inner self.");
        return make_pair(se, te);
    }

private:
    const Graph& _g;
    rng_t& _rng;
    typedef tr1::unordered_map<pair<size_t, size_t>, vector<typename graph_traits<Graph>::vertex_descriptor>, hash<pair<size_t, size_t> > > deg_vertices_t;
    deg_vertices_t _deg_source_vertices;
    deg_vertices_t _deg_target_vertices;
    vector_property_map<bool, EdgeIndexMap> _edge_is_new;
};

//==============================================================================
// RandomRewire
//==============================================================================
void GraphInterface::RandomRewire(rewire_strat_t strat, bool self_loops, bool parallel_edges, size_t seed)
{
    bool reversed = GetReversed();
    SetReversed(false);

/*
    if (strat == UNCORRELATED_STRAT)
        check_filter(*this, bind<void>(graph_rewire<RandomRewireStrategy>(), _1, _edge_index, seed, self_loops, parallel_edges),
                     never_reversed(), directed_check());
    else
*/
        check_filter(*this, bind<void>(graph_rewire<CorrelatedRewireStrategy>(), _1, _edge_index, seed, self_loops, parallel_edges),
                     never_reversed(), directed_check());

    SetReversed(reversed);
}
