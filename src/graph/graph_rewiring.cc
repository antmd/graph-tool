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
// along with this program. If not, see <http://www.gnu.org/licenses/>.

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
    size_t operator()(typename graph_traits<Graph>::vertex_descriptor v,
                      const Graph& g)
    {
        return get_degree(v, g, typename is_convertible
                          <typename graph_traits<Graph>::directed_category,
                           directed_tag>::type());
    }

    template <class Graph>
    size_t get_degree(typename graph_traits<Graph>::vertex_descriptor v,
                      const Graph& g, true_type)
    {
        return in_degree(v,g);
    }

    template <class Graph>
    size_t get_degree(typename graph_traits<Graph>::vertex_descriptor v,
                      const Graph& g, false_type)
    {
        return out_degree(v,g);
    }
};

struct source_in
{
    template <class Graph>
    typename graph_traits<Graph>::vertex_descriptor
    operator()(typename graph_traits<Graph>::edge_descriptor e, const Graph& g)
    {
        return get_source(e, g, typename is_convertible
                          <typename graph_traits<Graph>::directed_category,
                           directed_tag>::type());
    }

    template <class Graph>
    typename graph_traits<Graph>::vertex_descriptor
    get_source(typename graph_traits<Graph>::edge_descriptor e, const Graph& g,
               true_type)
    {
        return source(e, g);
    }

    template <class Graph>
    typename graph_traits<Graph>::vertex_descriptor
    get_source(typename graph_traits<Graph>::edge_descriptor e, const Graph& g,
               false_type)
    {
        return target(e, g);
    }
};

struct target_in
{
    template <class Graph>
    typename graph_traits<Graph>::vertex_descriptor
    operator()(typename graph_traits<Graph>::edge_descriptor e, const Graph& g)
    {
        return get_target(e, g, typename is_convertible
                          <typename graph_traits<Graph>::directed_category,
                           directed_tag>::type());
    }

    template <class Graph>
    typename graph_traits<Graph>::vertex_descriptor
    get_target(typename graph_traits<Graph>::edge_descriptor e, const Graph& g,
               true_type)
    {
        return target(e, g);
    }

    template <class Graph>
    typename graph_traits<Graph>::vertex_descriptor
    get_target(typename graph_traits<Graph>::edge_descriptor e, const Graph& g,
               false_type)
    {
        return source(e, g);
    }
};

template <class Graph>
bool is_adjacent(typename graph_traits<Graph>::vertex_descriptor u,
                 typename graph_traits<Graph>::vertex_descriptor v, Graph& g )
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
bool is_adjacent_in_new(typename graph_traits<Graph>::vertex_descriptor u,
                        typename graph_traits<Graph>::vertex_descriptor v,
                        Graph& g , EdgeIsNew& edge_is_new,
                        const typename graph_traits<Graph>::edge_descriptor& e1,
                        const typename graph_traits<Graph>::edge_descriptor& e2)
{
    typename graph_traits<Graph>::out_edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = out_edges(u, g); ei != ei_end; ++ei)
    {
        if ( (target(*ei, g) == v) && edge_is_new[*ei] &&
             (*ei != e1) && (*ei != e2) )
            return true;
    }
    return false;
}

template <class RandomAccessIterator, class RandomNumberGenerator>
class iterative_random_shuffle
{
public:
    iterative_random_shuffle( RandomAccessIterator first,
                              RandomAccessIterator last,
                              RandomNumberGenerator& rand )
        : _i(first), _last(last), _rand(rand)
    {
        std::iter_swap(_i, _i + _rand(_last - _i) );
    }
    typename RandomAccessIterator::value_type operator*()
    {
        return *_i;
    }
    iterative_random_shuffle operator++()
    {
        ++_i;
        if(_i != _last)
            std::iter_swap(_i, _i + _rand(_last - _i) );
        return *this;
    }
    bool operator==(const RandomAccessIterator& i)
    {
        return _i == i;
    }
    bool operator!=(const RandomAccessIterator& i)
    {
        return _i != i;
    }
private:
    RandomAccessIterator _i, _last;
    RandomNumberGenerator& _rand;
};

template <class Graph, class EdgeIndexMap, class EdgesType>
void swap_edge_triad( typename graph_traits<Graph>::edge_descriptor& e,
                      typename graph_traits<Graph>::edge_descriptor& se,
                      typename graph_traits<Graph>::edge_descriptor& te,
                      EdgesType& edges, EdgeIndexMap edge_index, Graph& g )
{
            typename graph_traits<Graph>::edge_descriptor ne, nse, nte;
            // split cases where different combinations of the three edges are
            // the same
            if( se != te )
            {
                ne = add_edge(source(se, g), target_in()(te, g), g).first;
                if(e != se)
                {
                    nse = add_edge(source(e, g), target(se, g), g).first;
                    edge_index[nse] = edge_index[se];
                    remove_edge(se, g);
                    edges[edge_index[nse]] = nse;
                }
                if(e != te)
                {
                    nte = add_edge(source_in()(te, g), target(e, g), g).first;
                    edge_index[nte] = edge_index[te];
                    remove_edge(te, g);
                    edges[edge_index[nte]] = nte;
                }
                edge_index[ne] = edge_index[e];
                remove_edge(e, g);
                edges[edge_index[ne]] = ne;
            }
            else
            {
                if(e != se)
                {
                    ne = se; nse = e;
                    tie(edge_index[ne], edge_index[nse]) =
                        make_pair(edge_index[e], edge_index[se]);
                    edges[edge_index[ne]] = ne; edges[edge_index[nse]] = nse;
                }
            }
}

template <template <class Graph, class EdgeIndexMap> class RewireStrategy>
struct graph_rewire
{
    template <class Graph, class EdgeIndexMap>
    void operator()(Graph& g, EdgeIndexMap edge_index, size_t seed,
                    bool self_loops, bool parallel_edges) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;

        rng_t rng(static_cast<rng_t::result_type>(seed));

        if (!self_loops)
        {
            // check the existence of self-loops
            bool has_self_loops = false;
            int i, N = num_vertices(g);
            #pragma omp parallel for default(shared) private(i) \
                schedule(dynamic)
            for (i = 0; i < N; ++i)
            {
                vertex_t v = vertex(i, g);
                if (v == graph_traits<Graph>::null_vertex())
                    continue;
                if (is_adjacent(v, v, g))
                    has_self_loops = true;
            }
            if (has_self_loops)
                throw GraphException("Self-loop detected. Can't rewire graph "
                                     "without self-loops if it already contains"
                                     " self-loops!");
        }

        if (!parallel_edges)
        {
            // check the existence of parallel edges
            bool has_parallel_edges = false;
            int i, N = num_vertices(g);
            #pragma omp parallel for default(shared) private(i) \
                schedule(dynamic)
            for (i = 0; i < N; ++i)
            {
                vertex_t v = vertex(i, g);
                if (v == graph_traits<Graph>::null_vertex())
                    continue;

                tr1::unordered_set<vertex_t> targets;
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
                throw GraphException("Parallel edge detected. Can't rewire "
                                     "graph without parallel edges if it "
                                     "already contains parallel edges!");
        }

        RewireStrategy<Graph, EdgeIndexMap> rewire(g, rng);

        vector<edge_t> edges(num_edges(g));
        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            vertex_t v = vertex(i, g);
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
            tie(se, te) = rewire(e, edges, self_loops, parallel_edges);
            swap_edge_triad(e, se, te, edges, edge_index, g);
        }
    }
};

template <class Graph, class EdgeIndexMap>
class RandomRewireStrategy
{
public:
    RandomRewireStrategy (const Graph& g, rng_t& rng): _g(g), _rng(rng)
    {
    }
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    typedef typename graph_traits<Graph>::edge_descriptor edge_t;

    template<class EdgesType>
    pair<edge_t,edge_t> operator()(const edge_t& e, const EdgesType& edges,
                                   bool self_loops, bool parallel_edges)
    {
        _edges_source = edges;
        _edges_target = edges;
        random_number_generator<rng_t, size_t> generator(_rng);
        bool pass=false;
        edge_t se, te;

        //try randomly drawn pairs of edges until one satisfies all the
        //consistency checks
        iterative_random_shuffle<typename edges_t::iterator,
                                 random_number_generator<rng_t, size_t> >
            esi(_edges_source.begin(), _edges_source.end(), generator);
        for(; esi != _edges_source.end() ; ++esi)
        {
            iterative_random_shuffle<typename edges_t::iterator,
                                     random_number_generator<rng_t, size_t> >
                eti(_edges_target.begin(), _edges_target.end(), generator);
            for(; eti != _edges_target.end() ; ++eti)
            {
                if(!self_loops)
                {
                    if((source(e, _g) == target(*esi, _g)))
                        break;
                    if((source(*esi, _g) == target_in()(*eti, _g)) ||
                       (source_in()(*eti, _g) == target(e, _g)))
                        continue;
                }
                if(!parallel_edges)
                {
                    //TODO: make this readable - count0
                    if(is_adjacent_in_new(source(*esi, _g),
                                          target_in()(*eti, _g),
                                          _g, _edge_is_new, *esi, *eti) ||
                        is_adjacent_in_new(source_in()(*eti, _g),
                                           target(e, _g), _g,
                                           _edge_is_new, *esi, *eti) ||
                        is_adjacent_in_new(source(e, _g), target(*esi, _g),
                                           _g, _edge_is_new, *esi, *eti) ||
                       (_edge_is_new[*eti] &&
                        (source_in()(*eti, _g) == source(*esi, _g)) &&
                        (target( e  , _g) == target_in()(*eti, _g)) ) ||
                       (_edge_is_new[*esi] && (source( e  , _g) ==
                                               source(*esi, _g))
                        && (target(*esi, _g) == target_in()(*eti, _g)) ) ||
                       ( _edge_is_new[*eti] && _edge_is_new[*esi]  &&
                         (*eti != *esi) &&
                         (source(e, _g) == source_in()(*eti, _g)) &&
                         (target(e, _g) == target(*esi, _g)) ) )
                        continue;
                }
                se = *esi;
                te = *eti;
                pass = true;
                break;
            }
            if(pass) break;
        }
        if(pass==false)
            throw GraphException("Bad things happen when you're not one with"
                                 "your inner self."); // also when you don't
                                                      // comment your code
                                                      // properly :-) - count0
        _edge_is_new[e]=true;
        return make_pair(se, te);
    }

private:
    const Graph& _g;
    rng_t& _rng;
    typedef vector<edge_t> edges_t;
    edges_t _edges_target, _edges_source;
    vector_property_map<bool, EdgeIndexMap> _edge_is_new;
};

template <class Graph, class EdgeIndexMap>
class CorrelatedRewireStrategy
{
public:
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    typedef typename graph_traits<Graph>::edge_descriptor edge_t;

    CorrelatedRewireStrategy (const Graph& g, rng_t& rng): _g(g), _rng(rng)
    {
        int i, N = num_vertices(_g);
        for (i = 0; i < N; ++i)
        {
            vertex_t v = vertex(i, _g);
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

    template<class EdgesType>
    pair<edge_t, edge_t> operator()(const edge_t& e, const EdgesType& edges,
                                    bool self_loops, bool parallel_edges)
    {
        vector<vertex_t>& source_vertices =
            _deg_source_vertices[make_pair(get_in_degree()(source(e, _g), _g),
                                           out_degree(source(e, _g), _g))];
        vector<vertex_t>& target_vertices =
            _deg_target_vertices[make_pair(get_in_degree()(target(e, _g), _g),
                                           out_degree(target(e, _g), _g))];

        random_number_generator<rng_t, size_t> generator(_rng);
        edge_t se, te;
        bool pass=false;

        //try new combinations until one satisfies all the consistency checks
        typedef iterative_random_shuffle
            <typename vector<vertex_t>::iterator,
             random_number_generator<rng_t, size_t> > viter_random_shuffle_t;
        viter_random_shuffle_t vs(source_vertices.begin(), 
                                  source_vertices.end(), generator);
        for(; vs != source_vertices.end(); ++vs)
        {
            viter_random_shuffle_t vt(target_vertices.begin(),
                                     target_vertices.end(), generator);
            for(; vt != target_vertices.end(); ++vt)
            {
                // set up the edges to draw from
                vector<edge_t> in_edges_vt, out_edges_vs;
                typename in_edge_iteratorS<Graph>::type ie, ie_end;
                typename graph_traits<Graph>::out_edge_iterator oe, oe_end;
                for(tie(ie, ie_end) = 
                        in_edge_iteratorS<Graph>::in_edges(*vt, _g);
                     ie != ie_end ; ++ie)
                    in_edges_vt.push_back(*ie);
                for(tie(oe, oe_end) = out_edges(*vs, _g); oe != oe_end ; ++oe)
                    out_edges_vs.push_back(*oe);

                // for combinations of in_vt and out_vs...
                typedef iterative_random_shuffle
                    <typename vector<edge_t>::iterator,
                     random_number_generator<rng_t, size_t> >
                    eiter_random_shuffle_t;
                eiter_random_shuffle_t out_vs_i(out_edges_vs.begin(),
                                                out_edges_vs.end(), generator);
                for(; out_vs_i != out_edges_vs.end(); ++out_vs_i)
                {
                    eiter_random_shuffle_t in_vt_i(in_edges_vt.begin(),
                                                   in_edges_vt.end(),
                                                   generator);
                    for(; in_vt_i != in_edges_vt.end(); ++in_vt_i)
                    {
                        if(!parallel_edges)  // check all break conditions
                                             // before continue conditions
                        {
                            if(_edge_is_new[*out_vs_i] &&
                               (source(e, _g)==*vs) &&
                               (target(*out_vs_i, _g)==*vt))
                                break;
                        }
                        if(!self_loops)
                        {
                            if((*vs == *vt) ||
                                (source(e, _g) == target(*out_vs_i, _g)))
                                break;
                            if((source_in()(*in_vt_i, _g) == target(e, _g)))
                                continue;
                        }
                        if(!parallel_edges)
                        {
                            if(is_adjacent_in_new(*vs, *vt, _g, _edge_is_new,
                                                  *out_vs_i, *in_vt_i) ||
                               is_adjacent_in_new(source_in()(*in_vt_i, _g),
                                                  target( e, _g), _g,
                                                  _edge_is_new, *out_vs_i,
                                                  *in_vt_i) ||
                               is_adjacent_in_new(source(e, _g),
                                                  target(*out_vs_i, _g), _g,
                                                  _edge_is_new, *out_vs_i,
                                                  *in_vt_i) ||
                               (_edge_is_new[*in_vt_i]
                                && (source_in()(*in_vt_i, _g) == *vs) &&
                                (target(e, _g)==*vt) ) ||
                               ( _edge_is_new[*in_vt_i] &&
                                 _edge_is_new[*out_vs_i] &&
                                 (*in_vt_i != *out_vs_i) &&
                                 (source(e, _g) ==
                                  source_in()(*in_vt_i, _g)) &&
                                 (target(e, _g) == target(*out_vs_i, _g))))
                                continue;
                        }
                        se = *out_vs_i;
                        te = *in_vt_i;
                        pass = true;
                        break;
                    }
                    if(pass) break;
                }
                if(pass) break;
            }
            if(pass) break;
        }
        if(pass==false)
            throw GraphException("Bad things happen when you're not one with"
                                 " your inner self.");
        _edge_is_new[e]=true;
        return make_pair(se, te);
    }

private:
    const Graph& _g;
    rng_t& _rng;
    typedef tr1::unordered_map<pair<size_t, size_t>, vector<vertex_t>,
                               hash<pair<size_t, size_t> > > deg_vertices_t;
    deg_vertices_t _deg_source_vertices;
    deg_vertices_t _deg_target_vertices;
    vector_property_map<bool, EdgeIndexMap> _edge_is_new;
};

void GraphInterface::RandomRewire(string strat, bool self_loops,
                                  bool parallel_edges, size_t seed)
{
    bool reversed = GetReversed();
    SetReversed(false);

    if (strat == "uncorrelated")
        check_filter(*this, bind<void>(graph_rewire<RandomRewireStrategy>(),
                                       _1, _edge_index, seed, self_loops,
                                       parallel_edges),
                     never_reversed(), directed_check());
    else if (strat == "correlated")
        check_filter(*this, bind<void>(graph_rewire<CorrelatedRewireStrategy>(),
                                       _1, _edge_index, seed, self_loops,
                                       parallel_edges),
                     never_reversed(), directed_check());
    else
        throw GraphException("invalid random rewire stategy: " + strat);
    SetReversed(reversed);
}
