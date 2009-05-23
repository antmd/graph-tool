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

#ifndef GRAPH_GENERATION_HH
#define GRAPH_GENERATION_HH

#include <tr1/unordered_set>
#include <tr1/unordered_map>
#include <boost/array.hpp>
#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
#include <boost/functional/hash.hpp>

#include <tr1/random>
#include <iostream>

namespace graph_tool
{
using namespace std;
using namespace boost;

typedef tr1::mt19937 rng_t;

// Graph Generation
// ================
//
// this graph generation routine is based on two inputs, the in-out degree
// probability matrix deg_matrix:
//
//   deg_matrix(j,k) = P(j,k)
//
// and the 4-dimensional in-out degree correlation matrix deg_corr:
//
//   deg_corr(j, k, j', k') = P((j,k) -> (j',k')), i.e. the probability that a
//   directed edge exists between two _specific_ vertices of degrees (j,k) and
//   (j',k')
//
// Both these structures are expensive to store, and having them would still
// mean that rejection sampling would have to be used to generate the graph,
// which can be quite slow, and unnecessary when the samples are easy to
// obtain. Therefore what is actually used are two functions, deg_sample() and
// deg_corr_sample(j,k) which should return a (j,k) degree sample which obeys
// the above probabilities.
//
// Furthermore we want also to generate _undirected_ graphs, which have the
// above probabilities changed to P(k) and deg_corr(k) respectively.

// Utilities
// =========

// desired vertex type, with desired j,k values and the index in the real graph
struct dvertex_t
{
    dvertex_t() {}
    dvertex_t(size_t in, size_t out): in_degree(in), out_degree(out) {}
    dvertex_t(const pair<size_t,size_t>& deg): in_degree(deg.first),
                                               out_degree(deg.second) {}
    size_t index, in_degree, out_degree;
    bool operator==(const dvertex_t& other) const {return other.index == index;}
};

// utility class to sample uniformly from a collection of values
template <class ValueType>
class Sampler
{
public:
    Sampler() {}

    template <class Iterator>
    Sampler(Iterator iter, Iterator end)
    {
        for(; iter != end; ++iter)
            _candidates.push_back(*iter);
    }

    void Insert(const ValueType& v) { _candidates.push_back(v); }

    void Remove(size_t index)
    {
        swap(_candidates[index], _candidates.back());
        _candidates.pop_back();
    }

    void RemoveLast()
    {
        if (_last != _candidates.size())
        {
            swap(_candidates[_last], _candidates.back());
            _candidates.pop_back();
        }
    }

    ValueType operator()(rng_t& rng, bool remove = true)
    {
        tr1::uniform_int<> sample(0, _candidates.size() - 1);
        int i = sample(rng);
        if (remove)
        {
            swap(_candidates[i], _candidates.back());
            ValueType ret = _candidates.back();
            _candidates.pop_back();
            _last = numeric_limits<size_t>::max();
            return ret;
        }
        else
        {
            _last = i;
            return _candidates[i];
        }
    }

private:
    vector<ValueType> _candidates;
    size_t _last;
};

// used for verbose display
void print_progress(size_t current, size_t total, stringstream& str)
{
    size_t atom = max(total/100, 1);
    if ( (current % atom == 0) || current == total)
    {
        for (size_t j = 0; j < str.str().length(); ++j)
            cout << "\b";
        str.str("");
        str << current+1 << " of " << total << " (" <<
            (current+1)*100/total << "%)";
        cout << str.str() << flush;
    }
}

void print_update(size_t current, stringstream& str)
{
    for (size_t j = 0; j < str.str().length(); ++j)
        cout << "\b";
    for (size_t j = 0; j < str.str().length(); ++j)
        cout << " ";
    for (size_t j = 0; j < str.str().length(); ++j)
        cout << "\b";
    str.str("");
    str << current;
    cout << str.str() << flush;
}

//
// Generation strategies
// =====================
//
// Directed and undirected graphs have different generation strategies, which
// are defined in the classes below.

//
// Directed graph generation strategy
//

class DirectedStrat
{
public:
    typedef pair<size_t, size_t> deg_t; // degree type

    DirectedStrat(size_t N, bool no_parallel, bool no_self_loops)
        : _N(N), _no_parallel(no_parallel), _no_self_loops(no_self_loops)
    {
        if (_no_parallel)
            _max_deg = (_no_self_loops) ? _N - 1 : N;
        else
            _max_deg = numeric_limits<size_t>::max();
    }

    // sample the degrees of all the vertices
    template <class Graph, class DegSample, class VertexIndex>
    size_t SampleDegrees(Graph& g, vector<dvertex_t>& vertices,
                         VertexIndex vertex_index, DegSample& deg_sample,
                         rng_t& rng, bool verbose)
    {
        stringstream str;
        size_t sum_j=0, sum_k=0;
        for(size_t i = 0; i < _N; ++i)
        {
            if (verbose)
                print_progress(i, _N, str);
            dvertex_t& v = vertices[i];
            v.index = vertex_index[add_vertex(g)];
            do
            {
                tie(v.in_degree, v.out_degree) = deg_sample();
            }
            while (_no_parallel &&
                   (v.in_degree > _max_deg || v.out_degree > _max_deg));
            sum_j += v.in_degree;
            sum_k += v.out_degree;
        }

        if (verbose)
        {
            cout << "\nfixing average degrees. Total degree difference: "
                 << flush;
            str.str("");
        }

        // <j> and <k> must be the same. Re-sample random pairs until this holds
        tr1::uniform_int<size_t> vertex_sample(0, _N-1);
        size_t count = 0;
        while (sum_j != sum_k)
        {
            dvertex_t& v = vertices[vertex_sample(rng)];
            sum_j -= v.in_degree;
            sum_k -= v.out_degree;
            do
            {
                tie(v.in_degree, v.out_degree) = deg_sample();
            }
            while (_no_parallel &&
                   (v.in_degree > _max_deg || v.out_degree > _max_deg));
            sum_j += v.in_degree;
            sum_k += v.out_degree;
            if (verbose && (count % 100 == 0 || sum_j == sum_k))
                print_update(min(sum_j-sum_k, sum_k-sum_j), str);
            count++;
        }
        return sum_k;
    }

    deg_t GetDegree(const dvertex_t& v) { return make_pair(v.out_degree,
                                                           v.out_degree); }
    bool IsTarget(const dvertex_t& v) { return v.in_degree > 0; }

    template <class Graph>
    bool IsStillTarget(Graph& g, const dvertex_t& v)
    {
        return (in_degree(vertex(v.index, g), g) != v.in_degree);
    }

private:
    size_t _N;
    bool _no_parallel;
    bool _no_self_loops;
    size_t _max_deg;
};

//
// Undirected graph generation strategy
//

class UndirectedStrat
{
public:
    typedef size_t deg_t;

    UndirectedStrat(size_t N, bool no_parallel, bool no_self_loops)
        : _N(N), _no_parallel(no_parallel), _no_self_loops(no_self_loops)
    {
        if (_no_parallel)
            _max_deg = (_no_self_loops) ? _N - 1 : N;
        else
            _max_deg = numeric_limits<size_t>::max();
    }

    // samples the degress of all vertices
    template <class Graph, class DegSample, class VertexIndex>
    size_t SampleDegrees(Graph& g, vector<dvertex_t>& vertices,
                         VertexIndex vertex_index, DegSample& deg_sample,
                         rng_t& rng, bool verbose)
    {
        stringstream str;
        size_t sum_k=0;
        for(size_t i = 0; i < _N; ++i)
        {
            if (verbose)
                print_progress(i, _N, str);
            dvertex_t& v = vertices[i];
            v.index = vertex_index[add_vertex(g)];
            do
            {
                v.out_degree = deg_sample(true);
            }
            while (_no_parallel && v.out_degree > _max_deg);
            sum_k += v.out_degree;
        }

        if (verbose)
        {
            cout << "\nFixing degree sum (must be even): "
                 << flush;
            str.str("");
        }

        // sum_k must be an even number (2*num_edges). Re-sample degrees until
        // this holds
        tr1::uniform_int<size_t> vertex_sample(0, _N-1);
        size_t count = 0;
        while (sum_k % 2 != 0)
        {
            dvertex_t& v = vertices[vertex_sample(rng)];
            sum_k -= v.out_degree;
            do
            {
                v.out_degree = deg_sample(true);
            }
            while (_no_parallel && (v.out_degree > _max_deg));
            sum_k +=  v.out_degree;
            if (verbose && (count % 100 || sum_k % 2 == 0))
                print_update(sum_k, str);
            count++;
        }
        return sum_k;
    }

    deg_t GetDegree(const dvertex_t& v) { return v.out_degree; }
    bool IsTarget(const dvertex_t& v) { return v.out_degree > 0; }

    template <class Graph>
    bool IsStillTarget(Graph& g, const dvertex_t& v)
    {
        return (out_degree(vertex(v.index, g), g) != v.out_degree);
    }

private:
    size_t _N;
    bool _no_parallel;
    bool _no_self_loops;
    size_t _max_deg;
};

//
// Correlation strategies
// ======================
//
// Graphs can be uncorrelated or two-point correlated. The first is a special
// case of the second, but it happens quite often and it can be optimized. Below
// are two classes which implement both strategies.

//
// Correlated graph strategy
//

template <class Deg>
class CorrelatedStrat
{
public:
    void InsertTarget(const dvertex_t& v, const Deg& deg)
    {
        _targets.insert(make_pair(deg, v));
    }

    template <class CorrSample>
    dvertex_t GetRandomTarget(const Deg& source_deg, CorrSample& corr_sample,
                              rng_t& rng)
    {
        // keep sampling until we find a vertex
        typename targets_t::iterator iter,end;
        Deg target_deg;
        do
        {
            //choose the target vertex according to correlation
            target_deg = corr_sample(source_deg);
            tie(iter, end) = _targets.equal_range(target_deg);
        }
        while (iter == end);
        // sample random target on the same class
        Sampler<typename targets_t::iterator> sampler;
        for(; iter != end; ++iter)
            sampler.Insert(iter);
        _last = sampler(rng);
        return _last->second;
    }

    void RemoveLast() { _targets.erase(_last); }

private:
    typedef tr1::unordered_multimap<Deg, dvertex_t, hash<Deg> > targets_t;
    targets_t _targets;
    typename targets_t::iterator _last;
};

//
// Uncorrelated graph strategy
//

template <class Deg>
class UncorrelatedStrat
{
public:
    void InsertTarget(const dvertex_t& v, const Deg& deg)
    {
        _targets.Insert(v);
    }

    template <class CorrSample>
    dvertex_t GetRandomTarget(const Deg& source_deg, CorrSample& corr_sample,
                              rng_t& rng)
    {
        return _targets(rng, false);
    }

    void RemoveLast() { _targets.RemoveLast(); }

private:
    Sampler<dvertex_t> _targets;
};

//
// Main Algorithm
// ==============
//
// generates a directed or undirected graph with given degree distribution and
// correlation as defined above

template <class IsCorrelated>
struct gen_random_graph
{
    gen_random_graph(size_t N): _N(N) {}

    template <class Graph, class DegSample, class CorrDegSample>
    void operator()(Graph& g, DegSample& deg_sample,
                    CorrDegSample& deg_corr_sample, bool no_parallel,
                    bool no_self_loops, bool undirected,
                    size_t seed, bool verbose)
        const
    {
        size_t N = _N;
        typename property_map<Graph,vertex_index_t>::type vertex_index =
            get(vertex_index_t(), g);
        rng_t rng(static_cast<rng_t::result_type>(seed));

        // figure out the necessary strategies
        typedef typename is_convertible
            <typename graph_traits<Graph>::directed_category,
            directed_tag>::type is_directed;
        typedef typename mpl::if_<is_directed, DirectedStrat,
                                  UndirectedStrat>::type gen_strat_t;
        typedef typename gen_strat_t::deg_t deg_t;
        typedef typename mpl::if_<IsCorrelated, CorrelatedStrat<deg_t>,
                                  UncorrelatedStrat<deg_t> >::type corr_strat_t;

        gen_strat_t gen_strat(N, no_parallel, no_self_loops);
        corr_strat_t corr_strat;

        stringstream str; // used for verbose status

        if (verbose)
            cout << "adding vertices: " << flush;

        // sample the N (j,k) pairs
        vector<dvertex_t> vertices(N);
        size_t E = gen_strat.SampleDegrees(g, vertices, vertex_index,
                                           deg_sample, rng, verbose);

        vector<dvertex_t> sources; // edge sources

        // fill up sources, targets and target_degrees
        sources.reserve(E);
        for(size_t i = 0; i < N; ++i)
        {
            for(size_t k = 0; k < vertices[i].out_degree; ++k)
                sources.push_back(vertices[i]);
            if (gen_strat.IsTarget(vertices[i]))
                corr_strat.InsertTarget(vertices[i],
                                        gen_strat.GetDegree(vertices[i]));
        }

        // shuffle sources
        for (size_t i = 0; i < sources.size(); ++i)
        {
            tr1::uniform_int<size_t> source_sample(i, sources.size()-1);
            swap(sources[i], sources[source_sample(rng)]);
        }

        if (verbose)
        {
            cout << "\nadding edges: " << flush;
            str.str("");
        }

        // connect the sources to targets
        for (int i = 0; i < int(sources.size()); ++i)
        {
            dvertex_t source = sources[i], target;
            deg_t source_deg = gen_strat.GetDegree(source);

            // in undirected graphs, there's no difference between source and
            // target
            if (!is_directed::value && !gen_strat.IsStillTarget(g, source))
                continue;

            tr1::unordered_set<size_t> old_targets;
            if (no_parallel)
            {
                typename graph_traits<Graph>::out_edge_iterator e, e_end;
                for (tie(e, e_end) = out_edges(vertex(source.index, g), g);
                     e != e_end; ++e)
                    old_targets.insert(boost::target(*e,g));
            }

            // get target
            target = corr_strat.GetRandomTarget(source_deg, deg_corr_sample,
                                                rng);

            // if unacceptable (because parallel and/or self_loop), throw it
            // back and disconnect a random previous (different) source
            if ((no_parallel && (old_targets.find(target.index)
                                 != old_targets.end())) ||
                (no_self_loops && (source.index == target.index)))
            {
                if (i == 0) // just try again
                {
                    i--;
                    continue;
                }

                // We need to check if all previous sources are not all equal to
                // the first one. In that case, we should just discard this one
                // and keep going.
                int j;
                for (j = i-1; j >= 0; --j)
                    if (!(sources[j] == source))
                        break;
                if (j < 0)
                {
                    i--;
                    continue; // just try again
                }

                size_t s_index;
                tr1::uniform_int<size_t> source_sample(0, j);
                // keep trying: we don't want the same source again, and no
                // sources with zero out-degree (which can happen when the graph
                // is undirected, when they're orphaned after a removed edge)
                do
                {
                    s_index = source_sample(rng);
                }
                while (sources[s_index] == source ||
                       out_degree(vertex(sources[s_index].index, g), g) == 0);

                dvertex_t rsource = sources[s_index]; // random source

                // get the random edge
                typename graph_traits<Graph>::out_edge_iterator e, e_end;
                tie(e, e_end) = out_edges(vertex(rsource.index, g), g);
                Sampler<typename graph_traits<Graph>::edge_descriptor>
                    esampler(e, e_end);
                typename graph_traits<Graph>::edge_descriptor re =
                    esampler(rng, false); // random edge

                // if edge's target was already full, put it back in target list
                dvertex_t rtarget =
                    vertices[vertex_index[boost::target(re, g)]];
                if (!gen_strat.IsStillTarget(g, rtarget))
                    corr_strat.InsertTarget(rtarget,
                                            gen_strat.GetDegree(rtarget));

                // for undirected graphs, we need also to check the source
                if (!is_directed::value)
                {
                    dvertex_t rsource =
                        vertices[vertex_index[boost::source(re, g)]];
                    if (!gen_strat.IsStillTarget(g, rsource))
                        corr_strat.InsertTarget(rsource,
                                                gen_strat.GetDegree(rsource));
                }

                // remove and swap with previous source and continue from then
                remove_edge(re, g);
                swap(sources[i-1], sources[s_index]);
                i -= 2;
                continue;
            }

            //add edge
            typename graph_traits<Graph>::edge_descriptor e;
            e = add_edge(vertex(source.index, g),
                         vertex(target.index, g), g).first;

            // if target received all the edges it should, remove it from target
            // list
            if (!gen_strat.IsStillTarget(g, target))
                corr_strat.RemoveLast();

            if (verbose)
                print_progress(i, E, str);
        }

        if (verbose)
            cout << "\n";
    }

    size_t _N;
};

} // graph_tool namespace

#endif // GRAPH_GENERATION_HH
