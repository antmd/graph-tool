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

#include "graph_util.hh"

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
// What is actually used are two functions, deg_sample() and
// corr_prob((j1,k1),(j2,k2)) which should return a (j,k) degree sample which
// obeys the above probabilities, and the given correlation probability.
//
// Furthermore we want also to generate _undirected_ graphs, which have the
// above probabilities changed to P(k) and deg_corr(k) respectively.

// Utilities
// =========

// desired vertex type, with desired j,k values and the index in the real graph
struct dvertex_t
{
    dvertex_t(): in_degree(0), out_degree(0) {}
    dvertex_t(size_t in, size_t out): in_degree(in), out_degree(out) {}
    dvertex_t(const pair<size_t,size_t>& deg): in_degree(deg.first),
                                               out_degree(deg.second) {}
    size_t index, in_degree, out_degree;
    bool operator==(const dvertex_t& other) const {return other.index == index;}
};

size_t hash_value(const dvertex_t& v)
{
    return v.index;
}


// utility class to sample uniformly from a collection of values
template <class ValueType>
class Sampler
{
public:
    Sampler(bool biased=false): _biased(biased), _erased_prob(0) {}

    template <class Iterator>
    Sampler(Iterator iter, Iterator end):
        _biased(false)
    {
        for(; iter != end; ++iter)
        {
            _candidates.push_back(*iter);
            _candidates_set.insert(make_pair(*iter, _candidates.size()-1));
        }
        assert(!_candidates.empty());
    }

    void Insert(const ValueType& v, double p = 0.0)
    {
        _candidates.push_back(v);
        _candidates_set.insert(make_pair(v, _candidates.size()-1));
        if (_biased)
        {
            if (_probs.size() > 0)
                _probs.push_back(_probs.back()+p);
            else
                _probs.push_back(p);
            _erased.push_back(false);
        }
    }

    void Remove(const ValueType& v)
    {
        typeof(_candidates_set.begin()) iter, end, temp;
        tie(iter, end) = _candidates_set.equal_range(v);
        assert(iter != end);

        if (_biased)
        {
            while(_erased[iter->second])
            {
                temp = iter;
                iter++;
                _candidates_set.erase(temp);
                if(iter == end)
                    break;
            }

            size_t index = iter->second;
            _erased[index] = true;
            _erased_prob += (index > 0) ?
                _probs[index]-_probs[index-1] : _probs[index];
        }
        else
        {
            size_t index = iter->second;
            temp = _candidates_set.find(_candidates.back());
            swap(_candidates[index], _candidates.back());
            _candidates.pop_back();
            if (!_candidates.empty() && temp != iter)
            {
                _candidates_set.erase(temp);
                _candidates_set.insert(make_pair(_candidates[index], index));
            }
        }
        _candidates_set.erase(iter);

        clean();
    }

    ValueType operator()(rng_t& rng, bool remove = false)
    {
        if (!_biased)
        {
            tr1::uniform_int<> sample(0, _candidates.size() - 1);
            int i = sample(rng);
            if (remove)
            {
                swap(_candidates[i], _candidates.back());
                ValueType ret = _candidates.back();
                _candidates.pop_back();
                return ret;
            }
            else
            {
                return _candidates[i];
            }
        }
        else
        {
            size_t i;
            do
            {
                if (_probs.back() > 0)
                {
                    tr1::variate_generator<rng_t&, tr1::uniform_real<> >
                        sample(rng, tr1::uniform_real<>(0.0, _probs.back()));
                    double r = sample();
                    i = upper_bound(_probs.begin(), _probs.end(),r) -
                        _probs.begin();
                }
                else
                {
                    // all probabilities are zero... sample randomly.
                    tr1::uniform_int<size_t>
                        sample(0, _candidates_set.size()-1);
                    size_t j = sample(rng), count = 0;
                    for (typeof(_candidates_set.begin()) iter =
                             _candidates_set.begin();
                         iter != _candidates_set.end(); ++iter)
                    {
                        if (count == j)
                        {
                            i = iter->second;
                            break;
                        }
                        count++;
                    }
                }
            } while (_erased[i]);

            if (remove)
            {
                _erased[i] = true;
                _erased_prob += (i > 0) ? _probs[i] - _probs[i-1] : _probs[i];
                clean();
            }

            return _candidates[i];
        }
    }

    void clean()
    {
        // if too many elements were erased, we need to make things less sparse
        if (_biased && !_candidates_set.empty() &&
            _erased_prob >= _probs.back()/3)
        {
            for (int i = _probs.size() - 1; i > 0; --i)
                _probs[i] -= _probs[i-1];

            for (int i = 0; i < _candidates.size(); ++i)
            {
                while (i < _erased.size() && _erased[i])
                {
                    swap(_candidates[i], _candidates.back());
                    _candidates.pop_back();

                    swap(_probs[i], _probs.back());
                    _probs.pop_back();

                    swap(_erased[i], _erased.back());
                    _erased.pop_back();
                }
            }

            for (size_t i = 1; i < _probs.size(); i++)
                _probs[i] += _probs[i-1];

            _candidates_set.clear();
            for (size_t i = 0; i < _candidates.size(); i++)
                _candidates_set.insert(make_pair(_candidates[i],i));
            _erased_prob = 0.0;
        }
    }

private:
    bool _biased;
    vector<ValueType> _candidates;
    tr1::unordered_multimap<ValueType, size_t, hash<ValueType> >
        _candidates_set;
    vector<double> _probs;
    vector<uint8_t> _erased;
    double _erased_prob;
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

class CorrelatedStrat
{
public:
    typedef pair<size_t,size_t> deg_t;
    typedef tr1::unordered_multimap<deg_t, dvertex_t, hash<deg_t> > targets_t;

    template <class CorrProb>
    void InsertSource(const dvertex_t& v, CorrProb& corr)
    {
        deg_t deg = make_pair(v.in_degree, v.out_degree);
        if (_sampler.find(deg) == _sampler.end())
        {
            _sampler[deg] = Sampler<deg_t>(true);
            for (typeof(_target_degs.begin()) iter = _target_degs.begin();
                 iter != _target_degs.end(); ++iter)
            {
                typeof(_corr_prob.begin()) prob =
                     _corr_prob.find(make_pair(deg, *iter));
                if (prob == _corr_prob.end())
                    prob = _corr_prob.insert(make_pair(make_pair(deg, *iter),
                                                       corr(deg, *iter)))
                        .first;
                _sampler[deg].Insert(*iter, prob->second);
            }
        }
    }

    template <class CorrProb>
    void InsertTarget(const dvertex_t& v, CorrProb& corr)
    {
        deg_t deg = make_pair(v.in_degree, v.out_degree);
        _targets.insert(make_pair(deg, v));

        // if it is a new target degree, include it on all previously included
        // samplers
        if (_target_degs.find(deg) == _target_degs.end())
        {
            for (typeof(_sampler.begin()) iter = _sampler.begin();
                 iter != _sampler.end(); ++iter)
            {
                typeof(_corr_prob.begin()) prob =
                     _corr_prob.find(make_pair(iter->first, deg));
                if (prob == _corr_prob.end())
                    prob =
                        _corr_prob.insert(make_pair(make_pair(iter->first, deg),
                                                    corr(iter->first, deg)))
                        .first;
                iter->second.Insert(deg, prob->second);
            }
            _target_degs.insert(deg);
        }
    }

    dvertex_t GetRandomTarget(const dvertex_t& source, rng_t& rng)
    {
        deg_t source_deg = make_pair(source.in_degree, source.out_degree),
            target_deg;

        //choose the target degree according to correlation
        target_deg = _sampler[source_deg](rng);

        // sample random target on the same class
        targets_t::iterator begin, end;
        tie(begin, end) = _targets.equal_range(target_deg);
        Sampler<dvertex_t> sampler;
        for(targets_t::iterator iter = begin; iter != end; ++iter)
            sampler.Insert(iter->second);
        dvertex_t v = sampler(rng);
        return v;
    }

    void Remove(const dvertex_t& v)
    {
        deg_t deg = make_pair(v.in_degree, v.out_degree);

        typeof(_targets.begin()) iter, end;
        tie(iter, end) = _targets.equal_range(deg);
        for (; iter != end; ++iter)
            if (iter->second == v)
            {
                Remove(iter);
                break;
            }
    }

    void Remove(const targets_t::iterator& i)
    {
        deg_t deg = i->first;

        _targets.erase(i);

        // remove degree from samplers if it no longer corresponds to any target
        if (_targets.find(deg) == _targets.end())
        {
            for (typeof(_sampler.begin()) iter = _sampler.begin();
                 iter != _sampler.end(); ++iter)
                iter->second.Remove(deg);
            _target_degs.erase(deg);
        }
    }

private:
    typedef tr1::unordered_map<pair<deg_t, deg_t>, double,
                               hash<pair<deg_t, deg_t> > > corr_prob_t;
    typedef tr1::unordered_map<deg_t, Sampler<deg_t>, hash<deg_t> > sampler_t;
    typedef tr1::unordered_set<deg_t, hash<deg_t> > deg_set_t;

    corr_prob_t _corr_prob;
    targets_t _targets;
    sampler_t _sampler;
    deg_set_t _target_degs;
};

//
// Uncorrelated graph strategy
//

class UncorrelatedStrat
{
public:
    typedef pair<size_t,size_t> deg_t;

    template <class CorrProb>
    void InsertSource(const dvertex_t& v, CorrProb& corr)
    {
    }

    template <class CorrProb>
    void InsertTarget(const dvertex_t& v, CorrProb& corr)
    {
        _targets.Insert(v);
    }

    dvertex_t GetRandomTarget(const dvertex_t& source, rng_t& rng)
    {
        return _targets(rng);
    }

    void Remove(const dvertex_t& v)
    {
        _targets.Remove(v);
    }

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

    template <class Graph, class DegSample, class CorrProb>
    void operator()(Graph& g, DegSample& deg_sample,
                    CorrProb& corr_prob, bool no_parallel,
                    bool no_self_loops, bool undirected,
                    size_t seed, bool verbose)
        const
    {
        size_t N = _N;
        typename property_map<Graph,vertex_index_t>::type vertex_index =
            get(vertex_index_t(), g);
        rng_t rng(static_cast<rng_t::result_type>(seed));

        // figure out the necessary strategies
        typedef typename mpl::if_<typename is_directed::apply<Graph>::type,
                                  DirectedStrat,
                                  UndirectedStrat>::type gen_strat_t;
        typedef typename gen_strat_t::deg_t deg_t;
        typedef typename mpl::if_<IsCorrelated, CorrelatedStrat,
                                  UncorrelatedStrat >::type corr_strat_t;

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
            corr_strat.InsertSource(vertices[i], corr_prob);
            if (gen_strat.IsTarget(vertices[i]))
                corr_strat.InsertTarget(vertices[i], corr_prob);
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
            if (!is_directed::apply<Graph>::type::value &&
                !gen_strat.IsStillTarget(g, source))
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
            target = corr_strat.GetRandomTarget(source, rng);

            // if unacceptable (because parallel and/or self_loop), throw it
            // back and disconnect a random previous (different) source
            if ((no_parallel && (old_targets.find(target.index)
                                 != old_targets.end())) ||
                (no_self_loops && (source.index == target.index)))
            {
                if (i < 2) // just try again
                {
                    i--;
                    continue;
                }

                // We need to check if all previous sources are not all equal to
                // the first one. In that case, we should just discard this one
                // and keep going.
                int j;
                for (j = i; j >= 0; --j)
                    if (!(sources[j] == source))
                        break;
                if (j < 0)
                {
                    i--;
                    continue; // just try again
                }

                size_t s_index;
                tr1::uniform_int<size_t> source_sample(0, i-1);
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
                tr1::uniform_int<size_t>
                    esample(0, out_degree(vertex(rsource.index, g),g)-1);
                size_t rei = esample(rng);
                typename graph_traits<Graph>::out_edge_iterator e, e_end;
                tie(e, e_end) = out_edges(vertex(rsource.index, g), g);
                for(size_t j = 0; j < rei; ++j)
                    e++;
                typename graph_traits<Graph>::edge_descriptor re = *e;

                // if edge's target was already full, put it back in target list
                dvertex_t rtarget =
                    vertices[vertex_index[boost::target(re, g)]];
                if (!gen_strat.IsStillTarget(g, rtarget))
                    corr_strat.InsertTarget(rtarget, corr_prob);

                // for undirected graphs, we need also to check the source
                if (!is_directed::apply<Graph>::type::value)
                {
                    dvertex_t rsource =
                        vertices[vertex_index[boost::source(re, g)]];
                    if (!gen_strat.IsStillTarget(g, rsource) &&
                        !(rsource == rtarget))
                        corr_strat.InsertTarget(rsource, corr_prob);
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
                corr_strat.Remove(target);

            if (!is_directed::apply<Graph>::type::value && !(source == target))
                if (!gen_strat.IsStillTarget(g, source))
                    corr_strat.Remove(source);

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
