// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2015 Tiago de Paula Peixoto <tiago@skewed.de>
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

#include <unordered_map>
#include <tuple>
#include <boost/functional/hash.hpp>
#include <map>
#include <set>
#include <iostream>

#include "graph_util.hh"
#include "random.hh"

namespace graph_tool
{
using namespace std;
using namespace boost;

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


// used for verbose display
void print_progress(size_t current, size_t total, stringstream& str)
{
    size_t atom = (total > 200) ? total/100 : 1;
    if ( ( (current+1) % atom == 0) || (current + 1) == total)
    {
        for (size_t j = 0; j < str.str().length(); ++j)
            cout << "\b";
        str.str("");
        str << current+1 << " of " << total << " ("
            << (current+1)*100 / total << "%)";
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
            _max_deg = _no_self_loops ? _N - 1 : N;
        else
            _max_deg = numeric_limits<size_t>::max();
    }

    // check whether a degree sequence is graphical
    template <class DegSequence>
    bool is_graphical(DegSequence& d)
    {
        size_t one = _no_self_loops ? 1 : 0;

        size_t sum_out_deg = 0;
        size_t j = 0;
        for(typeof(d.begin()) di = d.begin(); di != d.end(); ++di)
        {
            size_t out_deg = di->first.second;
            size_t count = di->second;

            j += count;
            sum_out_deg += out_deg * count;

            size_t sum_in_deg = 0;
            typeof(d.begin()) dj_end = di; ++dj_end;
            for(typeof(d.begin()) dj = d.begin(); dj != dj_end; ++dj)
                sum_in_deg += min(dj->first.first, j - one)*dj->second;

            size_t sum_rest = 0;
            typeof(d.begin()) dj = di;
            for(++dj; dj != d.end(); ++dj)
                sum_rest += min(dj->first.first, j)*dj->second;

            if (sum_out_deg > sum_in_deg + sum_rest)
                return false;
        }
        return true;
    }

    // check whether a degree sequence is graphical, when parallel loops are
    // allowed but self-loops are not
    template <class DegSequence>
    bool is_graphical_parallel(DegSequence& d)
    {
        size_t sum_in_deg = 0;
        for(typeof(d.begin()) di = d.begin(); di != d.end(); ++di)
            sum_in_deg += di->first.first * di->second;
        for(typeof(d.begin()) di = d.begin(); di != d.end(); ++di)
        {
            if (di->first.second > sum_in_deg - di->first.first)
                return false;
        }
        return true;
    }

    // sample the degrees of all the vertices
    template <class DegSample>
    size_t SampleDegrees(vector<dvertex_t>& vertices, DegSample& deg_sample,
                         rng_t& rng, bool verbose)
    {
        stringstream str;
        size_t sum_j=0, sum_k=0;
        for(size_t i = 0; i < _N; ++i)
        {
            if (verbose)
                print_progress(i, _N, str);
            dvertex_t& v = vertices[i];
            do
            {
                tie(v.in_degree, v.out_degree) = deg_sample(i);
            }
            while (_no_parallel &&
                   (v.in_degree > _max_deg || v.out_degree > _max_deg));
            sum_j += v.in_degree;
            sum_k += v.out_degree;
            if (_no_parallel || _no_self_loops)
                _deg_seq[make_pair(v.in_degree, v.out_degree)]++;
        }

        if (verbose)
        {
            cout << "\nfixing average degrees. Total degree difference: "
                 << flush;
                str.str("");
        }

        // Sequence must be graphical. Re-sample random pairs until this holds
        uniform_int_distribution<size_t> vertex_sample(0, _N-1);
        size_t count = 0;
        while(sum_j != sum_k || (_no_parallel && !is_graphical(_deg_seq)) ||
              (_no_self_loops && !_no_parallel &&
               !is_graphical_parallel(_deg_seq)))
        {
            size_t i = vertex_sample(rng);
            dvertex_t& v = vertices[i];
            if (_no_parallel || _no_self_loops)
            {
                typeof(_deg_seq.begin()) iter =
                    _deg_seq.find(make_pair(v.in_degree, v.out_degree));
                iter->second--;
                if (iter->second == 0)
                    _deg_seq.erase(iter);
            }

            sum_j -= v.in_degree;
            sum_k -= v.out_degree;
            do
            {
                tie(v.in_degree, v.out_degree) = deg_sample(i);
            }
            while (_no_parallel &&
                   (v.in_degree > _max_deg || v.out_degree > _max_deg));
            sum_j += v.in_degree;
            sum_k += v.out_degree;
            if (_no_parallel || _no_self_loops)
                _deg_seq[make_pair(v.in_degree, v.out_degree)]++;
            if (verbose && (count % 100 == 0 || sum_j == sum_k))
                print_update(min(sum_j-sum_k, sum_k-sum_j), str);
            count++;
        }
        return sum_k;
    }

    struct deg_cmp
    {
        bool operator()(const deg_t& d1, const deg_t& d2) const
        {
            if (d1.second == d2.second)
                return d1.first > d2.first;
            return d1.second > d2.second;
        }
    };

private:
    size_t _N;
    bool _no_parallel;
    bool _no_self_loops;
    size_t _max_deg;
    map<deg_t, size_t, deg_cmp> _deg_seq;
};

//
// Undirected graph generation strategy
//

class UndirectedStrat
{
public:
    typedef size_t deg_t; // degree type

    UndirectedStrat(size_t N, bool no_parallel, bool no_self_loops)
        : _N(N), _no_parallel(no_parallel), _no_self_loops(no_self_loops)
    {
        if (_no_parallel)
            _max_deg = (_no_self_loops) ? _N - 1 : N;
        else
            _max_deg = numeric_limits<size_t>::max();
    }

    // check whether a degree sequence is graphical
    template <class DegSequence>
    bool is_graphical(DegSequence& d)
    {
        size_t one = (_no_self_loops) ? 1 : 0;
        size_t sum_deg = 0;
        size_t j = 0;
        for(typeof(d.begin()) di = d.begin(); di != d.end(); ++di)
        {
            j += di->second;
            sum_deg += di->first * di->second;
            size_t sum_rest = 0;
            typeof(d.begin()) dj = di;
            for(++dj; dj != d.end(); ++dj)
                sum_rest += min(dj->first, j)*dj->second;
            if (sum_deg > j*(j-one) + sum_rest)
                return false;
        }
        return true;
    }

    // check whether a degree sequence is graphical, when parallel loops are
    // allowed but self-loops are not
    template <class DegSequence>
    bool is_graphical_parallel(DegSequence& d)
    {
        size_t sum_deg = 0;
        for(typeof(d.begin()) di = d.begin(); di != d.end(); ++di)
            sum_deg += di->first * di->second;
        for(typeof(d.begin()) di = d.begin(); di != d.end(); ++di)
        {
            if (di->first > sum_deg - di->first)
                return false;
        }
        return true;
    }

    // samples the degress of all vertices
    template <class DegSample>
    size_t SampleDegrees(vector<dvertex_t>& vertices, DegSample& deg_sample,
                         rng_t& rng, bool verbose)
    {
        stringstream str;
        size_t sum_k=0;
        sum_k = 0;
        for(size_t i = 0; i < _N; ++i)
        {
            if (verbose)
                print_progress(i, _N, str);
            dvertex_t& v = vertices[i];
            do
            {
                v.out_degree = deg_sample(i, true);
            }
            while (_no_parallel && v.out_degree > _max_deg);
            sum_k += v.out_degree;
            _deg_seq[v.out_degree]++;
        }

        if (verbose)
        {
            cout << "\nFixing degree sequence: "
                 << flush;
            str.str("");
        }

        // sum_k must be an even number (2*num_edges), and degree sequence must
        // be graphical, if multiple edges are not allowed. Re-sample degrees
        // until this holds.
        uniform_int_distribution<size_t> vertex_sample(0, _N-1);
        size_t count = 0;
        while (sum_k % 2 != 0 || (_no_parallel && !is_graphical(_deg_seq)) ||
               (_no_self_loops && !_no_parallel &&
                !is_graphical_parallel(_deg_seq)))
        {
            size_t i = vertex_sample(rng);
            dvertex_t& v = vertices[i];
            if (_no_parallel || _no_self_loops)
            {
                typeof(_deg_seq.begin()) iter = _deg_seq.find(v.out_degree);
                iter->second--;
                if(iter->second == 0)
                    _deg_seq.erase(iter);
            }
            sum_k -= v.out_degree;
            do
            {
                v.out_degree = deg_sample(i, true);
            }
            while (_no_parallel && (v.out_degree > _max_deg));
            sum_k +=  v.out_degree;
            if (_no_parallel || _no_self_loops)
                _deg_seq[v.out_degree]++;
            if (verbose && (count % 100 || sum_k % 2 == 0))
                print_update(sum_k, str);
            count++;
        }
        return sum_k/2;
    }

private:
    size_t _N;
    bool _no_parallel;
    bool _no_self_loops;
    size_t _max_deg;
    map<size_t, size_t, greater<size_t> > _deg_seq;
};

//
// Main Algorithm
// ==============
//
// generates a directed or undirected graph with given degree distribution

template <class Cmp>
struct cmp_in
{
    bool operator()(const pair<size_t,size_t>& v1,
                    const pair<size_t,size_t>& v2) const
    {
        if (v1.first == v2.first)
            return cmp(v1.second, v2.second);
        return cmp(v1.first, v2.first);
    }
    Cmp cmp;
};

template <class Cmp>
struct cmp_out
{
    bool operator()(const pair<size_t,size_t>& v1,
                    const pair<size_t,size_t>& v2) const
    {
        if (v1.second == v2.second)
            return cmp(v1.first, v2.first);
        return cmp(v1.second, v2.second);
    }
    Cmp cmp;
};

template <class Graph>
pair<size_t, size_t> get_deg(dvertex_t& v, Graph& g)
{
    return make_pair(v.in_degree - in_degreeS()(vertex(v.index, g), g),
                     v.out_degree - out_degree(vertex(v.index, g), g));
}

template <class Graph>
bool is_source(const pair<size_t, size_t>& deg)
{
    return deg.second > 0;
}

template <class Graph>
bool is_target(const pair<size_t, size_t>& deg)
{
    if (is_directed::apply<Graph>::type::value)
        return deg.first > 0;
    else
        return is_source<Graph>(deg);
}


template <class Vset, class Targets, class Sources, class Graph>
bool update_deg(size_t t_i, const pair<size_t, size_t>& deg, Vset& vset,
                Targets& targets, Sources& sources, Graph&)
{
    if (is_source<Graph>(deg))
        sources.insert(deg);
    if (is_target<Graph>(deg))
        targets.insert(deg);
    vset[deg].push_back(t_i);
    return true;
}

struct gen_graph
{
    template <class Graph, class DegSample>
    void operator()(Graph& g, size_t N, DegSample& deg_sample, bool no_parallel,
                    bool no_self_loops, rng_t& rng, bool verbose, bool verify)
        const
    {
        typename property_map<Graph,vertex_index_t>::type vertex_index =
            get(vertex_index_t(), g);

        // figure out the necessary strategy
        typedef typename mpl::if_<typename is_directed::apply<Graph>::type,
                                  DirectedStrat,
                                  UndirectedStrat>::type gen_strat_t;

        gen_strat_t gen_strat(N, no_parallel, no_self_loops);
        stringstream str; // used for verbose status

        if (verbose)
            cout << "adding vertices: " << flush;

        vector<dvertex_t> vertices(N);
        for(size_t i = 0; i < N; ++i)
            vertices[i].index = vertex_index[add_vertex(g)];

        // sample the N (j,k) pairs
        size_t E = gen_strat.SampleDegrees(vertices, deg_sample, rng, verbose);

        // source and target degree lists
        typedef pair<size_t, size_t> deg_t;
        set<deg_t, cmp_out<greater<size_t> > > sources;
        set<deg_t, cmp_in<greater<size_t> > > targets;

        // vertices with a given degree
        unordered_map<deg_t, vector<size_t>,
                      boost::hash<deg_t> > vset;

        size_t num_e = 0;
        for (size_t i = 0; i < vertices.size();  ++i)
        {
            deg_t deg = get_deg(vertices[i], g);

            if (is_source<Graph>(deg))
                sources.insert(deg);
            if (is_target<Graph>(deg))
                targets.insert(deg);
            if (is_target<Graph>(deg) || is_source<Graph>(deg))
                vset[deg].push_back(i);
        }

        if (verbose)
        {
            cout << endl << "adding edges: " << flush;
            str.str("");
        }

        vector<size_t> skip;

        // connect edges: from sources with the largest in-degree to the ones
        // with largest out-degree
        while (!sources.empty())
        {
            // find source. The out-degree must be non-zero, and there must be a
            // vertex with the chosen degree.
            deg_t s_deg = *sources.begin();
            typeof(vset.begin()) sv_iter = vset.find(s_deg);
            if (s_deg.second == 0 || sv_iter == vset.end() ||
                sv_iter->second.empty())
            {
                sources.erase(sources.begin());
                continue;
            }

            vector<size_t>& s_list = sv_iter->second;
            size_t s_i = s_list.front();
            typename graph_traits<Graph>::vertex_descriptor s =
                vertex(vertices[s_i].index, g);

            deg_t ns_deg = get_deg(vertices[s_i], g);
            if (ns_deg != s_deg)
            {
                swap(s_list.back(), s_list.front());
                s_list.pop_back();
                update_deg(s_i, ns_deg, vset, targets, sources, g);
                continue;
            }

            // find the targets.
            // we will keep an iterator to the current target degree
            typeof(targets.begin()) t_iter = targets.begin();
            typeof(vset.begin()) v_iter = vset.find(*t_iter);
            while (v_iter == vset.end() || v_iter->second.empty())
            {
                targets.erase(t_iter);
                t_iter = targets.begin();
                v_iter = vset.find(*t_iter);
            }

            skip.clear();
            skip.push_back(s_i);

            if (no_self_loops)
            {
                swap(s_list.back(), s_list.front());
                s_list.pop_back();
            }

            while (s_deg.second > 0)
            {
                //assert(!targets.empty());
                //assert(t_iter != targets.end());

                while (v_iter == vset.end() || v_iter->second.empty())
                {
                    ++t_iter;
                    v_iter = vset.find(*t_iter);
                }

                deg_t t_deg = *t_iter;

                vector<size_t>& v_list = v_iter->second;

                size_t t_i = v_list.front();

                deg_t nt_deg = get_deg(vertices[t_i], g);
                if (nt_deg != t_deg)
                {
                    swap(v_list.back(), v_list.front());
                    v_list.pop_back();
                    update_deg(t_i, nt_deg, vset, targets, sources, g);
                    //t_iter = targets.begin();
                    //v_iter = vset.find(*t_iter);
                    continue;
                }

                // remove target from vertex list, and get new t_i
                skip.push_back(t_i);

                swap(v_list.back(), v_list.front());
                v_list.pop_back();

                typename graph_traits<Graph>::vertex_descriptor t =
                    vertex(vertices[t_i].index, g);

                if ((s == t) && (!is_directed::apply<Graph>::type::value &&
                                 s_deg.second < 2))
                    continue;

                add_edge(s, t, g);
                s_deg = get_deg(vertices[s_i], g);

                // if parallel edges are allowed, we should update the target
                // list right away
                if (!no_parallel)
                {
                    for (size_t i = 0; i < skip.size(); ++i)
                    {
                        if (no_self_loops && skip[i] == s_i)
                            continue;
                        update_deg(skip[i],
                                   get_deg(vertices[skip[i]], g), vset,
                                   targets, sources, g);
                    }
                    skip.clear();
                    if (no_self_loops)
                        skip.push_back(s_i);
                }
                if (verbose)
                    print_progress(num_e++, E, str);
            }

            if (!s_list.empty() && s_list.front() == s_i)
            {
                swap(s_list.back(), s_list.front());
                s_list.pop_back();
            }

            // update modified degrees
            for (size_t i = 0; i < skip.size(); ++i)
                update_deg(skip[i],
                           get_deg(vertices[skip[i]], g),
                           vset, targets, sources, g);
       }
        if (verbose)
            cout << endl;

        if (verify)
        {
            for (size_t i = 0; i < vertices.size(); ++i)
            {
                deg_t dseq = make_pair(vertices[i].in_degree,
                                       vertices[i].out_degree);
                deg_t deg = make_pair(in_degreeS()(vertex(i, g), g),
                                      out_degree(vertex(i, g), g));
                if (deg != dseq)
                    throw GraphException("Graph does not match the desired "
                                         "sequence! Vertex " +
                                         lexical_cast<string>(i) +
                                         ", wanted: " +
                                         lexical_cast<string>(dseq.first) +
                                         " " +
                                         lexical_cast<string>(dseq.second) +
                                         ", got: " +
                                         lexical_cast<string>(deg.first) +
                                         " " +
                                         lexical_cast<string>(deg.second) +
                                         " This is a bug.");
            }
        }
    }
};

} // graph_tool namespace

#endif // GRAPH_GENERATION_HH
