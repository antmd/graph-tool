// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2013 Tiago de Paula Peixoto <tiago@skewed.de>
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

#ifndef SAMPLER_HH
#define SAMPLER_HH

#include "random.hh"
#include <iostream>

namespace graph_tool
{
using namespace std;
using namespace boost;

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
            Insert(*iter);
    }

    void Insert(const ValueType& v)
    {
        _candidates.push_back(v);
        _candidates_set.insert(make_pair(v, _candidates.size() - 1));
    }

    bool HasValue(const ValueType& v)
    {
        typeof(_candidates_set.begin()) iter, end;
        tie(iter, end) = _candidates_set.equal_range(v);
        return (iter != end);
    }

    void Remove(const ValueType& v)
    {
        typeof(_candidates_set.begin()) iter, back;
        iter = _candidates_set.find(v);

        if (iter == _candidates_set.end())
            return;

        back = _candidates_set.find(_candidates.back());

        size_t index = iter->second;
        swap(_candidates[index], _candidates.back());
        _candidates.pop_back();

        if (!_candidates.empty() && back != iter)
        {
            _candidates_set.erase(back);
            _candidates_set.insert(make_pair(_candidates[index], index));
        }
        _candidates_set.erase(iter);
    }

    bool Empty()
    {
        return _candidates.empty();
    }

    size_t Size()
    {
        return _candidates.size();
    }

    ValueType operator()(rng_t& rng, bool remove = false)
    {
        //assert(!_candidates.empty());
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

private:
    vector<ValueType> _candidates;
    tr1::unordered_multimap<ValueType, size_t, hash<ValueType> >
        _candidates_set;
};


template <class ValueType>
class WeightedSampler
{
public:

    void Insert(const ValueType& v, double p)
    {
        _candidates.push_back(make_pair(v, p));
        _candidates_set.insert(make_pair(v, _candidates.size() - 1));
        _erased.push_back(false);
        _rebuild = true;
    }

    bool HasValue(const ValueType& v)
    {
        typeof(_candidates_set.begin()) iter, end;
        tie(iter, end) = _candidates_set.equal_range(v);
        return (iter != end);
    }

    void Remove(const ValueType& v)
    {
        typeof(_candidates_set.begin()) iter, end, temp;
        tie(iter, end) = _candidates_set.equal_range(v);
     
        if (iter == end)
            return;

        while(_erased[iter->second])
        {
            temp = iter++;
            _candidates_set.erase(temp);
            if (iter == end)
                return;
        }

        size_t index = iter->second;
        _erased[index] = true;
        _erased_prob += _candidates[index].second;
        if (_erased_prob >= 0.3)
            _rebuild = true;
    }

    bool Empty()
    {
        return _candidates.empty();
    }

    size_t Size()
    {
        return _candidates.size();
    }

    void BuildTable()
    {
        // remove possibly erased elements
        size_t i = 0;
        while (i < _candidates.size())
        {
            if (_erased[i])
            {
                swap(_candidates[i], _candidates.back());
                swap(_erased[i], _erased.back());
                _candidates.pop_back();
                _erased.pop_back();
            }
            else
            {
                ++i;
            }
        }
        _erased_prob = 0;

        vector<pair<size_t, double> > remainder;
        _alias.resize(_candidates.size());

        double P_sum = 0;
        for (size_t i = 0; i < _candidates.size(); ++i)
            P_sum += _candidates[i].second;

        
        size_t N = _candidates.size();
        double P = 1.0 / N;
        for (size_t i = 0; i < _candidates.size(); ++i)
        {
            _candidates[i].second /= P_sum;
            double pi = _candidates[i].second;
            if (pi > P)
                remainder.push_back(make_pair(i, pi - P));
            _alias[i] = make_pair(i, .1);
        }


        for (size_t i = 0; i < _candidates.size(); ++i)
        {
            double pi = _candidates[i].second;
            if (pi < P)
            {
                for (size_t j = 0; j < remainder.size(); ++j)
                {
                    if (remainder[j].second >= P - pi)
                    {
                        _alias[i] = make_pair(remainder[j].first, pi * N);
                        remainder[j].second -= P - pi;
                        if (remainder[j].second <= 0)
                        {
                            swap(remainder[j], remainder.back());
                            remainder.pop_back();
                        }
                        break;
                    }
                }
            }
        }
        _rebuild = false;
    }


    ValueType operator()(rng_t& rng, bool remove = false)
    {
        if (_rebuild)
            BuildTable();

        tr1::variate_generator<rng_t&, tr1::uniform_real<> >
            sample(rng, tr1::uniform_real<>(0.0, 1.0));
        size_t i;
        do
        {
            double r = sample() * _candidates.size();
            i = floor(r);       // in [0, n-1]
            double x = r - i;   // in [0, 1)

            if (x > _alias[i].second)
                i = _alias[i].first;
        }
        while (_erased[i]);

        if (remove)
        {
            _erased[i] = true;
            _erased_prob += _candidates[i].second;
            if (_erased_prob >= 0.3)
                _rebuild = true;
        }
        return _candidates[i].first;
    }

private:
    vector<pair<ValueType, double> > _candidates;
    tr1::unordered_multimap<ValueType, size_t, hash<ValueType> >
        _candidates_set;
    vector<pair<size_t, double> > _alias;
    vector<uint8_t> _erased;
    bool _erased_prob;
    bool _rebuild;
};

} // namespace graph_tool

#endif // SAMPLER_HH
