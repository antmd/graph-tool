// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007-2010 Tiago de Paula Peixoto <tiago@skewed.de>
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

#if (GCC_VERSION >= 40400)
#   include <tr1/random>
#else
#   include <boost/tr1/random.hpp>
#endif
#include <iostream>

namespace graph_tool
{
using namespace std;
using namespace boost;

typedef tr1::mt19937 rng_t;

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
        //assert(!_candidates.empty());
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
        //assert(iter != end);

        if (_biased)
        {
            while(_erased[iter->second])
            {
                temp = iter++;
                _candidates_set.erase(temp);
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
            size_t i = 0;
            do
            {
                if (_probs.back() > 0)
                {
                    tr1::variate_generator<rng_t&, tr1::uniform_real<> >
                        sample(rng, tr1::uniform_real<>(0.0, _probs.back()));
                    double r = sample();
                    i = upper_bound(_probs.begin(), _probs.end(), r) -
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
            for (int i = int(_probs.size()) - 1; i > 0; --i)
                _probs[i] -= _probs[i-1];

            for (size_t i = 0; i < _candidates.size(); ++i)
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

} // namespace graph_tool

#endif // SAMPLER_HH
