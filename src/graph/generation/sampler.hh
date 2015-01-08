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

#ifndef SAMPLER_HH
#define SAMPLER_HH

#include "random.hh"
#include <functional>
#include <boost/mpl/if.hpp>

namespace graph_tool
{
using namespace std;
using namespace boost;

// Discrete sampling via vose's alias method.

// See http://www.keithschwarz.com/darts-dice-coins/ for a very clear
// explanation.

template <class Value, class KeepReference = mpl::true_>
class Sampler
{
public:
    Sampler(const vector<Value>& items,
            const vector<double>& probs)
        : _items(items), _probs(probs), _alias(items.size())
    {
        double S = 0;
        for (size_t i = 0; i < _probs.size(); ++i)
            S += _probs[i];

        for (size_t i = 0; i < _probs.size(); ++i)
        {
            _probs[i] *= _probs.size() / S;
            if (_probs[i] < 1)
                _small.push_back(i);
            else
                _large.push_back(i);
        }

        while (!(_small.empty() || _large.empty()))
        {
            size_t l = _small.back();
            size_t g = _large.back();
            _small.pop_back();
            _large.pop_back();

            _alias[l] = g;
            _probs[g] = (_probs[l] + _probs[g]) - 1;
            if (_probs[g] < 1)
                _small.push_back(g);
            else
                _large.push_back(g);
        }

        // fix numerical instability
        for (size_t i = 0; i < _large.size(); ++i)
            _probs[_large[i]] = 1;
        for (size_t i = 0; i < _small.size(); ++i)
            _probs[_small[i]] = 1;
        _large.clear();
        _small.clear();
    }

    Sampler() {}

    template <class RNG>
    const Value& sample(RNG& rng)
    {
        uniform_int_distribution<size_t> sample(0, _probs.size() - 1);
        size_t i = sample(rng);

        bernoulli_distribution coin(_probs[i]);
        if (coin(rng))
            return _items[i];
        else
            return _items[_alias[i]];
    }

    size_t size() const { return _items.size(); }

private:

    typedef typename mpl::if_<KeepReference,
                              const vector<Value>&,
                              vector<Value> >::type items_t;
    items_t _items;
    vector<double> _probs;
    vector<size_t> _alias;
    vector<size_t> _small;
    vector<size_t> _large;
};

// uniform sampling from containers

template <class Container, class RNG>
const typename Container::value_type& uniform_sample(const Container& v, RNG& rng)
{
    std::uniform_int_distribution<size_t> i_rand(0, v.size() - 1);
    return v[i_rand(rng)];
}


} // namespace graph_tool

#endif // SAMPLER_HH
