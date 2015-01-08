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

#ifndef DYNAMIC_SAMPLER_HH
#define DYNAMIC_SAMPLER_HH

#include "random.hh"
#include <functional>
#include <boost/mpl/if.hpp>

namespace graph_tool
{
using namespace std;
using namespace boost;

template <class Value>
class DynamicSampler
{
public:
    DynamicSampler() : _back(0) {}

    DynamicSampler(const vector<Value>& items,
                   const vector<double>& probs)
        : _back(0)
    {
        for (size_t i = 0; i < items.size(); ++i)
            insert(items[i], probs[i]);
    }

    typedef Value value_type;

    size_t get_left(size_t i)   { return 2 * i + 1;               }
    size_t get_right(size_t i)  { return 2 * i + 2;               }
    size_t get_parent(size_t i) { return i > 0 ? (i - 1) / 2 : 0; }

    template <class RNG>
    const Value& sample(RNG& rng)
    {
        uniform_real_distribution<> sample(0, 1);
        double u = _tree[0] * sample(rng), c = 0;

        size_t pos = 0;
        while (_idx[pos] < 0)
        {
            size_t l = get_left(pos);
            double a = _tree[l];
            if (u < a + c)
            {
                pos = l;
            }
            else
            {
                pos = get_right(pos);
                c += a;
            }
        }
        size_t i = _idx[pos];
        return _items[i];
    }

    size_t insert(const Value& v, double w)
    {
        size_t pos;
        if (_free.empty())
        {
            if (_back > 0)
            {
                // move parent to left leaf
                pos = get_parent(_back);
                size_t l = get_left(pos);
                _idx[l] = _idx[pos];
                _ipos[_idx[l]] = l;
                _tree[l] = _tree[pos];
                _idx[pos] = -1;

                // position new item to the right
                _back = get_right(pos);
            }

            pos = _back;
            check_size(pos);

            _idx[pos] = _items.size();
            _items.push_back(v);
            _ipos.push_back(pos);
            _tree[pos] = w;
            _back++;
            check_size(_back);
        }
        else
        {
            pos = _free.back();
            _items[_idx[pos]] = v;
            _tree[pos] = w;
            _free.pop_back();
        }

        insert_leaf_prob(pos);
        return _idx[pos];
    }

    void remove(size_t i)
    {
        size_t pos = _ipos[i];
        remove_leaf_prob(pos);
        _free.push_back(pos);
    }

    void reset()
    {
        _items.clear();
        _tree.clear();
        _idx.clear();
        _back = 0;
        _free.clear();
    }

    void rebuild()
    {
        vector<Value> items;
        vector<double> probs;

        for (size_t i = 0; i < _tree.size(); ++i)
        {
            if (_idx[i] < 0)
                continue;
            items.push_back(_items[_idx[i]]);
            probs.push_back(_tree[i]);
        }

        reset();

        for (size_t i = 0; i < items.size(); ++i)
            insert(items[i], probs[i]);
    }

private:

    void check_size(size_t i)
    {
        if (i >= _tree.size())
        {
            _idx.resize(i + 1, -1);
            _tree.resize(i + 1, 0);
        }
    }

    void remove_leaf_prob(size_t i)
    {
        size_t parent = i;
        double w = _tree[i];
        while (parent > 0)
        {
            parent = get_parent(parent);
            _tree[parent] -= w;
        }
        _tree[i] = 0;
    }

    void insert_leaf_prob(size_t i)
    {
        size_t parent = i;
        double w = _tree[i];

        while (parent > 0)
        {
            parent = get_parent(parent);
            _tree[parent] += w;
        }
    }


    vector<Value> _items;
    vector<size_t> _ipos;   // position of the item in the tree

    vector<double> _tree;  // tree nodes with weight sums
    vector<int>    _idx;   // index in _items
    int _back;             // last item in tree

    vector<size_t> _free; // empty leafs
};



} // namespace graph_tool

#endif // DYNAMIC_SAMPLER_HH
