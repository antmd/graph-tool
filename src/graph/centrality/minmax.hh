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

#ifndef MINMAX_HH
#define MINMAX_HH

#include <vector>
#include <boost/mpl/bool.hpp>
#include <iostream>
using namespace  std;

template <class Value, class Compare = std::less<Value> >
class double_priority_queue
{
public:
    double_priority_queue() {}
    double_priority_queue(const Compare& cmp): _cmp(cmp) {}

    size_t size() { return _queue.size(); }
    bool empty() { return _queue.empty(); }

    void push(const Value& v)
    {
        _queue.push_back(v);
        bubble_up(_queue.size()-1);
    }

    void pop_top()
    {
        // for (size_t i = 0; i < _queue.size(); ++i)
        //     if (_cmp(top(), _queue[i]))
        //     {
        //         cout << size_t(log2(i+1)) << endl;
        //         cout << "ops top" << endl;
        //         abort();
        //     }
        size_t i = &top() - &_queue[0];
        swap(_queue[i], _queue.back());
        _queue.pop_back();
        if (!_queue.empty())
            trickle_down(i);
    }

    void pop_bottom()
    {
        swap(_queue.front(), _queue.back());
        _queue.pop_back();
        if (!_queue.empty())
            trickle_down(0);
    }

    const Value& top()
    {
        if (_queue.size() < 2)
            return _queue[0];
        if (_queue.size() < 3)
            return _queue[1];
        if (_cmp(_queue[1], _queue[2]))
            return _queue[2];
        else
            return _queue[1];
    }

    const Value& bottom()
    {
        // for (size_t i = 0; i < _queue.size(); ++i)
        //     if (_cmp(_queue[i], _queue.front()))
        //         cout << "ops bottom" << endl;
        return _queue.front();
    }

    // general comparison: if IsMin is mpl::true_, then a < b, otherwise a > b
    // is computed
    template <class IsMin>
    bool cmp(size_t a, size_t b, IsMin)
    {
        if (IsMin::value)
            return _cmp(_queue[a], _queue[b]);
        else
            return _cmp(_queue[b], _queue[a]);
    }

    void trickle_down(size_t i)
    {
        if ((size_t(log2(i+1)) % 2) == 0)
            trickle_down(i, boost::mpl::true_());
        else
            trickle_down(i, boost::mpl::false_());
    }


    template <class IsMin>
    void trickle_down(size_t pos, IsMin is_min)
    {
        size_t i = pos;
        while (true)
        {
            size_t m = 2*i+1;
            if (m >= _queue.size())
                break;

            for (size_t j = 2*i+1; j < 2*i+3; ++j)
            {
                if (j >= _queue.size())
                    break;
                if (cmp(j, m, is_min))
                    m = j;
                for (size_t l = 2*j+1; l < 2*j+3; ++l)
                {
                    if (l >= _queue.size())
                        break;
                    if (cmp(l, m, is_min))
                        m = l;
                }
            }

            if (m > 2*i+2) // is grandchild
            {
                if (cmp(m, i, is_min))
                {
                    swap(_queue[m], _queue[i]);
                    if (!cmp(m, (m-1)/2, is_min))
                        swap(_queue[m], _queue[(m-1)/2]);
                    i = m;
                }
                else
                {
                    break;
                }
            }
            else
            {
                if (cmp(m, i, is_min))
                    swap(_queue[m], _queue[i]);
                break;
            }
        }
    }

    void bubble_up(size_t i)
    {
        if ((size_t(log2(i+1)) % 2) == 0)
        {
            if ((i > 0) && cmp(i, (i-1)/2, boost::mpl::false_()))
            {
                swap(_queue[i], _queue[(i-1)/2]);
                bubble_up((i-1)/2,  boost::mpl::false_());
            }
            else
            {
                bubble_up(i,  boost::mpl::true_());
            }
        }
        else
        {
            if ((i > 0) && cmp(i, (i-1)/2, boost::mpl::true_()))
            {
                swap(_queue[i], _queue[(i-1)/2]);
                bubble_up((i-1)/2, boost::mpl::true_());
            }
            else
            {
                bubble_up(i, boost::mpl::false_());
            }
        }
    }

    template <class IsMin>
    void bubble_up(size_t pos, IsMin is_min)
    {
        size_t i = pos;
        while (true)
        {
            if (i > 2) // has grandparent
            {
                size_t m = ((i-1)/2 - 1)/2;
                if (cmp(i, m, is_min))
                {
                    swap(_queue[m], _queue[i]);
                    i = m;
                }
                else
                {
                    break;
                }
            }
            else
            {
                break;
            }
        }
    }

private:
    std::vector<Value> _queue;
    Compare _cmp;
};


#endif // MINMAX_HH
