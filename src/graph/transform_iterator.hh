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

#ifndef TRANSFORM_ITERATOR_HH
#define TRANSFORM_ITERATOR_HH

#include <iterator>
#include <boost/iterator/transform_iterator.hpp>

template <class Predicate, class Iterator>
class transform_random_access_iterator:
    public boost::transform_iterator<Predicate, Iterator>
{
public:
    typedef Iterator iter_t;
    typedef boost::transform_iterator<Predicate, Iterator> base_t;
    transform_random_access_iterator() {}
    transform_random_access_iterator(const base_t& iter) :base_t(iter) {}
    transform_random_access_iterator(const Iterator& iter, const Predicate& pred = Predicate())
        : base_t(iter, pred) {}
};

namespace std
{
template <class Predicate, class Iterator>
struct iterator_traits<transform_random_access_iterator<Predicate, Iterator> >
{
    typedef transform_random_access_iterator<Predicate, Iterator> titer_t;
    typedef typename titer_t::base_t base_t;
    typedef typename iterator_traits<base_t>::difference_type difference_type;
    typedef typename iterator_traits<base_t>::value_type value_type;
    typedef typename iterator_traits<base_t>::value_type reference;
    typedef typename iterator_traits<base_t>::pointer pointer;
    typedef typename iterator_traits<typename titer_t::iter_t>::iterator_category iterator_category;
};
}

#endif // TRANSFORM_ITERATOR_HH
