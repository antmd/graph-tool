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

// Copyright (C) Vladimir Prus 2003.
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/graph/vector_property_map.html for
// documentation.
//

//
// This is a modification of boost's vector property map which optionally
// disables bound checking for better performance.
//

#ifndef FAST_VECTOR_PROPERTY_MAP_HH
#define FAST_VECTOR_PROPERTY_MAP_HH

#include <boost/property_map.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>

namespace boost {

template<typename T, typename IndexMap>
class unchecked_fast_vector_property_map;

template<typename T, typename IndexMap = identity_property_map>
class fast_vector_property_map
    : public boost::put_get_helper<
              typename std::iterator_traits<
                  typename std::vector<T>::iterator >::reference,
              fast_vector_property_map<T, IndexMap> >
{
public:
    typedef typename property_traits<IndexMap>::key_type  key_type;
    typedef T value_type;
    typedef typename std::iterator_traits<
        typename std::vector<T>::iterator >::reference reference;
    typedef boost::lvalue_property_map_tag category;

    template<typename Type, typename Index>
    friend class unchecked_fast_vector_property_map;

    typedef unchecked_fast_vector_property_map<T, IndexMap> unchecked_t;
    typedef IndexMap index_map_t;
    typedef fast_vector_property_map<T,IndexMap> self_t;

    fast_vector_property_map(const IndexMap& index = IndexMap())
        : store(new std::vector<T>()), index(index)
    {}

    fast_vector_property_map(unsigned initial_size,
                             const IndexMap& index = IndexMap())
        : store(new std::vector<T>(initial_size)), index(index)
    {}

    typename std::vector<T>::iterator storage_begin()
    {
        return store->begin();
    }

    typename std::vector<T>::iterator storage_end()
    {
        return store->end();
    }

    typename std::vector<T>::const_iterator storage_begin() const
    {
        return store->begin();
    }

    typename std::vector<T>::const_iterator storage_end() const
    {
        return store->end();
    }

    void reserve(size_t size) const
    {
        if (store->size() < size)
            store->resize(size);
    }

    std::vector<T>& get_storage() const { return (*store); }

    unchecked_t get_unchecked(size_t size = 0) const
    {
        reserve(size);
        return unchecked_t(*this, size);
    }

public:
    // Copy ctor absent, default semantics is OK.
    // Assignment operator absent, default semantics is OK.
    // CONSIDER: not sure that assignment to 'index' is correct.

    reference operator[](const key_type& v) const {
        typename property_traits<IndexMap>::value_type i = get(index, v);
        if (static_cast<unsigned>(i) >= store->size()) {
            store->resize(i + 1, T());
        }
        return (*store)[i];
    }
protected:
    // Conceptually, we have a vector of infinite size. For practical
    // purposes, we start with an empty vector and grow it as needed.
    // Note that we cannot store pointer to vector here -- we cannot
    // store pointer to data, because if copy of property map resizes
    // the vector, the pointer to data will be invalidated.
    // I wonder if class 'pmap_ref' is simply needed.
    shared_ptr< std::vector<T> > store;
    IndexMap index;
};

template<typename T, typename IndexMap = identity_property_map>
class unchecked_fast_vector_property_map
    : public boost::put_get_helper<
                typename std::iterator_traits<
                    typename std::vector<T>::iterator >::reference,
                unchecked_fast_vector_property_map<T, IndexMap> >
{
public:
    typedef typename property_traits<IndexMap>::key_type  key_type;
    typedef T value_type;
    typedef typename std::iterator_traits<
        typename std::vector<T>::iterator >::reference reference;
    typedef boost::lvalue_property_map_tag category;

    typedef fast_vector_property_map<T, IndexMap> vmap_t;

    unchecked_fast_vector_property_map(const vmap_t& vmap = vmap_t(),
                                       size_t size = 0)
        : _vmap(vmap)
    {
        if (size > 0 && _vmap.store->size() < size)
            _vmap.store->resize(size);
    }

    void reserve(size_t size) const { _vmap.reserve(size); }

    reference operator[](const key_type& v) const
    {
        typename property_traits<IndexMap>::value_type i =
            get(_vmap.index, v);
        return (*_vmap.store)[i];
    }

private:

    vmap_t _vmap;
};

template<typename T, typename IndexMap>
fast_vector_property_map<T, IndexMap>
make_fast_vector_property_map(IndexMap index)
{
    return fast_vector_property_map<T, IndexMap>(index);
}

}

#endif
