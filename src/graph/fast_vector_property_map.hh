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

#include <boost/version.hpp>
#if (BOOST_VERSION >= 104000)
#   include <boost/property_map/property_map.hpp>
#else
#   include <boost/property_map.hpp>
#endif
#include <memory>
#include <vector>

namespace boost {

template<typename T, typename IndexMap>
class unchecked_vector_property_map;

template<typename T, typename IndexMap = identity_property_map>
class checked_vector_property_map
    : public boost::put_get_helper<
              typename std::iterator_traits<
                  typename std::vector<T>::iterator >::reference,
              checked_vector_property_map<T, IndexMap> >
{
public:
    typedef typename property_traits<IndexMap>::key_type  key_type;
    typedef T value_type;
    typedef typename std::iterator_traits<
        typename std::vector<T>::iterator >::reference reference;
    typedef boost::lvalue_property_map_tag category;

    template<typename Type, typename Index>
    friend class unchecked_vector_property_map;

    typedef unchecked_vector_property_map<T, IndexMap> unchecked_t;
    typedef IndexMap index_map_t;
    typedef checked_vector_property_map<T,IndexMap> self_t;

    checked_vector_property_map(const IndexMap& idx = IndexMap())
        : store(std::make_shared<std::vector<T>>()), index(idx) {}

    checked_vector_property_map(unsigned initial_size,
                                const IndexMap& idx = IndexMap())
        : store(std::make_shared<std::vector<T>>(initial_size)), index(idx) {}

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

    // deep copy
    checked_vector_property_map copy() const
    {
        checked_vector_property_map pmap(index);
        *(pmap.store) = *store;
        return pmap;
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
    std::shared_ptr< std::vector<T> > store;
    IndexMap index;
};

template<typename T, typename IndexMap = identity_property_map>
class unchecked_vector_property_map
    : public boost::put_get_helper<
                typename std::iterator_traits<
                    typename std::vector<T>::iterator >::reference,
                unchecked_vector_property_map<T, IndexMap> >
{
public:
    typedef typename property_traits<IndexMap>::key_type  key_type;
    typedef T value_type;
    typedef typename std::iterator_traits<
        typename std::vector<T>::iterator >::reference reference;
    typedef boost::lvalue_property_map_tag category;

    typedef checked_vector_property_map<T, IndexMap> checked_t;

    unchecked_vector_property_map(const checked_t& checked = checked_t(),
                                  size_t size = 0)
        : _checked(checked)
    {
        if (size > 0 && _checked.store->size() < size)
            _checked.store->resize(size);
    }

    unchecked_vector_property_map(const IndexMap& index_map,
                                  size_t size = 0)
    {
        *this = unchecked_vector_property_map(checked_t(index_map), size);
    }

    void reserve(size_t size) const { _checked.reserve(size); }

    reference operator[](const key_type& v) const __attribute__((always_inline))
    {
        typename property_traits<IndexMap>::value_type i =
            get(_checked.index, v);
        return (*_checked.store)[i];
    }

    std::vector<T>& get_storage() const { return _checked.get_storage(); }

    checked_t get_checked() {return _checked;}

    // deep copy
    unchecked_vector_property_map copy() const
    {
        unchecked_vector_property_map pmap(_checked.index,
                                           _checked.store->size());
        *(pmap._checked.store) = *(_checked.store);
        return pmap;
    }

private:
    checked_t _checked;
};

template<typename T, typename IndexMap>
checked_vector_property_map<T, IndexMap>
make_checked_vector_property_map(IndexMap index)
{
    return checked_vector_property_map<T, IndexMap>(index);
}

template<typename T, typename IndexMap>
unchecked_vector_property_map<T, IndexMap>
make_unchecked_vector_property_map(IndexMap index)
{
    return unchecked_vector_property_map<T, IndexMap>(index);
}


template <class Type, class Index>
unchecked_vector_property_map<Type, Index>
get_unchecked(checked_vector_property_map<Type, Index> prop)
{
    return prop.get_unchecked();
}

template <class Prop>
Prop
get_unchecked(Prop prop)
{
    return prop;
}

template <class Type, class Index>
checked_vector_property_map<Type, Index>
get_checked(unchecked_vector_property_map<Type, Index> prop)
{
    return prop.get_checked();
}

template <class Prop>
Prop
get_checked(Prop prop)
{
    return prop;
}


}

#endif
