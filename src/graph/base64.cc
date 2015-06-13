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


#include "base64.hh"

#include <sstream>
#include <boost/archive/iterators/ostream_iterator.hpp>
#include <boost/archive/iterators/base64_from_binary.hpp>
#include <boost/archive/iterators/binary_from_base64.hpp>
#include <boost/archive/iterators/transform_width.hpp>

std::string base64_encode(const std::string& s)
{
    static const std::string base64_padding[] = {"", "==","="};
    namespace bai = boost::archive::iterators;
    std::stringstream os;
    typedef bai::base64_from_binary<bai::transform_width<const char *, 6, 8> >
        base64_enc;
    std::copy(base64_enc(s.c_str()), base64_enc(s.c_str() + s.size()),
              std::ostream_iterator<char>(os));
    os << base64_padding[s.size() % 3];
    return os.str();
}

std::string base64_decode(const std::string& s)
{
    namespace bai = boost::archive::iterators;
    std::stringstream os;

    typedef bai::transform_width<bai::binary_from_base64<const char *>, 8, 6> base64_dec;

    unsigned int size = s.size();

    // Remove the padding characters, cf. https://svn.boost.org/trac/boost/ticket/5629
    if (size && s[size - 1] == '=')
    {
        --size;
        if (size && s[size - 1] == '=')
            --size;
    }
    if (size == 0)
        return std::string();

    std::copy(base64_dec(s.data()), base64_dec(s.data() + size),
              std::ostream_iterator<char>(os));

    return os.str();
}
