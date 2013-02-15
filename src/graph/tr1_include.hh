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

#ifndef TR1_INCLUDE_HH
#define TR1_INCLUDE_HH

// include tr1 from libstdc++ only if it is not _really_ old. Otherwise use
// boost.

#if defined(__GLIBCXX__) && __GLIBCXX__ > 20070719
#   define TR1_HEADER(header) <tr1/header>
#else
#   define BOOST_TR1_USE_OLD_TUPLE
#   define TR1_HEADER(header) <boost/tr1/header.hpp>
#endif

#endif
