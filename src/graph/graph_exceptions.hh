// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007-2011 Tiago de Paula Peixoto <tiago@skewed.de>
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

#ifndef GRAPH_EXCEPTIONS_HH
#define GRAPH_EXCEPTIONS_HH
#include "config.h"
#include <string>

// Exceptions
// ==========
//
// This is the main exception which will be thrown the outside world, when
// things go wrong

namespace graph_tool{
using namespace std;

class GraphException : public std::exception
{
    string _error;
public:
    GraphException(const string& error) {_error = error;}
    virtual ~GraphException() throw () {}
    virtual const char * what () const throw () {return _error.c_str();}
protected:
    virtual void SetError(const string& error) {_error = error;}
};

class IOException : public GraphException
{
public:
    IOException(const string& error): GraphException(error) {}
};

class ValueException : public GraphException
{
public:
    ValueException(const string& error): GraphException(error) {}
};

} // namespace std

#endif // GRAPH_EXCEPTIONS_HH
