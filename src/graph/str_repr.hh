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

#ifndef STR_REPR_HH
#define STR_REPR_HH

#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

//
// Data type string representation
// ===============================
//
// String representation of individual data types. Among other things, we have
// to take care specifically that no information is lost with floating point
// I/O.
//

namespace boost
{

//
// "chars" should be printed as numbers, since they can be non-printable
//

template <>
string lexical_cast<string,uint8_t>(const uint8_t& val)
{
    return lexical_cast<std::string>(int(val));
}

template <>
uint8_t lexical_cast<uint8_t,string>(const string& val)
{
    return uint8_t(lexical_cast<int>(val));
}

// float, double and long double should be printed in hexadecimal format to
// preserve internal representation

template <>
string lexical_cast<string,float>(const float& val)
{
    char* str = 0;
    int retval = asprintf(&str, "%a", val);
    if (retval == -1)
        throw bad_lexical_cast();
    std::string ret = str;
    free(str);
    return ret;
}

template <>
float lexical_cast<float,string>(const string& val)
{
    float ret;
    int nc = sscanf(val.c_str(), "%a", &ret);
    if (nc != 1)
        throw bad_lexical_cast();
    return ret;
}

template <>
string lexical_cast<string,double>(const double& val)
{
    char* str = 0;
    int retval = asprintf(&str, "%la", val);
    if (retval == -1)
        throw bad_lexical_cast();
    std::string ret = str;
    free(str);
    return ret;
}

template <>
double lexical_cast<double,string>(const string& val)
{
    double ret;
    int nc = sscanf(val.c_str(), "%la", &ret);
    if (nc != 1)
        throw bad_lexical_cast();
    return ret;
}

template <>
string lexical_cast<string,long double>(const long double& val)
{
    char* str = 0;
    int retval = asprintf(&str, "%La", val);
    if (retval == -1)
        throw bad_lexical_cast();
    std::string ret = str;
    free(str);
    return ret;
}

template <>
long double lexical_cast<long double,string>(const string& val)
{
    long double ret;
    int nc = sscanf(val.c_str(), "%La", &ret);
    if (nc != 1)
        throw bad_lexical_cast();
    return ret;
}
} // namespace boost

//
// stream i/o of std::vector<>
//

namespace std
{

// string vectors need special attention, since separators must be properly
// escaped.
template <>
ostream& operator<<(ostream& out, const vector<string>& vec)
{
    for (size_t i = 0; i < vec.size(); ++i)
    {
        string s = vec[i];
        // escape separators
        boost::replace_all(s, "\\", "\\\\");
        boost::replace_all(s, ", ", ",\\ ");

        out << s;
        if (i < vec.size() - 1)
            out << ", ";
    }
    return out;
}

template <>
istream& operator>>(istream& in, vector<string>& vec)
{
    using namespace boost;
    using namespace boost::algorithm;
    using namespace boost::xpressive;

    vec.clear();
    string data;
    while (in.good())
    {
        string line;
        getline(in, line);
        data += line;
    }

    if (data == "")
        return in; // empty string is OK

    sregex re = sregex::compile(", ");
    sregex_token_iterator iter(data.begin(), data.end(), re, -1), end;
    for (; iter != end; ++iter)
    {
        vec.push_back(*iter);
        // un-escape separators
        boost::replace_all(vec.back(), ",\\ ", ", ");
        boost::replace_all(vec.back(), "\\\\", "\\");
    }
    return in;
}

} // std namespace

#endif // STR_REPR_HH
