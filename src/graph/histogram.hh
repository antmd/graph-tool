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
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#ifndef HISTOGRAM_HH
#define HISTOGRAM_HH

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

namespace boost{
template<typename T> void hash_combine(size_t & seed, T const & v);
//==============================================================================
// tuple hash function
//==============================================================================
inline std::size_t hash_value(const boost::tuple<double,double,double>& hist)
{
    std::size_t seed = 0;
    boost::hash_combine(seed, boost::get<0>(hist));
    boost::hash_combine(seed, boost::get<1>(hist));
    boost::hash_combine(seed, boost::get<2>(hist));
    return seed;
}
}

#include <boost/functional/hash.hpp>


namespace std 
{
    template <class T1, class T2>
    std::ostream & operator<<(std::ostream &o, const std::pair<T1,T2> &p)
    {
        o << p.first << " \t" << p.second;
        return o;
    }

    template <class T1, class T2>
    std::istream & operator>>(std::istream &i, std::pair<T1,T2> &p)
    {
        i >> p.first >> p.second;
        return i;
    }
    
    template <class T>
    std::ostream & operator<<(std::ostream &o, const std::vector<T> &v)
    {
        for(size_t i = 0; i < v.size(); ++i)
        {
            o << v[i];
            if (i != v.size() - 1)
                o << " \t";
        }
        return o;
    }

    template <class T>
    std::istream & operator>>(std::istream &s, std::vector<T> &v)
    {
        for(size_t i = 0; i < v.size(); ++i)
            s >> v[i];
        return s;
    }

}

//==============================================================================
// PrintHistogram()
// Prints a historgram to file or stream
//==============================================================================

template <class Histogram>
void PrintHistogram (const Histogram &m, std::ostream &stream)
{
    std::vector<typename Histogram::key_type> keys;
    for (typeof(m.begin()) iter = m.begin(); iter != m.end(); iter++)
        keys.push_back(iter->first);
    sort(keys.begin(), keys.end());
    stream.precision(20);
    stream.setf(std::ios_base::scientific);
    for (size_t i = 0; i < keys.size(); ++i)
    {
        typename Histogram::const_iterator val = m.find(keys[i]);
        stream << *val << std::endl;
    }
}

template <class Histogram> 
void PrintHistogram (const Histogram &m, std::string output_file)
{
    
    boost::trim(output_file);
    std::ofstream file(output_file.c_str(), std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
    file.exceptions(std::ios_base::badbit | std::ios_base::failbit);
    boost::iostreams::filtering_stream<boost::iostreams::output> stream;
    if (boost::ends_with(output_file,".gz"))
        stream.push(boost::iostreams::gzip_compressor());
    if (boost::ends_with(output_file,".bz2"))
        stream.push(boost::iostreams::bzip2_compressor());
    stream.push(file);
    stream.exceptions(std::ios_base::badbit);
    PrintHistogram(m, stream);
}

template <class Histogram>
void ReadHistogram(Histogram &m, std::istream &stream)
{
    typename Histogram::key_type key;

    while (stream && !stream.eof())
    {
        stream >> key;
        if (stream.eof() || !stream)
            continue;
        stream >> m[key];
    }
}

template <class Histogram>
void ReadHistogram(Histogram &m, std::string input_file)
{
    boost::trim(input_file);
    std::ifstream file(input_file.c_str(), std::ios_base::in | std::ios_base::binary);
    file.exceptions(std::ios_base::badbit | std::ios_base::failbit );
    boost::iostreams::filtering_stream<boost::iostreams::input> stream;
    if (boost::ends_with(input_file,".gz"))
        stream.push(boost::iostreams::gzip_decompressor());
    if (boost::ends_with(input_file,".bz2"))
        stream.push(boost::iostreams::bzip2_decompressor());
    stream.push(file);
    stream.exceptions(std::ios_base::badbit);
    ReadHistogram(m, stream);
}

//==============================================================================
// GetHistogramMean()
// Gets the mean value of a histogram
//==============================================================================
template <class Histogram> 
double GetHistogramMean (const Histogram &m)
{
    int total = 0;
    double mean = 0;
    for (typeof(m.begin()) iter = m.begin(); iter != m.end(); iter++)
    {
        mean += double(iter->first * iter->second);
        total += iter->second;
    }
    
    return (total > 0)?mean/total:0.0;
}

//==============================================================================
// GetHistogramDeviation()
// Gets the standard deviation of a histogram
//==============================================================================
template <class Histogram> 
double GetHistogramDeviation (const Histogram &m, double avg)
{  
    double dev = 0.0;
    int total = 0;
    for (typeof(m.begin()) iter = m.begin(); iter != m.end(); iter++)
    {
        dev += double( (iter->first - avg) * (iter->first - avg) * iter->second);
        total += iter->second;
    }
    return (total > 1)?sqrt(dev/(total-1)):0.0;
}


#endif //HISTOGRAM_HH
