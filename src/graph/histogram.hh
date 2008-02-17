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

#ifndef HISTOGRAM_HH
#define HISTOGRAM_HH

#include <tr1/unordered_map>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

namespace boost{
template<typename T> void hash_combine(size_t & seed, T const & v);

//
// tuple hash function
//

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

// histogram types
typedef std::tr1::unordered_map<double,double> hist_t;
typedef std::tr1::unordered_map<std::pair<double,double>,double,
                                boost::hash<std::pair<double,double> > >
    hist2d_t;
typedef std::tr1::unordered_map<boost::tuple<double,double,double>,double,
                                boost::hash<boost::tuple<double,
                                                         double,double> > >
    hist3d_t;
typedef std::tr1::unordered_map<double,std::pair<double,double> > avg_corr_t;


// gets the mean value of a histogram
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

// gets the standard deviation of a histogram
template <class Histogram>
double GetHistogramDeviation (const Histogram &m, double avg)
{
    double dev = 0.0;
    int total = 0;
    for (typeof(m.begin()) iter = m.begin(); iter != m.end(); iter++)
    {
        dev += double( (iter->first - avg) *
                       (iter->first - avg) * iter->second);
        total += iter->second;
    }
    return (total > 1)?sqrt(dev/(total-1)):0.0;
}

#endif //HISTOGRAM_HH
