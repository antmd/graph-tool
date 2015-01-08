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

#ifndef HISTOGRAM_HH
#define HISTOGRAM_HH

#include <vector>
#include <utility>
#include <algorithm>
#include <array>
#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>

#include <boost/type_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/int.hpp>

#include <boost/python/object.hpp>

//
// This is a generic multidimensional histogram type
//

template <class ValueType, class CountType, size_t Dim>
class Histogram
{
public:
    typedef std::array<ValueType,Dim> point_t; // point type to be
                                               // histogrammed
    typedef std::array<size_t,Dim> bin_t;      // bin type

    typedef boost::multi_array<CountType,Dim> count_t; // the histogram itself

    typedef boost::mpl::int_<Dim> dim;
    typedef CountType count_type;
    typedef ValueType value_type;

    // floating point type to calculate the mean
    typedef typename boost::mpl::if_<boost::is_floating_point<ValueType>,
                                     ValueType, double>::type mean_t;

    Histogram(const std::array<std::vector<ValueType>, Dim>& bins):
        _bins(bins)
    {
        bin_t new_shape;
        for (size_t j = 0; j < Dim; ++j)
        {
            if (_bins[j].size() < 1)
                throw std::range_error("invalid bin edge number < 1!");

            _data_range[j] = std::make_pair(0, 0);
            value_type delta = _bins[j][1] - _bins[j][0];

            if (_bins[j].size() == 2)
            {
                _data_range[j] = std::make_pair(_bins[j][0], _bins[j][0]);
                delta = _bins[j][1];
                _const_width[j] = true;
            }
            else
            {
                // detect whether the given bins are of constant width, for faster
                // binning
                _const_width[j] = true;
                for (size_t i = 2; i < _bins[j].size(); ++i)
                {
                    value_type d = _bins[j][i] - _bins[j][i-1];
                    if (delta != d)
                        _const_width[j] = false;
                }

                if (_const_width[j])
                    _data_range[j] = std::make_pair(_bins[j].front(),
                                                    _bins[j].back());
            }
            if (delta == 0)
                throw std::range_error("invalid bin size of zero!");

            new_shape[j] = _bins[j].size() - 1;
        }
        _counts.resize(new_shape);
    }

    void PutValue(const point_t& v, const CountType& weight = 1)
    {
        bin_t bin;
        for (size_t i = 0; i < Dim; ++i)
        {
            if (_const_width[i])
            {
                value_type delta;

                if (_data_range[i].first == _data_range[i].second)
                {
                    delta = _bins[i][1];

                    if (v[i] < _data_range[i].first)
                        return; // out of bounds
                }
                else
                {
                    delta = _bins[i][1] - _bins[i][0];

                    if (v[i] < _data_range[i].first ||
                        v[i] >= _data_range[i].second)
                        return; // out of bounds
                }

                bin[i] = size_t((v[i] - _data_range[i].first) / delta);
                if (bin[i] >= _counts.shape()[i]) // modify shape
                {
                    bin_t new_shape;
                    for (size_t j = 0; j < Dim; ++j)
                        new_shape[j] = _counts.shape()[j];
                    new_shape[i] = bin[i] + 1;
                    _counts.resize(new_shape);
                    while (_bins[i].size() < new_shape[i] + 1)
                        _bins[i].push_back(_bins[i].back() + delta);
                }
            }
            else // arbitrary bins widths. do a binary search
            {
                std::vector<ValueType>& bins = _bins[i];
                typeof(bins.begin()) iter = upper_bound(bins.begin(),
                                                        bins.end(), v[i]);
                if (iter == bins.end())
                {
                    return;  // falls off from last bin, do not count
                }
                else
                {
                    bin[i] = iter - bins.begin();
                    if (bin[i] == 0)
                        return; // falls off from fist bin, do not count
                    else
                       --bin[i];
                }
            }
        }
        _counts(bin) += weight;
    }

    boost::multi_array<CountType,Dim>& GetArray() { return _counts; }

    std::array<std::pair<ValueType,ValueType>,Dim>& GetDataRange()
    { return _data_range; }

    std::array<std::vector<ValueType>, Dim>& GetBins() { return _bins; }

protected:
    boost::multi_array<CountType,Dim> _counts;
    std::array<std::vector<ValueType>, Dim> _bins;
    std::array<std::pair<ValueType,ValueType>,Dim> _data_range;
    std::array<bool,Dim> _const_width;
};


// This class will encapsulate a histogram, and atomically sum it to a given
// resulting histogram (which is shared among all copies) after it is
// destructed, or when the Gather() member function is called. This enables, for
// instance, a histogram to be built in parallel.

template <class Histogram>
class SharedHistogram: public Histogram
{
public:

    SharedHistogram(Histogram& hist): Histogram(hist), _sum(&hist) {}
    ~SharedHistogram()
    {
        Gather();
    }

    void Gather()
    {
        if (_sum != 0)
        {
            #pragma omp critical
            {
                typename Histogram::bin_t idx;

                typename Histogram::bin_t shape;
                for (size_t i = 0; i <  this->_counts.num_dimensions(); ++i)
                    shape[i] = std::max(this->_counts.shape()[i],
                                        _sum->GetArray().shape()[i]);
                _sum->GetArray().resize(shape);
                for (size_t i = 0; i < this->_counts.num_elements(); ++i)
                {
                    size_t offset = 1;
                    for (size_t j = 0; j < this->_counts.num_dimensions(); ++j)
                    {
                        size_t L = this->_counts.shape()[j];
                        idx[j] = ((i / offset) % L);
                        offset *= L;
                    }
                    _sum->GetArray()(idx) += this->_counts(idx);
                }
                for (int i = 0; i < Histogram::dim::value; ++i)
                {
                    if (_sum->GetBins()[i].size() < this->_bins[i].size())
                        _sum->GetBins()[i] = this->_bins[i];
                }
            }
            _sum = 0;
        }
    }
private:
    Histogram* _sum;
};


//
// useful functions to get the the mean and standard deviations from simple
// map-based, non-binned histograms. Not to be used with the above type.
//

// gets the mean value of a histogram
template <class Map>
double GetMapMean (const Map &m)
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
template <class Map>
double GetMapDeviation (const Map &m, double avg)
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
