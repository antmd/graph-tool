// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007-2010 Tiago de Paula Peixoto <tiago@forked.de>
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
#include <boost/array.hpp>
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
    typedef boost::array<ValueType,Dim> point_t; // point type to be
                                                 // histogrammed
    typedef boost::array<size_t,Dim> bin_t;      // bin type

    typedef boost::multi_array<CountType,Dim> count_t; // the histogram itself

    typedef boost::mpl::int_<Dim> dim;
    typedef CountType count_type;
    typedef ValueType value_type;

    // floating point type to calculate the mean
    typedef typename boost::mpl::if_<boost::is_floating_point<ValueType>,
                                     ValueType, double>::type mean_t;

    Histogram(const boost::array<std::vector<ValueType>, Dim>& bins,
              const boost::array<std::pair<ValueType,ValueType>,Dim>&
              data_range)
        : _bins(bins), _data_range(data_range)
    {
        bin_t new_shape;
        for (size_t j = 0; j < Dim; ++j)
            new_shape[j] = _bins[j].size();
        _counts.resize(new_shape);
    }

    void PutValue(const point_t& v, const CountType& weight = 1)
    {
        bin_t bin;
        for (size_t i = 0; i < Dim; ++i)
        {
            if (_bins[i].size() == 1) // constant bin width
            {
                bin[i] = (v[i] - _data_range[i].first)/_bins[i][0];
                if (bin[i] >= _counts.shape()[i])
                {
                    boost::array<size_t, Dim> new_shape;
                    for (size_t j = 0; j < Dim; ++j)
                        new_shape[j] = _counts.shape()[j];
                    new_shape[i] = bin[i] + 1;
                    _counts.resize(new_shape);
                    if (v[i] > _data_range[i].second)
                        _data_range[i].second = v[i];
                }
            }
            else // arbitrary bins. do a binary search
            {
                std::vector<ValueType>& bins = _bins[i];
                typeof(bins.begin()) iter = upper_bound(bins.begin(),
                                                        bins.end(), v[i]);
                if (iter == bins.end()) // larger than any bin, thus belongs to
                {                       // the last one
                    bin[i] = bins.size() - 1;
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

    boost::array<std::pair<ValueType,ValueType>,Dim>& GetDataRange()
    { return _data_range; }


    boost::array<std::vector<ValueType>, Dim> GetBins()
    {
        boost::array<std::vector<ValueType>, Dim> bins;
        for (size_t j = 0; j < Dim; ++j)
            if (_bins[j].size() == 1) // constant bin width
            {
                for (ValueType i = _data_range[j].first;
                     i <= _data_range[j].second; i += _bins[j][0])
                    bins[j].push_back(i);
            }
            else
            {
                bins[j] = _bins[j];
            }
        return bins;
    }

protected:
    boost::multi_array<CountType,Dim> _counts;
    boost::array<std::vector<ValueType>, Dim> _bins;
    boost::array<std::pair<ValueType,ValueType>,Dim> _data_range;
};


// This class will encapsulate a histogram, and atomically sum it to a given
// resulting histogram (which is shared among all copies) after it is
// destructed, or when the Gather() member function is called. This enables, for
// instance, a histogram to be built in parallel.

template <class Histogram>
class SharedHistogram: public Histogram
{
public:

    SharedHistogram(Histogram &hist): Histogram(hist), _sum(&hist) {}
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
                typename Histogram::bin_t shape;
                for (size_t i = 0; i <  this->_counts.num_dimensions(); ++i)
                    shape[i] = max(this->_counts.shape()[i],
                                   _sum->GetArray().shape()[i]);
                _sum->GetArray().resize(shape);
                for (size_t i = 0; i < this->_counts.num_elements(); ++i)
                    _sum->GetArray().data()[i] += this->_counts.data()[i];
                for (size_t i = 0; i < size_t(Histogram::dim::value); ++i)
                {
                    _sum->GetDataRange()[i].first =
                        min(this->_data_range[i].first,
                            _sum->GetDataRange()[i].first);
                    _sum->GetDataRange()[i].second =
                        max(this->_data_range[i].second,
                            _sum->GetDataRange()[i].second);
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
