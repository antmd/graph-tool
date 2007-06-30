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

#ifndef SHARED_MAP_HH
#define SHARED_MAP_HH

#include <tr1/unordered_set>

//==============================================================================
// SharedMap
// This will encapsulate a map, and atomically sum it to a given resulting map 
// (which is shared among all copies) after it is destructed, or when the 
// Gather() member function is called. This enables, for instance, a histogram
// to built in parallel.
//==============================================================================

template <class Map>
class SharedMap: public Map
{
public:
    SharedMap(Map &map):_sum(&map) {}
    ~SharedMap()
    {
        Gather();
    }

    void Gather()
    {
        if (_sum != 0)
        {
            for (typeof(this->begin()) iter = this->begin(); iter != this->end(); ++iter)
            {
                #pragma omp critical
                {
                    (*_sum)[iter->first] += iter->second;
                }
            }
            _sum = 0;
        }
    }
private:
    Map* _sum;
};


#endif //SHARED_MAP_HH
