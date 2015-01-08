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

#ifndef GRAPH_LATTICE_HH
#define GRAPH_LATTICE_HH

#include <iostream>
#include "graph_util.hh"

namespace graph_tool
{
using namespace std;
using namespace boost;


void get_pos(size_t i, const vector<size_t>& shape, vector<int>& pos)
{
    size_t offset = 1;
    for (size_t j = 0; j < shape.size(); ++j)
    {
        size_t L = shape[j];
        pos[j] = ((i / offset) % L);
        offset *= L;
    }
}

size_t get_idx(vector<int>& pos, const vector<size_t>& shape)
{
    size_t offset = 1, idx = 0;
    for (size_t j = 0; j < shape.size(); ++j)
    {
        idx += pos[j] * offset;
        offset *= shape[j];
    }
    return idx;
}

void periodic(int& x, size_t size)
{
    if (x >= int(size))
        x -= size;
    if (x < 0)
        x += size;
}

struct get_lattice
{
    template <class Graph>
    void operator()(Graph& g, vector<size_t>& shape,
                    bool periodic_boundary) const
    {
        int N = 1;
        for (size_t i = 0; i < shape.size(); ++i)
            N *= shape[i];
        for (int i = 0; i < N; ++i)
            add_vertex(g);

        vector<int> pos(shape.size());
        //#pragma omp parallel for default(shared) private(i) 
        //    firstprivate(pos) schedule(runtime) if (N > 100)
        for (int i = 0; i < N; ++i)
        {
            get_pos(i, shape, pos);
            for (size_t j = 0; j < shape.size(); ++j)
            {
                for (int k = -1; k <= 1; k += 2)
                {
                    pos[j] += k;
                    if (periodic_boundary)
                        periodic(pos[j], shape[j]);
                    if (pos[j] > 0 && size_t(pos[j]) < shape[j])
                    {
                        int m = get_idx(pos, shape);
                        if (m > i)
                        {
                            #pragma omp critical
                            add_edge(vertex(i, g), vertex(m, g), g);
                        }
                    }
                    pos[j] -= k;
                    if (periodic_boundary)
                        periodic(pos[j], shape[j]);
                }
            }
        }
    }
};

} // namespace graph_tool

#endif // GRAPH_LATTICE_HH
