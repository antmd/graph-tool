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

#ifndef GRAPH_AVERAGE_HH
#define GRAPH_AVERAGE_HH

#include <algorithm>
#include <boost/numeric/conversion/bounds.hpp>
#include <boost/numeric/conversion/cast.hpp>

namespace graph_tool
{
using namespace std;
using namespace boost;

class VertexAverageTraverse
{
public:
    template <class Graph, class DegreeSelector, class ValueType>
    void operator()(Graph& g, typename graph_traits<Graph>::vertex_descriptor v,
                    DegreeSelector& deg, ValueType& a, ValueType& aa,
                    size_t& count)
    {
        ValueType x = deg(v, g);
        a += x;
        aa += x*x;
        count++;
    }
};

class EdgeAverageTraverse
{
public:
    template <class Graph, class EdgeProperty, class ValueType>
    void operator()(Graph& g, typename graph_traits<Graph>::vertex_descriptor v,
                    EdgeProperty& eprop, ValueType& a, ValueType& aa,
                    size_t& count)
    {
        typename graph_traits<Graph>::out_edge_iterator e, e_begin, e_end;
        tie(e_begin,e_end) = out_edges(v,g);
        for(e = e_begin; e != e_end; ++e)
        {
            ValueType x = eprop[*e];
            a += x;
            aa += x*x;
            count++;
        }
    }
};

// generalized functor to obtain average of different types of "degrees"
template <class AverageTraverse>
struct get_average
{
    get_average(long double& a, long double& dev)
        : _a(a), _dev(dev) {}

    template <class Graph, class DegreeSelector>
    void operator()(Graph& g, DegreeSelector deg) const
    {
        typedef typename DegreeSelector::value_type value_type;
        long double a = 0, aa = 0;
        size_t count = 0;

        AverageTraverse traverse;
        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) \
            reduction(+:a,aa,count) schedule(runtime) if (N > 100)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            traverse(g, v, deg, a, aa, count);
        }

        _a = a/count;
        _dev = sqrt((aa/count - _a*_a))/sqrt(count);
    }

    long double& _a;
    long double& _dev;
};

} // graph_tool namespace

#endif // GRAPH_AVERAGE_HH

