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

#ifndef GRAPH_GEOMETRIC_HH
#define GRAPH_GEOMETRIC_HH

#include <iostream>
#include <unordered_map>

#include <boost/functional/hash.hpp>
#include "graph_util.hh"

#ifndef __clang__
#include <ext/numeric>
using __gnu_cxx::power;
#else
template <class Value>
Value power(Value value, int n)
{
    return pow(value, n);
}
#endif

namespace graph_tool
{
using namespace std;
using namespace boost;


template <class Point, class Range>
void get_box(const Point& p, double w, vector<int>& box,
             const Range& ranges, bool periodic)
{
    if (box.size() != p.size())
        box.resize(p.size());
    for (size_t i = 0; i < p.size(); ++i)
    {
        box[i] = int(floor(p[i] / w));
        if (periodic && p[i] == ranges[i].second)
            box[i] -= 1;
    }
}

template <class Point, class Range>
double get_dist(const Point& p1, const Point& p2,
                const Range& ranges, bool periodic)
{
    double r = 0, diff, size;
    for (size_t i = 0; i < p1.size(); ++i)
    {
        diff = abs(p1[i] - p2[i]);
        if (periodic)
        {
            size = abs(ranges[i].second - ranges[i].first);
            diff = min(diff, abs(diff - size));
        }
        r += pow(diff, 2);
    }
    return sqrt(r);
}

void periodic(vector<int>& box, const vector<pair<int,int> >& ranges)
{
    for (size_t i = 0; i < box.size(); ++i)
    {
        if (box[i] >= ranges[i].second)
            box[i] = ranges[i].first;
        if (box[i] < ranges[i].first)
            box[i] = ranges[i].second - 1;
    }
}


template <class Graph>
bool is_adjacent(typename graph_traits<Graph>::vertex_descriptor v1,
                 typename graph_traits<Graph>::vertex_descriptor v2,
                 Graph& g)
{
    typename graph_traits<Graph>::out_edge_iterator e, e_end;
    for (tie(e, e_end) = out_edges(v1, g); e != e_end; ++e)
        if (target(*e, g) == v2)
            return true;
    return false;
}

struct get_geometric
{
    template <class Graph, class Pos>
    void operator()(Graph& g, Pos upos, vector<vector<double> >& points,
                    vector<pair<double, double> >& ranges,
                    double r, bool periodic_boundary) const
    {
        int N = points.size();
        typename Pos::checked_t pos = upos.get_checked();
        vector<int> box;
        vector<pair<int, int> > box_ranges;

        double w = 2 * r;
        if (periodic_boundary)
        {
            box_ranges.resize(ranges.size());
            for (size_t i = 0; i < ranges.size(); ++i)
            {
                box_ranges[i].first = floor(ranges[i].first / w);
                box_ranges[i].second = ceil(ranges[i].second / w);
            }
        }

        std::unordered_multimap<vector<int>,
                                typename graph_traits<Graph>::vertex_descriptor,
                                boost::hash<vector<int> > > boxes;

        for (int i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = add_vertex(g);
            pos[v].resize(points[i].size());
            copy(points[i].begin(), points[i].end(), pos[v].begin());
            get_box(points[i], w, box, ranges, periodic_boundary);
            boxes.insert(make_pair(box, v));
        }

        int i;
        #pragma omp parallel for default(shared) private(i, box) \
            schedule(runtime) if (N > 100)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);

            get_box(points[i], w, box, ranges, periodic_boundary);
            for (int k = 0; k < power(3, int(box.size())); ++k)
            {
                for (int j = 0; j < int(box.size()); ++j)
                  box[j] += ((k / power(3, j)) % 3) - 1;

                if (periodic_boundary)
                    periodic(box, box_ranges);

                typeof(boxes.begin()) iter, end;
                for (tie(iter, end) = boxes.equal_range(box);
                     iter != end; ++iter)
                {
                    typename graph_traits<Graph>::vertex_descriptor w =
                        iter->second;
                    double d = get_dist(pos[v], pos[w], ranges,
                                        periodic_boundary);

                    if (w > v && d <= r &&
                        (!periodic_boundary || !is_adjacent(v, w, g)))
                    {
                        #pragma omp critical
                        add_edge(v, w, g);
                    }
                }
                get_box(points[i], w, box, ranges, periodic_boundary);
            }
        }
    }
};

} // namespace graph_tool

#endif // GRAPH_GEOMETRIC_HH
