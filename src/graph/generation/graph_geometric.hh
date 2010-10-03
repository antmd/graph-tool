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

#ifndef GRAPH_GEOMETRIC_HH
#define GRAPH_GEOMETRIC_HH

#include <iostream>
#include <tr1/unordered_map>
#include <boost/functional/hash.hpp>
#include "graph_util.hh"

namespace graph_tool
{
using namespace std;
using namespace boost;


void get_box(const vector<double>& p, double w, vector<int>& box)
{
    if (box.size() != p.size())
        box.resize(p.size());
    for (size_t i = 0; i < p.size(); ++i)
        box[i] = int(p[i] / w);
}

template <class Point>
double get_dist(const Point& p1, const Point& p2)
{
    double r = 0;
    for (size_t i = 0; i < p1.size(); ++i)
        r += pow(p1[i] - p2[i], 2);
    return sqrt(r);
}

void periodic(int& x, int size)
{
    if (x >= int(size))
        x -= size;
    if (x < 0)
        x += size;
}


bool is_boundary(vector<double>& point, double r,
                 const vector<pair<double, double> > ranges)
{
    bool boundary = 0;
    for (size_t j = 0; j < ranges.size(); ++j)
    {
        if (point[j] - r < ranges[j].first)
        {
            point[j] += ranges[j].second - ranges[j].first;
            boundary = true;
        }
        if (point[j] + r >= ranges[j].second)
        {
            point[j] -= ranges[j].second - ranges[j].first;
            boundary = true;
        }
    }
    return boundary;
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

        tr1::unordered_multimap<vector<int>,
                                typename graph_traits<Graph>::vertex_descriptor,
                                boost::hash<vector<int> > > boxes;

        for (int i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = add_vertex(g);
            pos[v].resize(points[i].size());
            copy(points[i].begin(), points[i].end(), pos[v].begin());

            get_box(points[i], r, box);
            boxes.insert(make_pair(box, v));

            if (periodic_boundary)
            {
                //insert the same vertex at the opposite box
                vector<double> point = points[i];
                if (is_boundary(point, r, ranges))
                {
                    get_box(point, r, box);
                    boxes.insert(make_pair(box, v));
                }
            }
        }

        #pragma omp parallel for default(shared) private(i, box) \
            schedule(dynamic)
        for (int i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            get_box(points[i], r, box);

            for (size_t k = 0; k < pow(3, box.size()); ++k)
            {
                for (size_t j = 0; j < box.size(); ++j)
                    box[j] += (int(k / pow(3, j)) % 3) - 1;

                typeof(boxes.begin()) iter, end;
                for (tie(iter, end) = boxes.equal_range(box);
                     iter != end; ++iter)
                {
                    typename graph_traits<Graph>::vertex_descriptor w =
                        iter->second;
                    double d = get_dist(pos[v], pos[w]);

                    if (periodic_boundary)
                    {
                        if (is_adjacent(v, w, g))
                            continue;
                        //get correct distance for boundary
                        vector<double> point = pos[w];
                        if (is_boundary(point, r, ranges))
                            d = min(d, get_dist(pos[v], point));
                    }

                    if (w > v && d <= r)
                    {
                        #pragma omp critical
                        add_edge(v, w, g);
                    }
                }

                for (size_t j = 0; j < box.size(); ++j)
                    box[j] -= (int(k / pow(3, j)) % 3) - 1;
            }
        }
    }
};

} // namespace graph_tool

#endif // GRAPH_GEOMETRIC_HH
