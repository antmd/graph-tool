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

#include "graph.hh"

#include "graph_filtering.hh"

#include "graph_selectors.hh"
#include "graph_properties.hh"

#include <boost/mpl/quote.hpp>

#include <cmath>

using namespace std;
using namespace boost;
using namespace graph_tool;

typedef pair<double, double> point_t;

point_t interpolate(const point_t& p1, const point_t& p2, double r = 0.5)
{
    point_t ret;
    ret.first = (1 - r) * p1.first + p2.first * r;
    ret.second = (1 - r) * p1.second + p2.second * r;
    return ret;
}

void to_bezier(const vector<point_t> &x, vector<point_t>& ncp)
{
    vector<point_t> cp(x.size() + 6);
    for (size_t i = 0; i < 3; ++i)
        cp[i] = x[0];
    for (size_t i = 0; i < x.size(); ++i)
        cp[i + 3] = x[i];
    for (size_t i = cp.size() - 3; i < cp.size(); ++i)
        cp[i] = x.back();

    vector<point_t> one_thirds(cp.size() - 1);
    vector<point_t> two_thirds(cp.size() - 1);

    for (size_t i = 0; i < cp.size() - 1; ++i)
    {
        const point_t& p1 = cp[i];
        const point_t& p2 = cp[i + 1];
        one_thirds[i] = interpolate(p1, p2, 1./3);
        two_thirds[i] = interpolate(p2, p1, 1./3);
    }

    ncp.resize((cp.size() - 3) * 3);
    for (size_t i = 0; i < cp.size() - 3; ++i)
    {
        size_t pos = i * 3;
        ncp[pos] = one_thirds[i + 1];
        ncp[pos + 1] = two_thirds[i + 1];
        ncp[pos + 2] = interpolate(two_thirds[i + 1], one_thirds[i + 2]);
    }
}

void transform(vector<point_t>& cp)
{
    point_t origin = cp[0];
    for (size_t i = 0; i < cp.size(); ++i)
    {
        cp[i].first -= origin.first;
        cp[i].second -= origin.second;
    }

    double t = atan2(cp.back().second - cp.front().second,
                     cp.back().first - cp.front().first);

    for (size_t i = 0; i < cp.size(); ++i)
    {
        double x = cp[i].first;
        double y = cp[i].second;

        cp[i].first = cos(t) * x + sin(t) * y;
        cp[i].second = -sin(t) * x + cos(t) * y;
    }

    point_t d;
    d.first = cp.back().first - cp.front().first;
    d.second = cp.back().second - cp.front().second;
    double r = sqrt(d.first * d.first + d.second * d.second);

    for (size_t i = 0; i < cp.size(); ++i)
        cp[i].first /= r;
}

template <class PosProp>
void get_control_points(vector<size_t>& path, PosProp pos, double beta,
                        vector<point_t>& ncp)
{
    size_t L = path.size();
    vector<point_t> cp(L);
    for (size_t i = 0; i < L; ++i)
        cp[i] = make_pair(double(pos[path[i]][0]),
                          double(pos[path[i]][1]));
    ncp.resize(L);
    for (size_t i = 0; i < L; ++i)
    {
        ncp[i].first = beta * cp[i].first + (1 - beta) * (cp[0].first + (cp.back().first - cp[0].first) * i / (L - 1.));
        ncp[i].second = beta * cp[i].second + (1 - beta) * (cp[0].second + (cp.back().second - cp[0].second) * i / (L - 1.));
    }
}

template <class Graph>
void tree_path(Graph& g, size_t s, size_t t, vector<size_t>& path)
{
    vector<size_t> s_root;
    vector<size_t> t_root;
    s_root.push_back(s);
    t_root.push_back(t);

    size_t v = s;
    size_t u = t;

    while (v != u)
    {
        typename graph_traits<Graph>::in_edge_iterator e, e_end;
        tie(e, e_end) = in_edges(v, g);
        if (e == e_end)
            throw GraphException("Invalid hierarchical tree: No path from source to target.");
        v = source(*e, g);
        s_root.push_back(v);

        tie(e, e_end) = in_edges(u, g);
        if (e == e_end)
            throw GraphException("Invalid hierarchical tree: No path from source to target.");
        u = source(*e, g);
        if (u != v)
            t_root.push_back(u);
    }
    path = s_root;
    for (typeof(t_root.rbegin()) iter = t_root.rbegin();
         iter != t_root.rend(); ++iter)
        path.push_back(*iter);
}


template<class T>
void pack(vector<point_t>& cp, vector<T>& ncp)
{
    ncp.resize(cp.size() * 2);
    for (size_t i = 0; i < cp.size(); ++i)
    {
        ncp[2 * i] = cp[i].first;
        ncp[2 * i + 1] = cp[i].second;
    }
}

struct do_get_cts
{
    template <class Graph, class Tree, class PosProp, class CMap>
    void operator()(Graph& g, Tree* t, PosProp tpos, double beta, CMap cts) const
    {
        vector<size_t> path;
        vector<point_t> cp;
        vector<point_t> ncp;

        typename graph_traits<Graph>::edge_iterator e, e_end;
        for(tie(e, e_end) = edges(g); e != e_end; ++e)
        {
            typename graph_traits<Graph>::vertex_descriptor u, v;
            u = source(*e, g);
            v = target(*e, g);
            if (u == v)
                continue;
            path.clear();
            tree_path(*t, u, v, path);
            cp.clear();
            get_control_points(path, tpos, beta, cp);
            ncp.clear();
            to_bezier(cp, ncp);
            transform(ncp);
            pack(ncp, cts[*e]);
        }
    }
};

struct get_pointers
{
    template <class List>
    struct apply
    {
        typedef typename boost::mpl::transform<List,
                                               boost::mpl::quote1<std::add_pointer> >::type type;
    };
};

void get_cts(GraphInterface& gi, GraphInterface& tgi,
             boost::any otpos, double beta, boost::any octs)
{
    typedef property_map_type::apply<vector<double>,
                                     GraphInterface::edge_index_map_t>::type
        eprop_t;

    eprop_t cts = boost::any_cast<eprop_t>(octs);


    run_action<graph_tool::detail::always_directed, boost::mpl::true_>()
        (gi, std::bind(do_get_cts(), placeholders::_1, placeholders::_2,
                       placeholders::_3, beta, cts),
         get_pointers::apply<graph_tool::detail::always_directed>::type(),
         vertex_scalar_vector_properties())
        (tgi.GetGraphView(), otpos);
}
