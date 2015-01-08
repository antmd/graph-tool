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

#include "graph_filtering.hh"
#include "graph.hh"
#include "graph_properties.hh"

#include <boost/bind.hpp>

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

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/hypot.hpp>
#include <boost/math/constants/constants.hpp>

#include <boost/graph/fruchterman_reingold.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

namespace graph_tool
{
// convert point types
template <class Val>
struct convert<vector<Val>, typename convex_topology<2>::point>
{
    vector<Val> operator()(const typename convex_topology<2>::point& p) const
    {
        vector<Val> v(2);
        for (size_t i = 0; i < 2; ++i)
            v[i] = p[i];
        return v;
    }
};

template <class Val>
struct convert<typename convex_topology<2>::point, vector<Val> >
{
    typename convex_topology<2>::point operator()(const vector<Val>& v) const
    {
        typename convex_topology<2>::point p;
        for (size_t i = 0; i < min(size_t(2), v.size()); ++i)
            p[i] = v[i];
        return p;
    }
};
} // graph_tool namespace

template<class T>
struct anneal_cooling
{
    typedef T result_type;

    anneal_cooling(T ti, T tf, std::size_t iterations)
        : _ti(ti), _tf(tf), _iter(0), _n_iter(iterations) 
    {
        _beta = (log(_tf) - log(_ti)) / _n_iter;
    }

    T operator()()
    {
        T temp = _ti * exp(T(_iter) * _beta);
        ++_iter;
        if (_iter == _n_iter)
            temp = 0;
        return temp;
    }

private:
    T _ti, _tf;
    size_t _iter;
    size_t _n_iter;
    T _beta;
};

template <class Topology>
struct get_layout
{
    template <class WeightMap, class Value>
    struct attr_force
    {
        attr_force(WeightMap w, Value a): _w(w), _a(a) {}
        WeightMap _w;
        Value _a;

        template <class Graph, class Edge, class KVal, class DVal>
        Value operator()(Edge e, KVal k, DVal dist, const Graph&) const
        {
            return _a * get(_w, e) * power(dist, 2) / k;
        }
    };

    template <class Value>
    struct rep_force
    {
        rep_force(Value r): _r(r){}
        Value _r;

        template <class Graph, class Vertex, class KVal, class DVal>
        Value operator()(Vertex v1, Vertex v2, KVal k, DVal dist,
                         const Graph&) const
        {
            return _r * power(k, 2) / dist;
        }
    };


    template <class Graph, class PosMap, class WeightMap>
    void operator()(Graph& g, PosMap pos, WeightMap weight,
                    pair<double, double> f, double scale, bool grid,
                    pair<double, double> temp, size_t n_iter) const
    {
        typedef typename property_traits<PosMap>::value_type::value_type pos_t;
        anneal_cooling<pos_t> cool(temp.first, temp.second, n_iter);
        attr_force<WeightMap, pos_t> af(weight, f.first);
        rep_force<pos_t> rf(f.second);
        Topology topology(scale);
        ConvertedPropertyMap<PosMap, typename Topology::point_type> cpos(pos);
        if (grid)
            fruchterman_reingold_force_directed_layout
                (g, cpos,
                 topology,
                 attractive_force(af).
                 repulsive_force(rf).
                 cooling(cool));
        else
            fruchterman_reingold_force_directed_layout
                (g, cpos,
                 topology,
                 attractive_force(af).
                 repulsive_force(rf).
                 cooling(cool).
                 force_pairs(all_force_pairs()));
    }


};


void fruchterman_reingold_layout(GraphInterface& g, boost::any pos,
                                 boost::any weight, double a, double r,
                                 bool square, double scale, bool grid,
                                 double ti, double tf, size_t max_iter)
{
    typedef ConstantPropertyMap<double,GraphInterface::edge_t> weight_map_t;
    typedef boost::mpl::push_back<edge_scalar_properties, weight_map_t>::type
        edge_props_t;

    if(weight.empty())
        weight = weight_map_t(1.0);
    if (square)
        run_action<graph_tool::detail::never_directed>()
            (g,
             std::bind(get_layout<square_topology<> >(), placeholders::_1,
                       placeholders::_2, placeholders::_3, make_pair(a, r), scale,
                       grid, make_pair(ti, tf), max_iter),
             vertex_floating_vector_properties(), edge_props_t())
            (pos, weight);
    else
        run_action<graph_tool::detail::never_directed>()
            (g,
             std::bind(get_layout<circle_topology<> >(), placeholders::_1,
                       placeholders::_2, placeholders::_3, make_pair(a, r),
                       scale, grid, make_pair(ti, tf), max_iter),
             vertex_floating_vector_properties(), edge_props_t()) (pos, weight);
}

#include <boost/python.hpp>

void export_fruchterman_reingold()
{
    python::def("fruchterman_reingold_layout", &fruchterman_reingold_layout);
}
