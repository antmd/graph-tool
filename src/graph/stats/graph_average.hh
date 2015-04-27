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


template <class Val1, class Val2>
void operator+=(std::vector<Val1>& a, const std::vector<Val2>& b)
{
    a.resize(std::max(a.size(), b.size()));
    for (size_t i = 0; i < std::min(a.size(), b.size()); ++i)
        a[i] += b[i];
}

template <class Val1, class Val2>
std::vector<Val1> operator*(const std::vector<Val1>& a, const std::vector<Val2>& b)
{
    std::vector<Val1> c(std::max(a.size(), b.size()));
    for (size_t i = 0; i < std::min(a.size(), b.size()); ++i)
        c[i] = a[i] * b[i];
    return c;
}

template <class Val1, class Val2>
std::vector<Val1> operator-(const std::vector<Val1>& a, const std::vector<Val2>& b)
{
    std::vector<Val1> c(std::max(a.size(), b.size()));
    for (size_t i = 0; i < std::min(a.size(), b.size()); ++i)
        c[i] = a[i] - b[i];
    for (size_t i = a.size(); i < std::max(a.size(), b.size()); ++i)
        c[i] = -b[i];
    return c;
}

struct get_avg_type
{
    template <class Type>
    struct apply
    {
        typedef typename mpl::if_<typename std::is_same<Type, python::object>::type,
                                  python::object, long double>::type
            type;
    };

    template <class Type>
    struct apply<vector<Type>>
    {
        typedef vector<long double> type;
    };
};

class VertexAverageTraverse
{
public:
    template <class Graph, class DegreeSelector, class ValueType>
    void operator()(Graph& g, typename graph_traits<Graph>::vertex_descriptor v,
                    DegreeSelector& deg, ValueType& a, ValueType& aa,
                    size_t& count)
    {
        const auto& x = deg(v, g);
        a += x;
        aa += x * x;
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
        for (auto e : out_edges_range(v, g))
        {
            const auto& x = eprop[e];
            a += x;
            aa += x * x;
            count++;
        }
    }
};

// explicit special initialization to get around python::object
template <class Val>
void init_avg(Val& v)
{
    v = Val(0.);
};

template <class Val>
void init_avg(std::vector<Val>&)
{
};

// generalized functor to obtain average of different types of "degrees"
template <class AverageTraverse>
struct get_average
{
    get_average(boost::python::object& a, boost::python::object& dev,
                size_t& count)
        : _a(a), _dev(dev), _count(count) {}

    template <class Graph, class DegreeSelector>
    void operator()(Graph& g, DegreeSelector deg) const
    {
        typedef typename DegreeSelector::value_type val_t;
        dispatch(g, deg, typename std::is_pod<val_t>::type());
    }

    template <class Graph, class DegreeSelector>
    void dispatch(Graph& g, DegreeSelector deg, std::true_type) const
    {
        typedef typename get_avg_type::apply<typename DegreeSelector::value_type>::type val_t;
        val_t a = 0., aa = 0.;
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

        _a = boost::python::object(a);
        _dev = boost::python::object(aa);
        _count = count;
    }

    template <class Graph, class DegreeSelector>
    void dispatch(Graph& g, DegreeSelector deg, std::false_type) const
    {
        typedef typename get_avg_type::apply<typename DegreeSelector::value_type>::type val_t;
        val_t a, aa;
        init_avg(a);
        init_avg(aa);
        size_t count = 0;

        AverageTraverse traverse;
        for (auto v : vertices_range(g))
            traverse(g, v, deg, a, aa, count);

        _a = boost::python::object(a);
        _dev = boost::python::object(aa);
        _count = count;
    }

    boost::python::object& _a;
    boost::python::object& _dev;
    size_t& _count;
};

} // graph_tool namespace

#endif // GRAPH_AVERAGE_HH
