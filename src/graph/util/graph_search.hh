// Copyright (C) 2008  Tiago de Paula Peixoto <tiago@skewed.de>
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

#ifndef GRAPH_SEARCH_HH
#define GRAPH_SEARCH_HH

#include "graph_python_interface.hh"
#include "graph_util.hh"

#if (GCC_VERSION >= 40400)
#   include <tr1/unordered_set>
#else
#   include <boost/tr1/unordered_set.hpp>
#endif

#ifdef USING_OPENMP
#include <omp.h>
#include <boost/type_traits.hpp>
#endif


namespace graph_tool
{
using namespace std;
using namespace boost;

// less than or equal for vectors. Compare sizes, or distance from origin.
template <class ValueType>
bool operator<=(const vector<ValueType>& v1, const vector<ValueType>& v2)
{
    if (v1.size() < v2.size())
        return true;
    if (v1.size() > v2.size())
        return false;
    long double d1 = 0, d2 = 0;
    for (size_t i = 0; i < v1.size(); ++i)
    {
        d1 += v1[i]*v1[i];
        d2 += v2[i]*v2[i];
    }
    return (d1 <= d2);
}

// sort strings in alphabetical (ASCII) order
bool operator<=(const string s1, const string s2)
{
    for (size_t i = 0; i < s1.size(); ++i)
    {
        if (s1[i] != s1[i])
            return (s1[i] < s1[i]);
    }
    return (s1.size() == s2.size());
}

// vector of strings. compare element by element
bool operator<=(const vector<string>& v1, const vector<string>& v2)
{
    if (v1.size() < v2.size())
        return true;
    if (v1.size() > v2.size())
        return false;
    for (size_t i = 0; i < v1.size(); ++i)
    {
        if (v1[i] > v2[i])
            return false;
    }
    return true;
}

// find vertices which match a certain (inclusive) property range
struct find_vertices
{
    template <class Graph, class DegreeSelector>
    void operator()(Graph& g, const python::object& pg, DegreeSelector deg,
                    python::tuple& prange, python::list& ret) const
    {
        typedef typename DegreeSelector::value_type value_type;
        pair<value_type,value_type> range;
        range.first = python::extract<value_type>(prange[0]);
        range.second = python::extract<value_type>(prange[1]);

        #ifdef USING_OPENMP
        size_t nt = omp_get_num_threads();
        if (is_convertible<value_type,python::object>::value)
            nt = 1; // python is not thread-safe
        #endif

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) schedule(dynamic) \
            num_threads(nt)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            value_type val = deg(v, g);
            if (val >= range.first && val <= range.second)
            {
                PythonVertex pv(pg, v);
                {
                    #pragma omp critical
                    ret.append(pv);
                }
            }
        }
    }
};

// find edges which match a certain (inclusive) property range
struct find_edges
{
    template <class Graph, class EdgeIndex, class EdgeProperty>
    void operator()(Graph& g, const python::object& pg, EdgeIndex eindex,
                    EdgeProperty prop, python::tuple& prange, python::list& ret)
        const
    {
        typedef typename property_traits<EdgeProperty>::value_type value_type;
        pair<value_type,value_type> range;
        range.first = python::extract<value_type>(prange[0]);
        range.second = python::extract<value_type>(prange[1]);

        tr1::unordered_set<size_t> edge_set;

        #ifdef USING_OPENMP
        size_t nt = omp_get_num_threads();
        if (is_convertible<value_type,python::object>::value)
            nt = 1; // python is not thread-safe
        #endif

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) schedule(dynamic) \
            num_threads(nt)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            typename graph_traits<Graph>::out_edge_iterator e, e_end;
            for (tie(e, e_end) = out_edges(v, g); e != e_end; ++e)
            {
                if (!is_directed::apply<Graph>::type::value)
                {
                    if (edge_set.find(eindex[*e]) == edge_set.end())
                        edge_set.insert(eindex[*e]);
                    else
                        continue;
                }

                value_type val = get(prop, *e);
                if (val >= range.first && val <= range.second)
                {
                    PythonEdge<Graph> pe(pg, *e);
                    {
                        #pragma omp critical
                        ret.append(pe);
                    }
                }
            }
        }
    }
};


} // graph_tool namespace

#endif // GRAPH_SEARCH_HH
