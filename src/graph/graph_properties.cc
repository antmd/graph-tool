// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2013 Tiago de Paula Peixoto <tiago@skewed.de>
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

#include "graph_python_interface.hh"
#include "graph.hh"
#include "graph_properties.hh"
#include "graph_filtering.hh"
#include "graph_selectors.hh"
#include "graph_util.hh"

#if (GCC_VERSION >= 40400)
#   include <tr1/unordered_set>
#else
#   include <boost/tr1/unordered_set.hpp>
#endif

#include <boost/mpl/for_each.hpp>

#include <boost/python/extract.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

namespace graph_tool
{

// global property types' names
const char* type_names[] =
    {"bool", "int16_t", "int32_t", "int64_t", "double", "long double",
     "string", "vector<bool>", "vector<int16_t>", "vector<int32_t>",
     "vector<int64_t>", "vector<double>", "vector<long double>",
     "vector<string>", "python::object"};


struct shift_vertex_property
{
    template <class PropertyMap>
    void operator()(PropertyMap, const GraphInterface::multigraph_t& g,
                    boost::any map, size_t vi, bool& found) const
    {
        try
        {
            PropertyMap pmap = any_cast<PropertyMap>(map);
            for (size_t i = vi; i < num_vertices(g)-1; ++i)
                pmap[vertex(i,g)] = pmap[vertex(i+1,g)];
            found = true;
        }
        catch (bad_any_cast&) {}
    }
};

// this function will shift all the properties when a vertex is to be deleted
void GraphInterface::ShiftVertexProperty(boost::any prop, size_t index) const
{
    bool found = false;
    mpl::for_each<writable_vertex_properties>
        (bind<void>(shift_vertex_property(), _1, ref(_state->_mg),
                    prop, index, ref(found)));
    if (!found)
        throw GraphException("invalid writable property map");
}

struct reindex_vertex_property
{
    template <class PropertyMap, class IndexMap>
    void operator()(PropertyMap, const GraphInterface::multigraph_t& g,
                    boost::any map, IndexMap old_index, bool& found) const
    {
        try
        {
            PropertyMap pmap = any_cast<PropertyMap>(map);
            for (size_t i = 0; i < num_vertices(g); ++i)
            {
                GraphInterface::vertex_t v = vertex(i, g);
                if (old_index[v] != int(i))
                    pmap[v] = pmap[vertex(old_index[v], g)];
            }
            found = true;
        }
        catch (bad_any_cast&) {}
    }
};


void GraphInterface::ReIndexVertexProperty(boost::any map,
                                           boost::any aold_index) const
{
    typedef property_map_type::apply<int32_t,
                                     GraphInterface::vertex_index_map_t>::type
        index_prop_t;
    index_prop_t old_index = any_cast<index_prop_t>(aold_index);

    bool found = false;
    mpl::for_each<writable_vertex_properties>
        (bind<void>(reindex_vertex_property(), _1, ref(_state->_mg),
                    map, old_index, ref(found)));
    if (!found)
        throw GraphException("invalid writable property map");

}

} // graph_tool namespace


struct do_infect_vertex_property
{
    template <class Graph, class IndexMap, class PropertyMap>
    void operator()(Graph& g, IndexMap index, PropertyMap prop,
                    python::object oval) const
    {
        typedef typename property_traits<PropertyMap>::value_type val_t;
        bool all = false;

        tr1::unordered_set<val_t, boost::hash<val_t> > vals;
        if (oval == python::object())
        {
            all = true;
        }
        else
        {
            for (int i = 0; i < len(oval); ++i)
            {
                val_t val = python::extract<val_t>(oval[i]);
                vals.insert(val);
            }
        }

        unchecked_vector_property_map<uint8_t, IndexMap>
            marked(index, num_vertices(g));

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            bool skip;
            {
                #pragma omp critical
                skip = marked[v];
            }
            if (skip)
                continue;
            if (!all && vals.find(prop[v]) == vals.end())
                continue;
            typename graph_traits<Graph>::adjacency_iterator a, a_end;
            for (tie(a, a_end) = adjacent_vertices(v, g); a != a_end; ++a)
            {
                if (prop[*a] == prop[v])
                    continue;
                {
                    #pragma omp critical
                    marked[*a] = true;
                }
                prop[*a] = prop[v];
            }
        }
    }
};

void infect_vertex_property(GraphInterface& gi, boost::any prop,
                            python::object val)
{
    run_action<>()(gi, bind<void>(do_infect_vertex_property(), _1,
                                  gi.GetVertexIndex(), _2, val),
                   writable_vertex_properties())(prop);
}

template <class Value>
vector<Value> operator-(const vector<Value>& a, const vector<Value>& b)
{
    vector<Value> c(a);
    c.resize(max(a.size(), b.size()), Value(0));
    for (size_t i = 0; i < b.size(); ++i)
        c[i] = a[i] - b[i];
    return c;
}

struct do_edge_difference
{
    template <class Graph, class EdgeIndexMap, class VertexPropertyMap>
    void operator()(Graph& g, EdgeIndexMap edge_index, VertexPropertyMap prop,
                    boost::any eprop) const
    {
        typedef typename property_traits<VertexPropertyMap>::value_type vval_t;
        typedef typename mpl::if_<is_same<vval_t, size_t>, int32_t, vval_t>::type
            val_t;
        typedef typename property_map_type::apply<val_t, EdgeIndexMap>::type
            eprop_t;
        eprop_t ediff = any_cast<eprop_t>(eprop);
        ediff.reserve(num_edges(g));

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            typename graph_traits<Graph>::out_edge_iterator e, e_end;
            for (tie(e, e_end) = out_edges(v, g); e != e_end; ++e)
                ediff[*e] = prop[target(*e, g)] - prop[source(*e, g)];
        }
    }
};

void edge_difference(GraphInterface& gi, boost::any prop,
                     boost::any eprop)
{
    typedef mpl::insert_range<vertex_scalar_properties,
                              mpl::end<vertex_scalar_properties>::type,
                              vertex_scalar_vector_properties>::type vprops_t;
    run_action<>()(gi, bind<void>(do_edge_difference(), _1,
                                  gi.GetEdgeIndex(), _2, eprop),
                   vprops_t())(prop);
}
