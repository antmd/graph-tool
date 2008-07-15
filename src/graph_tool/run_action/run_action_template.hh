// Copyright (C) 2008  Tiago de Paula Peixoto <tiago@forked.de>
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

#include <map>
#include <set>
#include <list>
#include <tr1/unordered_set>
#include <tr1/unordered_map>
#include <boost/tuple/tuple.hpp>
#include <boost/type_traits.hpp>
#include <boost/python.hpp>
#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_properties.hh"
#include "histogram.hh"
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

using namespace boost;
using namespace boost::tuples;
using namespace std;
using namespace graph_tool;

namespace graph_tool
{
// arbitrary code execution
template <class Action>
void RunAction(GraphInterface &g, const Action& a)
{
    run_action<>()(g, lambda::bind<void>(a, lambda::_1, g._vertex_index,
                                         g._edge_index))();
}
}

// metafunction to get the correct property map
template <class IndexMap>
struct prop_bind_t
{
    template <class Value>
    struct as
    {
        typedef typename mpl::if_<is_same<Value,bool>,
                                  uint8_t, Value>::type val_t;
        typedef vector_property_map<val_t,IndexMap> type;
    };
};

// the action function object
template <class Args>
struct action_${code_hash}
{
    action_${code_hash}(const Args& args, py::object& return_val)
                           : _args(args), _return_val(return_val) {}

    template <class Graph, class VertexIndex, class EdgeIndex>
    void operator()(Graph* __gp, VertexIndex vertex_index,
                    EdgeIndex edge_index) const
    {
        Graph& g = *__gp;
        // convenience typedefs

        // descriptors and iterators
        typedef typename graph_traits<Graph>::vertex_descriptor
            vertex_t;
        typedef typename graph_traits<Graph>::vertex_iterator
            vertex_iter_t;
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;
        typedef typename graph_traits<Graph>::edge_iterator edge_iter_t;
        typedef typename graph_traits<Graph>::out_edge_iterator
            out_edge_iter_t;
        typedef typename in_edge_iteratorS<Graph>::type in_edge_iter_t;

        // property maps of all types
        typedef prop_bind_t<VertexIndex> vertex_prop_t;
        typedef prop_bind_t<EdgeIndex> edge_prop_t;
        typedef prop_bind_t<ConstantPropertyMap<size_t,graph_property_tag> >
            graph_prop_t;

        // for all types this should include a
        //  typedef typename vertex_prop_t::template as<type>::type
        //      vprop_type_t;
        // both for vertex, edges and graph properties
        ${property_map_types}

        // the arguments will be expanded below
        ${arg_expansion}

        bool __exception_thrown = false;
        string __exception_error;
        try
        {
            // where the actual code is included
            ${code}
        }
        catch (const GraphException& e)
        {
            __exception_error = e.what();
            __exception_thrown = true;
        }
        catch (const bad_any_cast& e)
        {
            __exception_error = e.what();
            __exception_error += " (wrong property map type?)";
            __exception_thrown = true;
        }
        catch (const std::exception& e)
        {
            __exception_error = "unknown exception thrown: ";
            __exception_error += e.what();
            __exception_thrown = true;
        }

        python::dict return_vals;
        return_vals["__exception_error"] = __exception_error;
        return_vals["__exception_thrown"] = __exception_thrown;

        // updated values will be inseted in return_vals below
        ${return_vals}

        _return_val = py::object(return_vals.ptr());
    }

    const Args& _args;
    py::object& _return_val;
};

// convenience function
template <class Args>
action_${code_hash}<Args> make_action(const Args& args, py::object& return_val)
{
    return action_${code_hash}<Args>(args, return_val);
}

