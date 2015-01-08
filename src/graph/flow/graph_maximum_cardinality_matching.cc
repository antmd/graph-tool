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
#include "graph_selectors.hh"
#include "graph_properties.hh"
#include "graph.hh"

#include "graph_augment.hh"

#include <boost/graph/max_cardinality_matching.hpp>

#include <boost/bind.hpp>

using namespace graph_tool;
using namespace boost;

struct get_max_cardinality_matching
{
    template <class Graph, class VertexIndex, class MatchMap>
    void operator()(Graph& g, VertexIndex vertex_index, MatchMap match,
                    bool &check) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        unchecked_vector_property_map<vertex_t,VertexIndex>
            mate(vertex_index, num_vertices(g));

        typename graph_traits<Graph>::edge_iterator e, e_end;
        for (tie(e, e_end) = edges(g); e != e_end; ++e)
            match[*e] = false;

        check =
            checked_edmonds_maximum_cardinality_matching(g, mate, vertex_index);

        for (tie(e, e_end) = edges(g); e != e_end; ++e)
        {
            if (mate[source(*e,g)] != graph_traits<Graph>::null_vertex() &&
                mate[source(*e,g)] == target(*e,g))
            {
                // search for an already matched parallel edge
                bool matched = false;
                typename graph_traits<Graph>::out_edge_iterator oe, oe_end;
                for (tie(oe, oe_end) = out_edges(source(*e, g),g); oe != oe_end;
                     ++oe)
                {
                    if (match[*oe])
                    {
                        matched = true;
                        break;
                    }
                }
                if (!matched)
                    match[*e] = true;
            }
        }
    }
};


bool max_cardinality_matching(GraphInterface& gi, boost::any match)
{
    bool check;
    run_action<graph_tool::detail::never_directed>()
        (gi, std::bind(get_max_cardinality_matching(),
                        placeholders::_1, gi.GetVertexIndex(),
                       placeholders::_2, std::ref(check)),
         writable_edge_scalar_properties()) (match);
    return check;
}
