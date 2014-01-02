// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2014 Tiago de Paula Peixoto <tiago@skewed.de>
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

#ifndef GRAPH_INCIDENCE_HH
#define GRAPH_INCIDENCE_HH

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_util.hh"

namespace graph_tool
{
using namespace boost;

struct get_incidence
{
    template <class Graph, class VIndex, class EIndex>
    void operator()(Graph& g, VIndex vindex, EIndex eindex,
                    multi_array_ref<double,1>& data,
                    multi_array_ref<int32_t,1>& i,
                    multi_array_ref<int32_t,1>& j) const
    {
        int pos = 0;
        typename graph_traits<Graph>::vertex_iterator v, v_end;
        for(tie(v, v_end) = vertices(g); v != v_end; ++v)
        {

            typename graph_traits<Graph>::out_edge_iterator e, e_end;
            for(tie(e, e_end) = out_edges(*v, g); e != e_end; ++e)
            {
                if (is_directed::apply<Graph>::type::value)
                    data[pos] = -1;
                else
                    data[pos] = 1;
                i[pos] = get(vindex, *v);
                j[pos] = get(eindex, *e);
                ++pos;
            }

            typename in_edge_iteratorS<Graph>::type ie, ie_end;
            for(tie(ie, ie_end) = in_edge_iteratorS<Graph>::get_edges(*v, g);
                ie != ie_end; ++ie)
            {
                data[pos] = 1;
                i[pos] = get(vindex, *v);
                j[pos] = get(eindex, *ie);
                ++pos;
            }
        }
    }
};

} // namespace graph_tool

#endif // GRAPH_INCIDENCE_HH
