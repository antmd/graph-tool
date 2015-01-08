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

#include <boost/graph/topological_sort.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

struct get_topological_sort
{
    template <class Graph>
    void operator()(const Graph& g, vector<int32_t>& sort) const
    {
        sort.clear();
        topological_sort(g, std::back_inserter(sort));
    }
};

bool topological_sort(GraphInterface& gi, vector<int32_t>& sort)
{
    try
    {
        run_action<>()
            (gi, std::bind(get_topological_sort(), placeholders::_1,
                           std::ref(sort)))();
        return true;
    }
    catch (not_a_dag& e)
    {
        return false;
    }
}
