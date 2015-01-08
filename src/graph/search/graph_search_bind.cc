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

#include <boost/python.hpp>

using namespace boost;

void export_bfs();
void export_dfs();
void export_dijkstra();
void export_bellman_ford();
void export_astar();
void export_astar_implicit();

BOOST_PYTHON_MODULE(libgraph_tool_search)
{
    export_bfs();
    export_dfs();
    export_dijkstra();
    export_bellman_ford();
    export_astar();
    export_astar_implicit();
}
