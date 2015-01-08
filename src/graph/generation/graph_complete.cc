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

#include "graph_complete.hh"

#include <boost/python.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

void complete(GraphInterface& gi, size_t N, bool directed,
              bool self_loops)
{
    get_complete()(gi.GetGraph(), N, directed, self_loops);
}

void circular(GraphInterface& gi, size_t N, size_t k, bool directed,
              bool self_loops)
{
    get_circular()(gi.GetGraph(), N, k, directed, self_loops);
}
