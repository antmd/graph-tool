// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007-2010 Tiago de Paula Peixoto <tiago@forked.de>
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

#include <boost/python.hpp>
#include <boost/lambda/bind.hpp>

#include "graph.hh"
#include "graph_selectors.hh"
#include "graph_pagerank.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

size_t pagerank(GraphInterface& g, boost::any rank, double d, double epslon,
                size_t max_iter)
{
    if (!belongs<writable_vertex_scalar_properties>()(rank))
        throw ValueException("vertex property must be writable");

    size_t iter;
    run_action<>()
        (g, bind<void>(get_pagerank(),
                       _1, g.GetVertexIndex(),  _2, d,
                       epslon, max_iter, ref(iter)),
         writable_vertex_scalar_properties())(rank);
    return iter;
}


void export_pagerank()
{
    using namespace boost::python;
    def("get_pagerank", &pagerank);
}
