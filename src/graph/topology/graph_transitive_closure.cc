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
#include <boost/graph/transitive_closure.hpp>

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_properties.hh"



using namespace graph_tool;
using namespace boost;

struct get_transitive_closure
{
    template <class Graph,  class TCGraph>
    void operator()(Graph& g, TCGraph& tcg) const
    {
        boost::transitive_closure(g, tcg);
    }
};


void transitive_closure(GraphInterface& gi, GraphInterface& tcgi)
{
    run_action<graph_tool::detail::always_directed>()
        (gi, std::bind(get_transitive_closure(), placeholders::_1,
                       std::ref(tcgi.GetGraph())))();
}
