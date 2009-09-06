// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007  Tiago de Paula Peixoto <tiago@forked.de>
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

#include "graph_union.hh"

#include <boost/bind.hpp>

#include <boost/python/extract.hpp>

using namespace graph_tool;
using namespace boost;

typedef property_map_type::apply<GraphInterface::vertex_t,
                                 GraphInterface::vertex_index_map_t>::type
    vprop_t;

typedef property_map_type::apply<GraphInterface::edge_t,
                                 GraphInterface::edge_index_map_t>::type
    eprop_t;

struct get_pointers
{
    template <class List>
    struct apply
    {
        typedef typename mpl::transform<List,
                                        mpl::quote1<add_pointer> >::type type;
    };
};

python::tuple graph_union(GraphInterface& ugi, GraphInterface& gi)
{
    vprop_t vprop(gi.GetVertexIndex());
    eprop_t eprop(gi.GetEdgeIndex());
    run_action<graph_tool::detail::always_directed,mpl::true_>()
        (ugi, bind<void>(graph_tool::graph_union(),
                         _1, _2, vprop, eprop),
         get_pointers::apply<graph_tool::detail::always_directed>::type())
        (gi.GetGraphView());
    return python::make_tuple(boost::any(vprop),boost::any(eprop));
}
