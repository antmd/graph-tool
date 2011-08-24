// Copyright (C) 2008  Tiago de Paula Peixoto <tiago@skewed.de>
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
#include "graph_properties.hh"
#include "graph_selectors.hh"

#include "graph_similarity.hh"

#include <boost/python.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

struct get_pointers
{
    template <class List>
    struct apply
    {
        typedef typename mpl::transform<List,
                                        mpl::quote1<add_pointer> >::type type;
    };
};


size_t similarity(GraphInterface& gi1, GraphInterface& gi2, boost::any label1,
                  boost::any label2)
{
    size_t s = 0;
    run_action<>()
        (gi1, bind<void>(get_similarity(), _1, _2, _3, label2, ref(s)),
         get_pointers::apply<graph_tool::detail::all_graph_views>::type(),
         vertex_scalar_properties())(gi2.GetGraphView(), label1);
    return s;
}

void export_similarity()
{
    python::def("similarity", &similarity);
};
