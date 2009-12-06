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

//  As a special exception, you have permission to link this program
//  with the CGAL library and distribute executables, as long as you
//  follow the requirements of the GNU GPL in regard to all of the
//  software in the executable aside from CGAL.

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Triangulation_3<Kernel> SimpleTriangulation;
typedef CGAL::Delaunay_triangulation_3<Kernel> DelaunayTriangulation;

namespace std
{
bool operator==(const SimpleTriangulation::Vertex& a,
                const SimpleTriangulation::Vertex& b)
{
    return a.point() == b.point();
}
}

#include "graph.hh"
#include "graph_util.hh"
#include "graph_filtering.hh"
#include "graph_triangulation.hh"
#include "numpy_bind.hh"


using namespace std;
using namespace boost;
using namespace graph_tool;

void triangulation(GraphInterface& gi, python::object points, boost::any pos,
                   string type)
{
    UndirectedAdaptor<GraphInterface::multigraph_t> g(gi.GetGraph());
    multi_array_ref<double,2> points_array = get_array<double,2>(points);
    typedef property_map_type::apply
        <vector<double>, GraphInterface::vertex_index_map_t>::type pos_type_t;
    pos_type_t pos_map = any_cast<pos_type_t>(pos);

    if (type == "simple")
        get_triangulation<SimpleTriangulation>()(g, points_array, pos_map);
    else if (type == "delaunay")
        get_triangulation<DelaunayTriangulation>()(g, points_array,
                                                   pos_map);
    gi.ReIndexEdges();
}
