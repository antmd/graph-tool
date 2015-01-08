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

//  As a special exception, you have permission to link this program
//  with the CGAL library and distribute executables, as long as you
//  follow the requirements of the GNU GPL in regard to all of the
//  software in the executable aside from CGAL.

#include "graph.hh"
#include "graph_util.hh"
#include "graph_filtering.hh"

#if (GCC_VERSION < 40400 || defined __clang__)
#   define CGAL_CFG_NO_TR1_ARRAY
#   define CGAL_CFG_NO_TR1_TUPLE
#endif
#ifdef __clang__
#   define CGAL_CFG_ARRAY_MEMBER_INITIALIZATION_BUG
#endif


#include <CGAL/version.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Triangulation_3<Kernel> SimpleTriangulation;
typedef CGAL::Delaunay_triangulation_3<Kernel> DelaunayTriangulation;

namespace CGAL
{
bool operator==(const SimpleTriangulation::Vertex& a,
                const SimpleTriangulation::Vertex& b)
{
    return a.point() == b.point();
}
}

// periodic triangulation is only available in more recent versions of CGAL
#if (CGAL_VERSION_NR >= 1030500000)
#include <CGAL/Periodic_3_triangulation_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
typedef CGAL::Periodic_3_triangulation_traits_3<Kernel> GT;
typedef CGAL::Periodic_3_Delaunay_triangulation_3<GT>
    PeriodicDelaunayTriangulation;
namespace CGAL
{
bool operator==(const PeriodicDelaunayTriangulation::Vertex& a,
                const PeriodicDelaunayTriangulation::Vertex& b)
{
    return a.point() == b.point();
}
}
#endif

#include "graph_triangulation.hh"
#include "numpy_bind.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

void triangulation(GraphInterface& gi, boost::python::object points,
                   boost::any pos, string type, bool periodic)
{
    UndirectedAdaptor<GraphInterface::multigraph_t> g(gi.GetGraph());
    multi_array_ref<double,2> points_array = get_array<double,2>(points);
    typedef property_map_type::apply
        <vector<double>, GraphInterface::vertex_index_map_t>::type pos_type_t;
    pos_type_t pos_map = any_cast<pos_type_t>(pos);

    if (type == "simple")
    {
        get_triangulation<SimpleTriangulation>()(g, points_array, pos_map);
    }
    else if (type == "delaunay")
    {
        if (!periodic)
        {
            get_triangulation<DelaunayTriangulation>()(g, points_array,
                                                       pos_map);
        }
        else
        {
#if (CGAL_VERSION_NR >= 1030500000)
            get_triangulation<PeriodicDelaunayTriangulation>()(g, points_array,
                                                               pos_map);
#else
            throw ValueException("Periodic Delaunay triangulation is only "
                                 "available with versions of CGAL newer than "
                                 "3.5.0.");
#endif
        }
    }
}
