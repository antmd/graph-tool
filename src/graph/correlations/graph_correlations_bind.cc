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

void export_assortativity();
void export_vertex_correlations();
void export_combined_vertex_correlations();
void export_avg_correlations();
void export_avg_combined_correlations();

BOOST_PYTHON_MODULE(libgraph_tool_correlations)
{
    export_assortativity();
    export_vertex_correlations();
    export_combined_vertex_correlations();
    export_avg_correlations();
    export_avg_combined_correlations();
}
