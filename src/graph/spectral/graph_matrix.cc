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
#include "graph.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

void adjacency(GraphInterface& g, boost::any index, boost::any weight,
               python::object odata, python::object oi,
               python::object oj);


void laplacian(GraphInterface& g, boost::any index, boost::any weight,
               string sdeg,
               python::object odata, python::object oi,
               python::object oj);


void norm_laplacian(GraphInterface& g, boost::any index, boost::any weight,
                    string sdeg,
                    python::object odata, python::object oi,
                    python::object oj);

void incidence(GraphInterface& g, boost::any vindex, boost::any eindex,
               python::object odata, python::object oi,
               python::object oj);

void transition(GraphInterface& g, boost::any index, boost::any weight,
                python::object odata, python::object oi,
                python::object oj);

BOOST_PYTHON_MODULE(libgraph_tool_spectral)
{
    using namespace boost::python;
    def("adjacency", &adjacency);
    def("laplacian", &laplacian);
    def("norm_laplacian", &norm_laplacian);
    def("incidence", &incidence);
    def("transition", &transition);
}
