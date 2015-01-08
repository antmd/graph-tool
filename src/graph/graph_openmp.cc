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

#include <boost/python.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

// OpenMP control
#ifdef USING_OPENMP
#include <omp.h>
#endif

bool openmp_enabled()
{
#ifdef USING_OPENMP
    return true;
#else
    return false;
#endif
}

size_t openmp_get_num_threads()
{
#ifdef USING_OPENMP
    return omp_get_max_threads();
#else
    throw GraphException("OpenMP was not enabled during compilation");
#endif
}

void openmp_set_num_threads(int n)
{
#ifdef USING_OPENMP
    omp_set_num_threads(n);
#else
    throw GraphException("OpenMP was not enabled during compilation");
#endif
}

python::tuple openmp_get_schedule()
{
#ifdef USING_OPENMP
    omp_sched_t kind;
    int chunk;
    omp_get_schedule(&kind, &chunk);
    string skind;
    switch (kind)
    {
        case omp_sched_static:
            skind = "static";
            break;
        case omp_sched_dynamic:
            skind = "dynamic";
            break;
        case omp_sched_guided:
            skind = "guided";
            break;
        case omp_sched_auto:
            skind = "auto";
            break;
        default:
            throw GraphException("Unknown schedule type");
    }

    return python::make_tuple(skind, chunk);
#else
    throw GraphException("OpenMP was not enabled during compilation");
#endif
}

void openmp_set_schedule(string skind, int chunk)
{
#ifdef USING_OPENMP
    omp_sched_t kind;
    if (skind == "static")
        kind = omp_sched_static;
    else if (skind == "dynamic")
        kind = omp_sched_dynamic;
    else if (skind == "guided")
        kind = omp_sched_guided;
    else if (skind == "auto")
        kind = omp_sched_auto;
    else
        throw GraphException("Unknown schedule type: " + skind);
    omp_set_schedule(kind, chunk);
#else
    throw GraphException("OpenMP was not enabled during compilation");
#endif
}


void export_openmp()
{
    using namespace boost::python;

    def("openmp_enabled", &openmp_enabled);
    def("openmp_get_num_threads", &openmp_get_num_threads);
    def("openmp_set_num_threads", &openmp_set_num_threads);
    def("openmp_get_schedule", &openmp_get_schedule);
    def("openmp_set_schedule", &openmp_set_schedule);
};
