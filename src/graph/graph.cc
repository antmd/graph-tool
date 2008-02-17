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

#include <iostream>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include <unistd.h>    /* standard unix functions, like getpid()         */
#include <sys/types.h> /* various type definitions, like pid_t           */
#include <signal.h>    /* signal name macros, and the signal() prototype */

#include "graph.hh"
#include "histogram.hh"
#include "graph_filtering.hh"
#include "graph_properties.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;


// this is the constructor for the graph interface
GraphInterface::GraphInterface()
    :_mg(),
     _reversed(false),
     _directed(true),
     _vertex_index(get(vertex_index,_mg)),
     _edge_index(get(edge_index_t(),_mg)),
     _vertex_filter_map(_vertex_index),
     _vertex_filter_invert(false),
     _vertex_filter_active(false),
     _edge_filter_map(_edge_index),
     _edge_filter_invert(false),
     _edge_filter_active(false)
{

}

// the destructor
GraphInterface::~GraphInterface()
{
}

// this will get the number of vertices, either the "soft" O(1) way, or the hard
// O(V) way, which is necessary if the graph is filtered
size_t GraphInterface::GetNumberOfVertices() const
{
    size_t n = 0;
    if (IsVertexFilterActive())
        run_action<>()(*this, var(n)=bind<size_t>(HardNumVertices(),_1))();
    else
        run_action<>()(*this, var(n)=bind<size_t>(SoftNumVertices(),_1))();
    return n;
}

// this will get the number of edges, either the "soft" O(E) way, or the hard
// O(E) way, which is necessary if the graph is filtered. Both cases are of
// linear complexity, since num_edges() is O(E) in Boost's adjacency_list
size_t GraphInterface::GetNumberOfEdges() const
{
    size_t n = 0;
    if (IsEdgeFilterActive() || IsVertexFilterActive())
        run_action<>()(*this, var(n)=bind<size_t>(HardNumEdges(),_1))();
    else
        run_action<>()(*this, var(n)=bind<size_t>(SoftNumEdges(),_1))();
    return n;
}

// signal handling

void catch_sig_stop(int sig_num)
{
    std::cerr << "graph-tool: received signal ";
    switch (sig_num)
    {
    case SIGINT:
        std::cerr << "SIGINT (Interrupt).";
        break;
    case SIGTERM:
        std::cerr << "SIGTERM (Termination).";
        break;
    case SIGQUIT:
        std::cerr << "SIGQUIT (Terminal quit).";
        break;
    case SIGHUP:
        std::cerr << "SIGHUP (Hangup).";
        break;
    case SIGPWR:
        std::cerr << "SIGPWR (Power failure restart).";
        break;
    case SIGSEGV:
        std::cerr << "SIGSEGV (Segmentation fault). "
                  << "There's a bug somewhere in the program. Go fix it.";
        break;
    case SIGBUS:
        std::cerr << "SIGBUS (BUS error). The bus is broken, I guess...";
        break;
    case SIGFPE:
        std::cerr << "SIGFPE (Floating-point exception). "
                  << "INFs and NANs are wreaking havoc.";
        break;
    case SIGILL:
        std::cerr << "SIGILL (Illegal instruction). Did you compile it right?";
        break;
    case SIGXCPU:
        std::cerr << "SIGXCPU (CPU limit exceeded). Time's over.";
        break;
    case SIGXFSZ:
        std::cerr << "SIGXFSZ (File size limit exceeded). "
                  << "The fascist sysadmin doesn't let you write more.";
        break;
    }
    std::cerr << " Bailing out cowardly and calling abort()." << std::endl;
    abort();
}

// initialize signal handling

void GraphInterface::InitSignalHandling()
{
    signal(SIGINT, catch_sig_stop);
    signal(SIGTERM, catch_sig_stop);
    signal(SIGQUIT, catch_sig_stop);
    signal(SIGHUP, catch_sig_stop);
    signal(SIGPWR, catch_sig_stop);
    signal(SIGSEGV, catch_sig_stop);
    signal(SIGBUS, catch_sig_stop);
    signal(SIGFPE, catch_sig_stop);
    signal(SIGILL, catch_sig_stop);
    signal(SIGXCPU, catch_sig_stop);
    signal(SIGXFSZ, catch_sig_stop);
}
