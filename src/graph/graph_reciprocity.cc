// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006  Tiago de Paula Peixoto <tiago@forked.de>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include "graph.hh"
#include "histogram.hh"
#include "graph_filtering.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

struct get_reciprocity
{
    template <class Graph>
    void operator()(Graph& g, double& reciprocity) const
    {        
        size_t L = 0;
        double Lbd = 0.0;

        typename graph_traits<Graph>::edge_iterator e,e_end;
        for (tie(e,e_end) = edges(g); e != e_end; ++e)
        {
            typename graph_traits<Graph>::vertex_descriptor s,t;
            s = source(*e,g);
            t = target(*e,g);

            size_t o = 0;
            typename graph_traits<Graph>::adjacency_iterator a, a_end;
            for (tie(a,a_end) = adjacent_vertices(s, g); a != a_end; ++a)
                if (*a == t)
                    o++;

            size_t i = 0;
            for (tie(a, a_end) = adjacent_vertices(t, g); a != a_end; ++a)
                if (*a == s)
                    i++;

            Lbd += min(i/double(o),1.0);
            L++;
        }

        size_t N = HardNumVertices()(g);
        double a = L/double(N*(N-1));

        reciprocity = (Lbd/L - a)/(1-a);
    }
};


double GraphInterface::GetReciprocity() const
{
    double reciprocity;
    check_filter(*this, bind<void>(get_reciprocity(), _1, var(reciprocity)), reverse_check(), directed_check()); 
    return reciprocity;
}
