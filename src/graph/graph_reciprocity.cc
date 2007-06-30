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

        int i, NV = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) reduction(+:L,Lbd) schedule(dynamic)
        for (i = 0; i < NV; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            typename graph_traits<Graph>::out_edge_iterator e, e_begin, e_end;
            tie(e_begin,e_end) = out_edges(v,g);
            for(e = e_begin; e != e_end; ++e)
            {
                typename graph_traits<Graph>::vertex_descriptor s,t;
                s = v;
                t = target(*e,g);

                size_t o = 0;
                typename graph_traits<Graph>::adjacency_iterator a, a_end;
                for (tie(a,a_end) = adjacent_vertices(s, g); a != a_end; ++a)
                    if (*a == t)
                        o++;

                size_t j = 0;
                for (tie(a, a_end) = adjacent_vertices(t, g); a != a_end; ++a)
                    if (*a == s)
                        j++;

                Lbd += min(j/double(o),1.0);
                L++;
            }
        }

        if(is_convertible<typename graph_traits<Graph>::directed_category, undirected_tag>::value)
        {
            L /= 2;
            Lbd /= 2;
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
