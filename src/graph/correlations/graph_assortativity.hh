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

#ifndef GRAPH_ASSORTATIVITY_HH
#define GRAPH_ASSORTATIVITY_HH

#include <unordered_map>

#include "shared_map.hh"
#include "graph_util.hh"

namespace graph_tool
{
using namespace std;
using namespace boost;


// this will calculate the assortativity coefficient, based on the property
// pointed by 'deg'

struct get_assortativity_coefficient
{
    template <class Graph, class DegreeSelector>
    void operator()(const Graph& g, DegreeSelector deg, double& r,
                    double& r_err) const
    {
        typedef typename mpl::if_<typename is_directed::apply<Graph>::type,
                                  size_t, double>::type count_t;

        count_t c = (is_directed::apply<Graph>::type::value) ? count_t(1) : count_t(0.5);
        count_t n_edges = 0;
        count_t e_kk = 0;

        typedef typename DegreeSelector::value_type val_t;
        typedef unordered_map<val_t, count_t> map_t;
        map_t a, b;

        SharedMap<map_t> sa(a), sb(b);
        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) firstprivate(sa,sb)\
            schedule(runtime) if (N > 100) reduction(+:e_kk, n_edges)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            val_t k1 = deg(v, g);
            typename graph_traits<Graph>::out_edge_iterator e, e_end;
            for (tie(e,e_end) = out_edges(v, g); e != e_end; ++e)
            {
                val_t k2 = deg(target(*e, g), g);
                if (k1 == k2)
                    e_kk += c;
                sa[k1] += c;
                sb[k2] += c;
                n_edges += c;
            }
        }

        sa.Gather();
        sb.Gather();

        double t1 = double(e_kk) / n_edges, t2 = 0.0;

        for (typeof(a.begin()) iter = a.begin(); iter != a.end(); ++iter)
            if (b.find(iter->second) != b.end())
                t2 += double(iter->second * b[iter->first]);
        t2 /= n_edges*n_edges;

        r = (t1 - t2)/(1.0 - t2);

        // "jackknife" variance
        double err = 0.0;
        #pragma omp parallel for default(shared) private(i) schedule(runtime) if (N > 100)\
            reduction(+:err)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            val_t k1 = deg(v, g);
            typename graph_traits<Graph>::out_edge_iterator e, e_end;
            for (tie(e,e_end) = out_edges(v, g); e != e_end; ++e)
            {
                val_t k2 = deg(target(*e,g), g);
                double tl2 = (t2*(n_edges*n_edges) - b[k1] - a[k2])/
                    ((n_edges-1.)*(n_edges-1.));
                double tl1 = t1*n_edges;
                if (k1 == k2)
                    tl1 -= 1;
                tl1 /= n_edges - 1;
                double rl = (tl1 - tl2)/(1.0 - tl2);
                err += (r-rl)*(r-rl)*c;
            }
        }
        r_err = sqrt(err);
    }
};

// this will calculate the _scalar_ assortativity coefficient, based on the
// scalar property pointed by 'deg'

struct get_scalar_assortativity_coefficient
{
    template <class Graph, class DegreeSelector>
    void operator()(const Graph& g, DegreeSelector deg, double& r,
                    double& r_err) const
    {
        typedef typename mpl::if_<typename is_directed::apply<Graph>::type,
                                  size_t, double>::type count_t;

        count_t c = (is_directed::apply<Graph>::type::value) ? count_t(1) : count_t(0.5);
        count_t n_edges = 0;
        double e_xy = 0.0;
        double a = 0.0, b = 0.0, da = 0.0, db = 0.0;
        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) \
            schedule(runtime) if (N > 100) reduction(+:e_xy,n_edges,a,b,da,db)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            double k1 = double(deg(v, g));
            typename graph_traits<Graph>::out_edge_iterator e, e_end;
            for (tie(e,e_end) = out_edges(v, g); e != e_end; ++e)
            {
                double k2 = double(deg(target(*e,g),g));
                a += k1*c;
                da += k1*k1*c;
                b += k2*c;
                db += k2*k2*c;
                e_xy += k1*k2*c;
                n_edges += c;
            }
        }

        double t1 = e_xy/n_edges;
        a /= n_edges;
        b /= n_edges;
        double stda = sqrt(da/n_edges - a*a);
        double stdb = sqrt(db/n_edges - b*b);

        if (stda*stdb > 0)
            r = (t1 - a*b)/(stda*stdb);
        else
            r = (t1 - a*b);

        // "jackknife" variance
        r_err = 0.0;

        double err = 0.0;
        #pragma omp parallel for default(shared) private(i) schedule(runtime) if (N > 100)\
            reduction(+:err)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            double k1 = double(deg(v, g));
            double al = (a*n_edges - k1)/(n_edges-1);
            double dal = sqrt((da - k1*k1)/(n_edges-1) - al*al);

            typename graph_traits<Graph>::out_edge_iterator e, e_end;
            for (tie(e,e_end) = out_edges(v, g); e != e_end; ++e)
            {
                double k2 = double(deg(target(*e, g), g));
                double bl = (b*n_edges - k2)/(n_edges-1);
                double dbl = sqrt((db - k2*k2)/(n_edges-1) - bl*bl);
                double t1l = (e_xy - k1*k2)/(n_edges-1);
                double rl;
                if (dal*dbl > 0)
                    rl = (t1l - al*bl)/(dal*dbl);
                else
                    rl = (t1l - al*bl);
                err += (r-rl)*(r-rl)*c;
            }
        }
        r_err = sqrt(err);
    }
};

} // graph_tool namespace

#endif //GRAPH_ASSORTATIVITY_HH
