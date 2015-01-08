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

#ifndef GRAPH_FDP_HH
#define GRAPH_FDP_HH

#include <limits>
#include <iostream>

#ifndef __clang__
#include <ext/numeric>
using __gnu_cxx::power;
#else
template <class Value>
Value power(Value value, int n)
{
    return pow(value, n);
}
#endif

namespace graph_tool
{
using namespace std;
using namespace boost;

template <class Pos, class Weight>
class QuadTree
{
public:
    QuadTree(const Pos& ll, const Pos& ur, int max_level)
        :_ll(ll), _ur(ur), _cm(2, 0), _count(0),
         _max_level(max_level)
    {
        _w = sqrt(power(_ur[0] - _ll[0], 2) +
                  power(_ur[1] - _ll[1], 2));
    }

    QuadTree& get_leaf(size_t i)
    {
        if (_max_level > 0 && _leafs.empty())
        {
            _leafs.reserve(4);
            for (size_t i = 0; i < 4; ++i)
            {
                Pos lll = _ll, lur = _ur;
                if (i % 2)
                    lll[0] += (_ur[0] - _ll[0]) / 2;
                else
                    lur[0] -= (_ur[0] - _ll[0]) / 2;
                if (i / 2)
                    lll[1] += (_ur[1] - _ll[1]) / 2;
                else
                    lur[1] -= (_ur[1] - _ll[1]) / 2;
                _leafs.push_back(QuadTree(lll, lur, _max_level - 1));
            }
        }

        return _leafs[i];
    }

    vector<std::tuple<Pos,Weight> >& get_dense_leafs()
    {
        return _dense_leafs;
    }

    size_t put_pos(Pos& p, Weight w)
    {
        _count += w;
        _cm[0] += p[0] * w;
        _cm[1] += p[1] * w;

        if (_max_level == 0)
        {
            if (_count == w)
                _dense_leafs.reserve(10);
            _dense_leafs.push_back(std::make_tuple(p, w));
            return 1;
        }
        else
        {
            int i = p[0] > (_ll[0] + (_ur[0] - _ll[0]) / 2);
            int j = p[1] > (_ll[1] + (_ur[1] - _ll[1]) / 2);
            size_t sc = (_max_level > 0 && _leafs.empty()) ? 4 : 0;
            return sc + get_leaf(i + 2 * j).put_pos(p, w);
        }
        return 0;
    }

    void get_cm(Pos& cm)
    {
        for (size_t i = 0; i < 2; ++i)
            cm[i] = _cm[i] / _count;
    }

    double get_w() const
    {
        return _w;
    }

    Weight get_count()
    {
        return _count;
    }

    int max_level()
    {
        return _max_level;
    }

private:
    Pos _ll, _ur;
    vector<QuadTree> _leafs;
    vector<std::tuple<Pos,Weight> > _dense_leafs;
    Pos _cm;
    Weight _count;
    int _max_level;
    double _w;
};

template <class Pos>
inline double dist(const Pos& p1, const Pos& p2)
{
    double r = 0;
    for (size_t i = 0; i < 2; ++i)
        r += power(double(p1[i] - p2[i]), 2);
    return sqrt(r);
}

template <class Pos>
inline double f_r(double C, double K, double p, const Pos& p1, const Pos& p2)
{
    double d = dist(p1, p2);
    if (d == 0)
        return 0;
    if (round(p) == p)
        return -C * power(K, int(1 + p)) / power(d, int(p));
    else
        return -C * pow(K, 1 + p) / pow(d, p);
}

template <class Pos>
inline double f_a(double K, const Pos& p1, const Pos& p2)
{
    return power(dist(p1, p2), 2) / K;
}

template <class Pos>
inline double get_diff(const Pos& p1, const Pos& p2, Pos& r)
{
    double abs = 0;
    for (size_t i = 0; i < 2; ++i)
    {
        r[i] = p1[i] - p2[i];
        abs += r[i] * r[i];
    }
    if (abs == 0)
        abs = 1;
    abs = sqrt(abs);
    for (size_t i = 0; i < 2; ++i)
        r[i] /= abs;
    return abs;
}

template <class Pos>
inline double norm(Pos& x)
{
    double abs = 0;
    for (size_t i = 0; i < 2; ++i)
        abs += power(x[i], 2);
    for (size_t i = 0; i < 2; ++i)
        x[i] /= sqrt(abs);
    return sqrt(abs);
}

struct get_sfdp_layout
{
    get_sfdp_layout(double C, double K, double p, double theta, double gamma,
                    double mu, double mu_p, double init_step,
                    double step_schedule, size_t max_level, double epsilon,
                    size_t max_iter, bool simple)
        : C(C), K(K), p(p), theta(theta), gamma(gamma), mu(mu), mu_p(mu_p),
          init_step(init_step), step_schedule(step_schedule),
          epsilon(epsilon), max_level(max_level), max_iter(max_iter),
          simple(simple) {}

    double C, K, p, theta, gamma, mu, mu_p, init_step, step_schedule, epsilon;
    size_t max_level, max_iter;
    bool simple;

    template <class Graph, class VertexIndex, class PosMap, class VertexWeightMap,
              class EdgeWeightMap, class PinMap, class GroupMap, class RNG>
    void operator()(Graph& g, VertexIndex vertex_index, PosMap pos,
                    VertexWeightMap vweight, EdgeWeightMap eweight, PinMap pin,
                    GroupMap group, bool verbose, RNG& rng) const
    {
        typedef typename property_traits<PosMap>::value_type pos_t;
        typedef typename property_traits<PosMap>::value_type::value_type val_t;

        typedef typename property_traits<VertexWeightMap>::value_type vweight_t;

        pos_t ll(2, numeric_limits<val_t>::max()),
            ur(2, -numeric_limits<val_t>::max());

        vector<pos_t> group_cm;
        vector<vweight_t> group_size;
        vector<size_t> vertices;

        int i, N = num_vertices(g), HN=0;
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v =
                vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            if (pin[v] == 0)
                vertices.push_back(v);
            pos[v].resize(2, 0);
            size_t s = group[v];

            if (s >= group_cm.size())
            {
                group_cm.resize(s + 1);
                group_size.resize(s + 1, 0);
            }
            group_cm[s].resize(2, 0);
            group_size[s] += get(vweight, v);

            for (size_t j = 0; j < 2; ++j)
            {
                ll[j] = min(pos[v][j], ll[j]);
                ur[j] = max(pos[v][j], ur[j]);
                group_cm[s][j] += pos[v][j] * get(vweight, v);
            }
            HN++;
        }

        for (size_t s = 0; s < group_size.size(); ++s)
        {
            if (group_size[s] == 0)
                continue;
            group_cm[s].resize(2, 0);
            for (size_t j = 0; j < 2; ++j)
                group_cm[s][j] /= group_size[s];
        }

        val_t delta = epsilon * K + 1, E = 0, E0;
        E0 = numeric_limits<val_t>::max();
        size_t n_iter = 0;
        val_t step = init_step;
        size_t progress = 0;

        while (delta > epsilon * K && (max_iter == 0 || n_iter < max_iter))
        {
            delta = 0;
            E0 = E;
            E = 0;

            pos_t nll(2, numeric_limits<val_t>::max()),
                nur(2, -numeric_limits<val_t>::max());

            QuadTree<pos_t, vweight_t> qt(ll, ur, max_level);
            for (i = 0; i < N; ++i)
            {
                typename graph_traits<Graph>::vertex_descriptor v =
                    vertex(i, g);
                if (v == graph_traits<Graph>::null_vertex())
                    continue;
                qt.put_pos(pos[v], vweight[v]);
            }

            std::shuffle(vertices.begin(), vertices.end(), rng);

            size_t nmoves = 0;
            N = vertices.size();
            vector<std::reference_wrapper<QuadTree<pos_t, vweight_t>>> Q;
            #pragma omp parallel for default(shared) private(i, Q) \
                reduction(+:E, delta, nmoves) schedule(runtime) if (N > 100)
            for (i = 0; i < N; ++i)
            {
                auto v = vertex(vertices[i], g);

                pos_t diff(2, 0), pos_u(2, 0), ftot(2, 0), cm(2, 0);

                // global repulsive forces
                Q.push_back(std::ref(qt));
                while (!Q.empty())
                {
                    auto& q = Q.back().get();
                    Q.pop_back();

                    if (q.max_level() == 0)
                    {
                        auto& dleafs = q.get_dense_leafs();
                        for(size_t j = 0; j < dleafs.size(); ++j)
                        {
                            val_t d = get_diff(get<0>(dleafs[j]), pos[v],
                                               diff);
                            if (d == 0)
                                continue;
                            val_t f = f_r(C, K, p, pos[v], get<0>(dleafs[j]));
                            f *= get<1>(dleafs[j]) * get(vweight, v);
                            for (size_t l = 0; l < 2; ++l)
                                ftot[l] += f * diff[l];
                        }
                    }
                    else
                    {
                        double w = q.get_w();
                        q.get_cm(cm);
                        double d = get_diff(cm, pos[v], diff);
                        if (w > theta * d)
                        {
                            for(size_t j = 0; j < 4; ++j)
                            {
                                auto& leaf = q.get_leaf(j);
                                if (leaf.get_count() > 0)
                                    Q.push_back(std::ref(leaf));
                            }
                        }
                        else
                        {
                            if (d > 0)
                            {
                                val_t f = f_r(C, K, p, cm, pos[v]);
                                f *= q.get_count() * get(vweight, v);
                                for (size_t l = 0; l < 2; ++l)
                                    ftot[l] += f * diff[l];
                            }
                        }
                    }
                }

                // local attractive forces
                for (auto e : out_edges_range(v, g))
                {
                    auto u = target(e, g);
                    if (u == v)
                        continue;
                    #pragma omp critical
                    {
                        pos_u = pos[u];
                    }
                    get_diff(pos_u, pos[v], diff);
                    val_t f = f_a(K, pos_u, pos[v]);
                    f *= get(eweight, e) * get(vweight, u) * get(vweight, v);
                    for (size_t l = 0; l < 2; ++l)
                        ftot[l] += f * diff[l];
                }

                // inter-group attractive forces
                if (gamma > 0)
                {
                    for (size_t s = 0; s < group_cm.size(); ++s)
                    {
                        if (group_size[s] == 0)
                            continue;
                        if (s == size_t(group[v]))
                            continue;
                        val_t d = get_diff(group_cm[s], pos[v], diff);
                        if (d == 0)
                            continue;
                        double Kp = K * power(HN, 2);
                        val_t f = f_a(Kp, group_cm[s], pos[v]) * gamma * \
                            group_size[s] * get(vweight, v);
                        for (size_t l = 0; l < 2; ++l)
                            ftot[l] += f * diff[l];
                    }
                }

                // inter-group repulsive forces
                if (gamma < 0)
                {
                    for (size_t s = 0; s < group_cm.size(); ++s)
                    {
                        if (group_size[s] == 0)
                            continue;
                        if (s == size_t(group[v]))
                            continue;
                        val_t d = get_diff(group_cm[s], pos[v], diff);
                        if (d == 0)
                            continue;
                        val_t f = f_r(C, K, p, cm, pos[v]);
                        f *= group_size[s] * get(vweight, v) * abs(gamma);
                        for (size_t l = 0; l < 2; ++l)
                            ftot[l] += f * diff[l];
                    }
                }

                // intra-group attractive forces
                if (group_size[group[v]] > 1 && mu > 0)
                {
                    val_t d = get_diff(group_cm[group[v]], pos[v], diff);
                    if (d > 0)
                    {
                        double Kp = K * pow(double(group_size[group[v]]), mu_p);
                        val_t f = f_a(Kp, group_cm[group[v]], pos[v]) * mu * \
                            group_size[group[v]] * get(vweight, v);
                        for (size_t l = 0; l < 2; ++l)
                            ftot[l] += f * diff[l];
                    }
                }

                E += power(norm(ftot), 2);

                #pragma omp critical
                {
                    for (size_t l = 0; l < 2; ++l)
                    {
                        group_cm[group[v]][l] *= group_size[group[v]];
                        group_cm[group[v]][l] -= pos[v][l];

                        ftot[l] *= step;
                        pos[v][l] += ftot[l];

                        nll[l] = min(pos[v][l], nll[l]);
                        nur[l] = max(pos[v][l], nur[l]);

                        group_cm[group[v]][l] += pos[v][l];
                        group_cm[group[v]][l] /= group_size[group[v]];
                    }
                }
                delta += norm(ftot);
                nmoves++;
            }
            n_iter++;
            ll = nll;
            ur = nur;
            delta /= nmoves;

            if (verbose)
                cout << n_iter << " " << E << " " << step << " "
                     << delta << " " << max_level << endl;

            if (simple)
            {
                step *= step_schedule;
            }
            else
            {
                if (E < E0)
                {
                    ++progress;
                    if (progress >= 5)
                    {
                        progress = 0;
                        step /= step_schedule;
                    }
                }
                else
                {
                    progress = 0;
                    step *= step_schedule;
                }
            }
        }
    }
};

} // namespace graph_tool


#endif // GRAPH_FDP_HH
