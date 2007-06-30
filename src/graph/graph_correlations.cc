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

#include <algorithm>
#include <tr1/unordered_set>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/random.hpp>

#include "graph.hh"
#include "histogram.hh"
#include "graph_filtering.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"
#include "shared_map.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

//==============================================================================
// get_correlation_histogram
// retrieves the generalized vertex-vertex correlation histogram                
//==============================================================================

template <class DegreeSelector1, class DegreeSelector2>
struct get_correlation_histogram
{
    get_correlation_histogram(DegreeSelector1& deg1, DegreeSelector2& deg2): _deg1(deg1), _deg2(deg2) {}

    template <class Graph, class WeightMap, class Hist>
    void operator()(Graph& g, WeightMap weight, Hist& hist) const
    {
        SharedMap<Hist> s_hist(hist);

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) firstprivate(s_hist) schedule(dynamic)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;

            typename Hist::key_type key;
            key.first = _deg1(v, g);
            typename graph_traits<Graph>::out_edge_iterator e, e_end;
            for (tie(e,e_end) = out_edges(v, g); e != e_end; ++e)
            {
                key.second = _deg2(target(*e,g),g);
                s_hist[key] += typename Hist::value_type::second_type(get(weight, *e));   
            }
        }
        s_hist.Gather();
    }
    DegreeSelector1& _deg1;
    DegreeSelector2& _deg2;
};

//==============================================================================
// GetDegreeCorrelationHistogram(deg1, deg2)
// retrieves the degree correlation histogram                
//==============================================================================

template <class WeightMap, class SecondDegreeSelectors>
struct choose_vertex_correlation_histogram
{
    choose_vertex_correlation_histogram(const GraphInterface &g, WeightMap weight, GraphInterface::deg_t deg1, 
                                        GraphInterface::deg_t deg2, GraphInterface::hist2d_t &hist)
        : _g(g), _weight(weight), _hist(hist) 
    {
        tie(_deg1, _deg_name1) = get_degree_type(deg1);
        tie(_deg2, _deg_name2) = get_degree_type(deg2);
    }

    template <class DegreeSelector1>
    struct check_second_degree
    {
        check_second_degree(choose_vertex_correlation_histogram<WeightMap,SecondDegreeSelectors> &parent):_parent(parent) {}
        template <class DegreeSelector2>
        void operator()(DegreeSelector2)
        {
            if (mpl::at<degree_selector_index,DegreeSelector2>::type::value == _parent._deg2)
            {
                DegreeSelector1 deg1(_parent._deg_name1, _parent._g);
                DegreeSelector2 deg2(_parent._deg_name2, _parent._g);
                check_filter(_parent._g, bind<void>(get_correlation_histogram<DegreeSelector1,DegreeSelector2>(deg1,deg2), _1, var(_parent._weight), var(_parent._hist)), 
                             reverse_check(),directed_check());
            }
        }
        choose_vertex_correlation_histogram<WeightMap,SecondDegreeSelectors> _parent;
    };

    template <class DegreeSelector>
    void operator()(DegreeSelector)
    {
        if (mpl::at<degree_selector_index,DegreeSelector>::type::value == _deg1)
            mpl::for_each<SecondDegreeSelectors>(check_second_degree<DegreeSelector>(*this));
    }

    const GraphInterface &_g;
    WeightMap _weight;
    GraphInterface::hist2d_t &_hist;
    GraphInterface::degree_t _deg1;
    string _deg_name1;
    GraphInterface::degree_t _deg2;
    string _deg_name2;
    
};

GraphInterface::hist2d_t 
GraphInterface::GetVertexCorrelationHistogram(GraphInterface::deg_t deg1, GraphInterface::deg_t deg2, string weight) const
{
    hist2d_t hist;
    typedef mpl::vector<in_degreeS, out_degreeS, total_degreeS, scalarS> degree_selectors;

    try 
    {
        if(weight != "")
        {
            try 
            {
                // FIXME: it would be good also to have a version for a static map (vector_property_map), 
                //        but adding this makes GCC use more than 1 GB of RAM in my system.

                dynamic_property_map& weight_prop = find_property_map(_properties, weight, typeid(graph_traits<multigraph_t>::edge_descriptor));
                typedef DynamicPropertyMapWrap<double,graph_traits<multigraph_t>::edge_descriptor> weight_map_t;
                weight_map_t weight_map(weight_prop);
                mpl::for_each<degree_selectors>(choose_vertex_correlation_histogram<weight_map_t, degree_selectors>(*this, weight_map, deg1, deg2, hist));
            }
            catch (property_not_found& e)
            {
                throw GraphException("error getting scalar property: " + string(e.what()));
            }
        }
        else
        {
            typedef ConstantPropertyMap<double,graph_traits<multigraph_t>::edge_descriptor>  weight_map_t;
            weight_map_t weight_map(1.0);
            mpl::for_each<degree_selectors>(choose_vertex_correlation_histogram<weight_map_t, degree_selectors>(*this, weight_map, deg1, deg2, hist));
        }
    }
    catch (dynamic_get_failure &e)
    {
        throw GraphException("error getting scalar property: " + string(e.what()));
    }

    return hist;
}
