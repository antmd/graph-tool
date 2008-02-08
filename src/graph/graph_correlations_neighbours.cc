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

#include "graph_filtering.hh"
#include "graph.hh"
#include "histogram.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"

#include "graph_correlations_neighbours.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

// implementations spread across different compile units to minimize memory
// usage during compilation
void graph_correlations_neighbours_imp1(const GraphInterface& g, 
                                        avg_corr_t& avg_corr,
                                        boost::any origin_deg,
                                        boost::any neighbours_deg,
                                        boost::any weight);
void graph_correlations_neighbours_imp2(const GraphInterface& g, 
                                        avg_corr_t& avg_corr,
                                        boost::any origin_deg,
                                        boost::any neighbours_deg,
                                        boost::any weight);
void graph_correlations_neighbours_imp3(const GraphInterface& g, 
                                        avg_corr_t& avg_corr,
                                        boost::any origin_deg,
                                        boost::any neighbours_deg,
                                        boost::any weight);
void graph_correlations_neighbours_imp4(const GraphInterface& g, 
                                        avg_corr_t& avg_corr,
                                        boost::any origin_deg,
                                        boost::any neighbours_deg,
                                        boost::any weight);
void graph_correlations_neighbours_imp5(const GraphInterface& g, 
                                        avg_corr_t& avg_corr,
                                        boost::any origin_deg,
                                        boost::any neighbours_deg,
                                        boost::any weight);
void graph_correlations_neighbours_imp6(const GraphInterface& g, 
                                        avg_corr_t& avg_corr,
                                        boost::any origin_deg,
                                        boost::any neighbours_deg,
                                        boost::any weight);


avg_corr_t
GraphInterface::GetAverageNearestNeighboursCorrelation(deg_t origin_deg,
                                                       deg_t neighbours_deg,
                                                       std::string weight)
    const
{
    avg_corr_t avg_corr;

    typedef ConstantPropertyMap<int,GraphInterface::edge_t> cweight_map_t;

    try
    {        
        any weight_prop;
        if (weight != "")
            weight_prop = prop(weight, _edge_index, _properties);
        else
            weight_prop = cweight_map_t(1);

        run_action<>()
            (*this, get_average_nearest_neighbours_correlation<avg_corr_t>
             (avg_corr), all_selectors(), all_selectors(), 
             mpl::vector<cweight_map_t>())
            (degree_selector(origin_deg, _properties),
             degree_selector(neighbours_deg, _properties), weight_prop);
        graph_correlations_neighbours_imp1
            (*this, avg_corr, degree_selector(origin_deg, _properties),
             degree_selector(neighbours_deg, _properties), weight_prop);
        graph_correlations_neighbours_imp2
            (*this, avg_corr, degree_selector(origin_deg, _properties),
             degree_selector(neighbours_deg, _properties), weight_prop);
        graph_correlations_neighbours_imp3
            (*this, avg_corr, degree_selector(origin_deg, _properties),
             degree_selector(neighbours_deg, _properties), weight_prop);
        graph_correlations_neighbours_imp4
            (*this, avg_corr, degree_selector(origin_deg, _properties),
             degree_selector(neighbours_deg, _properties), weight_prop);
        graph_correlations_neighbours_imp5
            (*this, avg_corr, degree_selector(origin_deg, _properties),
             degree_selector(neighbours_deg, _properties), weight_prop);
        graph_correlations_neighbours_imp6
            (*this, avg_corr, degree_selector(origin_deg, _properties),
             degree_selector(neighbours_deg, _properties), weight_prop);
    }
    catch (dynamic_get_failure &e)
    {
        throw GraphException("error getting scalar property: " +
                             string(e.what()));
    }
    return avg_corr;
}
