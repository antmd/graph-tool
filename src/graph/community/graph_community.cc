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

#include "graph_python_interface.hh"
#include "graph_filtering.hh"
#include "graph.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"

#include <boost/mpl/push_back.hpp>
#include <boost/python.hpp>

#include "graph_community.hh"

using namespace std;
using namespace boost;

using namespace graph_tool;


void community_structure(GraphInterface& g, double gamma, string corr_name,
                         size_t n_iter, double Tmin, double Tmax, size_t Nspins,
                         rng_t& rng, bool verbose, string history_file,
                         boost::any weight, boost::any property)
{
    typedef property_map_types::apply<boost::mpl::vector<int32_t,int64_t>,
                                      GraphInterface::vertex_index_map_t,
                                      boost::mpl::bool_<false> >::type
        allowed_spin_properties;

    if (!belongs<allowed_spin_properties>()(property))
        throw ValueException("vertex property is not of integer type int32_t "
                             "or int64_t");

    typedef DynamicPropertyMapWrap<double,GraphInterface::edge_t> weight_map_t;
    typedef ConstantPropertyMap<double,GraphInterface::edge_t> no_weight_map_t;
    typedef boost::mpl::vector<weight_map_t,no_weight_map_t> weight_properties;

    if (weight.empty())
        weight = no_weight_map_t(1.0);
    else
        weight = weight_map_t(weight, edge_scalar_properties());

    comm_corr_t corr;
    if (corr_name ==  "erdos")
        corr = ERDOS_REYNI;
    else if (corr_name == "uncorrelated")
        corr = UNCORRELATED;
    else if (corr_name == "correlated")
        corr = CORRELATED;
    else
        throw ValueException("invalid correlation type: " + corr_name);

    run_action<graph_tool::detail::never_directed>()
        (g, std::bind(get_communities_selector(corr, g.GetVertexIndex()),
                      placeholders::_1, placeholders::_2, placeholders::_3, gamma, n_iter,
                      make_pair(Tmin, Tmax), Nspins,
                      std::ref(rng), make_pair(verbose,history_file)),
         weight_properties(), allowed_spin_properties())
        (weight, property);
}


double modularity(GraphInterface& g, boost::any weight, boost::any property)
{
    double modularity = 0;

    typedef ConstantPropertyMap<int32_t,GraphInterface::edge_t> weight_map_t;
    typedef boost::mpl::push_back<edge_scalar_properties, weight_map_t>::type
        edge_props_t;

    if(weight.empty())
        weight = weight_map_t(1);

    run_action<graph_tool::detail::never_directed>()
        (g, std::bind(get_modularity(), placeholders::_1, placeholders::_2,
                      placeholders::_3, std::ref(modularity)),
         edge_props_t(), vertex_scalar_properties())
        (weight, property);
    return modularity;
}

using namespace boost::python;


extern void community_network(GraphInterface& gi, GraphInterface& cgi,
                              boost::any community_property,
                              boost::any condensed_community_property,
                              boost::any vertex_count,
                              boost::any edge_count, boost::any vweight,
                              boost::any eweight, bool self_loops);

void community_network_vavg(GraphInterface& gi, GraphInterface& cgi,
                            boost::any community_property,
                            boost::any condensed_community_property,
                            boost::any vertex_count,
                            boost::any vweight,
                            boost::python::list avprops);

void community_network_eavg(GraphInterface& gi, GraphInterface& cgi,
                            boost::any community_property,
                            boost::any condensed_community_property,
                            boost::any edge_count,
                            boost::any eweight,
                            boost::python::list aeprops,
                            bool self_loops);


extern void export_blockmodel();
extern void export_blockmodel_overlap();

BOOST_PYTHON_MODULE(libgraph_tool_community)
{
    def("community_structure", &community_structure);
    def("modularity", &modularity);
    def("community_network", &community_network);
    def("community_network_vavg", &community_network_vavg);
    def("community_network_eavg", &community_network_eavg);

    export_blockmodel();
    export_blockmodel_overlap();
}
