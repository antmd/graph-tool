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
#include "graph_selectors.hh"
#include "graph_properties.hh"

#include <boost/lambda/bind.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/python.hpp>

#include "graph_community.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;

using namespace graph_tool;


void community_structure(GraphInterface& g, double gamma, string corr_name,
                         size_t n_iter, double Tmin, double Tmax, size_t Nspins,
                         bool new_spins, size_t seed, bool verbose,
                         string history_file, boost::any weight,
                         boost::any property)
{
    using boost::lambda::bind;
    using boost::lambda::_1;
    using boost::lambda::_2;
    using boost::lambda::_3;

    typedef property_map_types::apply<mpl::vector<int32_t,int64_t>,
                                      GraphInterface::vertex_index_map_t,
                                      mpl::bool_<false> >::type
        allowed_spin_properties;

    if (!belongs<allowed_spin_properties>()(property))
        throw ValueException("vertex property is not of integer type int32_t "
                             "or int64_t");

    typedef DynamicPropertyMapWrap<double,GraphInterface::edge_t> weight_map_t;
    typedef ConstantPropertyMap<double,GraphInterface::edge_t> no_weight_map_t;
    typedef mpl::vector<weight_map_t,no_weight_map_t> weight_properties;

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

    bool directed = g.GetDirected();
    g.SetDirected(false);
    run_action<graph_tool::detail::never_directed>()
        (g, bind<void>(get_communities_selector(corr, g.GetVertexIndex()),
                       _1, _2, _3, gamma, n_iter,
                       make_pair(Tmin, Tmax), Nspins,
                       seed, make_pair(verbose,history_file)),
         weight_properties(), allowed_spin_properties())
        (weight, property);
    g.SetDirected(directed);
}


double modularity(GraphInterface& g, boost::any weight, boost::any property)
{
    using boost::lambda::bind;
    using boost::lambda::_1;
    using boost::lambda::_2;
    using boost::lambda::_3;

    double modularity = 0;

    typedef ConstantPropertyMap<int32_t,GraphInterface::edge_t> weight_map_t;
    typedef mpl::push_back<edge_scalar_properties, weight_map_t>::type
        edge_props_t;

    if(weight.empty())
        weight = weight_map_t(1);

    bool directed = g.GetDirected();
    g.SetDirected(false);
    run_action<graph_tool::detail::never_directed>()
        (g, bind<void>(get_modularity(), _1, _2, _3, var(modularity)),
         edge_props_t(), vertex_scalar_properties())
        (weight, property);
    g.SetDirected(directed);

    return modularity;
}

using namespace boost::python;


extern void community_network(GraphInterface& gi, GraphInterface& cgi,
                              boost::any community_property,
                              boost::any vertex_count,
                              boost::any edge_count, boost::any weight);

BOOST_PYTHON_MODULE(libgraph_tool_community)
{
    def("community_structure", &community_structure);
    def("modularity", &modularity);
    def("community_network", &community_network);
}
