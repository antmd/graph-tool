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

#include "graph_community.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

void GraphInterface::GetCommunityStructure(double gamma, comm_corr_t corr, 
                                           size_t n_iter, double Tmin, 
                                           double Tmax, size_t Nspins, 
                                           size_t seed, bool verbose, 
                                           string history_file, string weight, 
                                           string property)
{
    typedef property_map_types::apply<mpl::vector<int32_t,int64_t>,
                                      vertex_index_map_t,
                                      mpl::bool_<false> >::type 
        allowed_spin_properties;

    boost::any vertex_prop;
    bool new_spins;
    try
    {
        find_property_map(_properties, property, typeid(vertex_t));
        vertex_prop = prop(property, _vertex_index, _properties);
        if (!belongs<allowed_spin_properties>()(vertex_prop))
            throw GraphException("vertex property " + property + 
                                 " is not of integer type int32_t or int64_t");
        new_spins = false;
    }
    catch (property_not_found)
    {
        typedef vector_property_map<int32_t,vertex_index_map_t> comm_map_t;
        comm_map_t comm_map(_vertex_index);
        _properties.property(property, comm_map);
        vertex_prop = comm_map;
        new_spins = true;
    }

    boost::any edge_prop;
    typedef DynamicPropertyMapWrap<double,edge_t> weight_map_t;
    typedef ConstantPropertyMap<double,edge_t> no_weight_map_t;
    typedef mpl::vector<weight_map_t,no_weight_map_t> weight_properties;

    if(weight != "")
    {
        try 
        {
            weight_map_t wrap(find_property_map(_properties, weight,
                                                typeid(edge_t)));
            edge_prop = wrap;
        }
        catch (property_not_found)
        {
            throw GraphException("edge property " + weight + " not found");
        }
    }
    else
    {
        edge_prop = no_weight_map_t(1.0);
    }

    bool directed = _directed;
    _directed = false;
    try
    {
        run_action<detail::never_directed>()
            (*this, bind<void>(get_communities_selector(corr),
                               _1, _2, _3, gamma, n_iter, 
                               make_pair(Tmin, Tmax), 
                               make_pair(Nspins, new_spins), 
                               seed, make_pair(verbose,history_file)),
             weight_properties(), allowed_spin_properties())
            (edge_prop, vertex_prop);
        }
    catch (property_not_found& e)
    {
        _directed = directed;
        throw GraphException("error getting scalar property: " +
                             string(e.what()));
    }
    _directed = directed;
}


double GraphInterface::GetModularity(string weight, string property)
{
    double modularity = 0;

    boost::any vertex_prop = prop(property, _vertex_index, _properties);

    boost::any edge_prop;
    typedef ConstantPropertyMap<int32_t,edge_t> weight_map_t;

    if(weight != "")
        edge_prop = prop(weight, _edge_index, _properties);
    else
        edge_prop = weight_map_t(1);
    
    typedef mpl::push_back<edge_scalar_properties, weight_map_t>::type
        weight_properties;

    bool directed = _directed;
    _directed = false;
    try
    {
        run_action<detail::never_directed>()
            (*this, bind<void>(get_modularity(), _1, _2, _3, var(modularity)),
             weight_properties(), vertex_scalar_properties())
            (edge_prop, vertex_prop);
    }
    catch (property_not_found& e)
    {
        _directed = directed;
        throw GraphException("error getting scalar property: " + 
                             string(e.what()));
    }
    _directed = directed;

    return modularity;
}
