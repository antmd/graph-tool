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

// based on code written by Alexandre Hannud Abdo <abdo@member.fsf.org>

#include "graph_filtering.hh"
#include "graph.hh"
#include "histogram.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"

#include <boost/graph/breadth_first_search.hpp>

#include "graph_extended_clustering.hh"

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

template <class PropertySequence>
struct prop_vector
{
    boost::any operator()(const vector<boost::any>& props, bool& found) const
    {
        boost::any prop_vec;
        mpl::for_each<PropertySequence>(bind<void>(get_prop_vector(), _1,
                                                   var(props), var(prop_vec),
                                                   var(found)));
        return prop_vec;
    }

    struct get_prop_vector
    {
        template <class Property>
            void operator()(Property, const vector<boost::any>& props,
                            boost::any& prop_vec, bool& found) const
        {
            if (typeid(Property) == props[0].type())
            {
                try
                {
                    vector<Property> vec;
                    vec.resize(props.size());
                    for (size_t i = 0; i < props.size(); ++i)
                        vec[i] = any_cast<Property>(props[i]);
                    prop_vec = vec;
                    found = true;
                }
                catch (bad_any_cast)
                {
                    found = false;
                }
            }
        }
    };
};


struct get_property_vector_type
{
    template <class Property>
    struct apply
    {
        typedef vector<Property> type;
    };
};

void GraphInterface::SetExtendedClusteringToProperty(string property_prefix,
                                                     size_t max_depth)
{
    typedef vector_property_map<double, vertex_index_map_t> cmap_t;
    vector<any> cmaps(max_depth);
    for (size_t i = 0; i < cmaps.size(); ++i)
    {
        string name = property_prefix + lexical_cast<string>(i+1);
        try
        {
            cmaps[i] = prop(name, _vertex_index, _properties);
        }
        catch (property_not_found)
        {
            cmap_t cmap(num_vertices(_mg), _vertex_index);
            _properties.property(name, cmap);
            cmaps[i] = cmap;
        }
    }

    typedef mpl::transform<vertex_floating_properties,
                           get_property_vector_type>::type
        property_vectors;

    bool found = false;

    run_action<>()
        (*this, bind<void>(get_extended_clustering(), _1, _vertex_index,_2),
         property_vectors())
        (prop_vector<vertex_scalar_properties>()(cmaps, found));

    if (!found)
        throw GraphException("All vertex properties " + property_prefix +
                             "* must be of the same floating point type!");
}
