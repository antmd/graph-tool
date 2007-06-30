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
#include <tr1/unordered_set>
#include <iostream>
#include <iomanip>

#include "graph.hh"
#include "histogram.hh"
#include "graph_filtering.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"
#include <boost/graph/graphviz.hpp>
#include <boost/graph/graphml.hpp>

using namespace std;
using namespace boost;
using namespace boost::lambda;
using namespace graph_tool;

//==============================================================================
// GetCommunityNetwork()
// retrieves the network of communities given a community structure
//==============================================================================

struct get_community_network
{
    template <class Graph, class CommunityMap>
    void operator()(const Graph& g, CommunityMap s, string property, dynamic_properties& edge_properties, string file, string format) const
    {
        tr1::unordered_set<typename boost::property_traits<CommunityMap>::value_type> comms;

        typename graph_traits<Graph>::vertex_iterator v, v_end;
        for (tie(v, v_end) = vertices(g); v != v_end; ++v)
            comms.insert(get(s, *v));

        typedef boost::property<edge_index_t, size_t> EdgeProperty;

        typedef adjacency_list <vecS, // edges
                                vecS, // vertices
                                typename mpl::if_<typename is_convertible<typename graph_traits<Graph>::directed_category, undirected_tag>::type,
                                                  bidirectionalS,
                                                  undirectedS>::type,
                                no_property,
                                EdgeProperty > comm_graph_t; 

        comm_graph_t comm_graph;

        typedef typename property_map<comm_graph_t, vertex_index_t>::type comm_vertex_index_map_t;
        comm_vertex_index_map_t comm_vertex_index(get(vertex_index, comm_graph));

        vector_property_map<typename property_traits<CommunityMap>::value_type, comm_vertex_index_map_t> comm_map(comm_vertex_index);
        tr1::unordered_map<typename property_traits<CommunityMap>::value_type, typename graph_traits<comm_graph_t>::vertex_descriptor> comm_vertex_map;
        for (typeof(comms.begin()) iter = comms.begin(); iter != comms.end(); ++iter)
        {
            comm_vertex_map[*iter] = add_vertex(comm_graph);
            comm_map[comm_vertex_map[*iter]] = *iter;
        }

        dynamic_properties dp;
        dp.property(property, comm_map);

        typedef typename property_map<comm_graph_t,edge_index_t>::type comm_edge_index_map_t;
        comm_edge_index_map_t comm_edge_index(get(edge_index_t(), comm_graph));

        typedef HashedDescriptorMap<comm_edge_index_map_t, typename graph_traits<Graph>::edge_descriptor> edge_map_t;
        edge_map_t edge_map(comm_edge_index);

        size_t e_index = 0;
        typename graph_traits<Graph>::edge_iterator e, e_end;
        for (tie(e, e_end) = edges(g); e != e_end; ++e)
        {
            typename graph_traits<comm_graph_t>::edge_descriptor edge;
            edge  = add_edge(comm_vertex_map[get(s,source(*e, g))], comm_vertex_map[get(s,target(*e, g))], comm_graph).first;
            edge_map[edge] = *e;
            comm_edge_index[edge] = e_index++;
        }

        for (typeof(edge_properties.begin()) iter = edge_properties.begin(); iter != edge_properties.end(); ++iter)
            dp.insert(iter->first, auto_ptr<dynamic_property_map>(new dynamic_property_map_wrap<edge_map_t>(edge_map, *iter->second)));

        bool graphviz = false;
        if (format == "")
            graphviz = ends_with(file,".dot") || ends_with(file,".dot.gz") || ends_with(file,".dot.bz2");
        else if (format == "dot")
            graphviz = true;
        else if (format != "xml")
            throw GraphException("error writing to file '" + file + "': requested invalid format '" + format + "'");
        try
        {
            iostreams::filtering_stream<iostreams::output> stream;
            ofstream file_stream;
            if (file == "-")
                stream.push(cout);
            else
            {
                file_stream.open(file.c_str(), ios_base::out | ios_base::binary);
                file_stream.exceptions(ios_base::badbit | ios_base::failbit);
                if (ends_with(file,".gz"))
                    stream.push(iostreams::gzip_compressor());
                if (ends_with(file,".bz2"))
                    stream.push(iostreams::bzip2_compressor());
                stream.push(file_stream);
            }
            stream.exceptions(ios_base::badbit | ios_base::failbit);
            
            if (graphviz)
            {
                dp.property("vertex_id", comm_vertex_index);
                write_graphviz(stream, comm_graph, dp, string("vertex_id"));
            }
            else
            {
                write_graphml(stream, comm_graph, comm_vertex_index, dp, true);
            }
            stream.reset();
        }
        catch (ios_base::failure &e)
        {
            throw GraphException("error writing to file '" + file + "':" + e.what());
        }
    }

    template <class EdgeMap>
    class dynamic_property_map_wrap: public dynamic_property_map
    {
    public:
        dynamic_property_map_wrap(EdgeMap& edge_map, dynamic_property_map& dm): _edge_map(edge_map), _dm(dm) {}
        any get(const any& key)
        {
            return _dm.get(_edge_map[any_cast<typename property_traits<EdgeMap>::key_type>(key)]);
        }

        string get_string(const any& key)
        {
            return _dm.get_string(_edge_map[any_cast<typename property_traits<EdgeMap>::key_type>(key)]);
        }

        void put(const any& key, const any& value)
        {
            return _dm.put(_edge_map[any_cast<typename property_traits<EdgeMap>::key_type>(key)], value);
        }

        const type_info& key() const
        {
            return typeid(typename property_traits<EdgeMap>::key_type);
        }

        const type_info& value() const
        {
            return _dm.value();
        }
    private:
        EdgeMap& _edge_map;
        dynamic_property_map& _dm;
    };

};

void GraphInterface::GetCommunityNetwork(string property, string out_file, string format) const
{
    try 
    {
        typedef DynamicPropertyMapWrap<double, graph_traits<multigraph_t>::vertex_descriptor> comm_map_t;
        comm_map_t comm_map(find_property_map(_properties, property, typeid(graph_traits<multigraph_t>::vertex_descriptor)));
        
        dynamic_properties_copy edge_properties;
        for (typeof(_properties.begin()) iter = _properties.begin(); iter != _properties.end(); ++iter)
            if (iter->second->key() == typeid(graph_traits<multigraph_t>::edge_descriptor))
                edge_properties.insert(iter->first, auto_ptr<dynamic_property_map>(iter->second));

        check_filter(*this, bind<void>(get_community_network(), _1, var(comm_map), property, var(edge_properties), out_file, format), reverse_check(), directed_check());
    }
    catch (property_not_found) 
    {
        throw GraphException("edge property " + property + " not found");
    }
}
