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
// GetLineGraph()
// retrieves the line graph
//==============================================================================

struct get_line_graph
{
    template <class Graph, class EdgeIndexMap>
    void operator()(const Graph& g, EdgeIndexMap edge_index, dynamic_properties& properties,  string file, string format) const
    {

        typedef boost::property<edge_index_t, size_t> EdgeProperty;

        typedef adjacency_list <vecS, // edges
                                vecS, // vertices
                                undirectedS,
                                no_property,
                                EdgeProperty > line_graph_t; 

        line_graph_t line_graph;

        typedef typename property_map<line_graph_t, vertex_index_t>::type line_vertex_index_map_t;
        line_vertex_index_map_t line_vertex_index(get(vertex_index, line_graph));

        typedef HashedDescriptorMap<line_vertex_index_map_t, typename graph_traits<typename Graph::original_graph_t>::edge_descriptor> edge_map_t;
        edge_map_t edge_map(line_vertex_index);

        typedef HashedDescriptorMap<EdgeIndexMap, typename graph_traits<line_graph_t>::vertex_descriptor> edge_to_vertex_map_t;
        edge_to_vertex_map_t edge_to_vertex_map(edge_index);

        typename graph_traits<Graph>::edge_iterator e, e_end;
        for (tie(e, e_end) = edges(g); e != e_end; ++e)
        {
            typename graph_traits<line_graph_t>::vertex_descriptor v;
            v = add_vertex(line_graph);
            edge_to_vertex_map[*e] = v;
            edge_map[v] = *e;
        }

        typedef typename property_map<line_graph_t,edge_index_t>::type line_edge_index_map_t;
        line_edge_index_map_t line_edge_index(get(edge_index_t(), line_graph));

        typedef HashedDescriptorMap<line_edge_index_map_t, typename graph_traits<Graph>::vertex_descriptor> vertex_map_t;
        vertex_map_t vertex_map(line_edge_index);

        size_t e_index = 0;
        typename graph_traits<Graph>::vertex_iterator v, v_end;
        for (tie(v, v_end) = vertices(g); v != v_end; ++v)
        {
            typename graph_traits<Graph>::out_edge_iterator e1, e2, e_end;
            for (tie(e1, e_end) = out_edges(*v, g); e1 != e_end; ++e1)
                for (e2 = e1; e2 != e_end; ++e2)
                    if (*e1 != *e2) 
                    {
                        typename graph_traits<line_graph_t>::edge_descriptor new_edge;
                        new_edge = add_edge(edge_to_vertex_map[*e1], edge_to_vertex_map[*e2], line_graph).first;
                        line_edge_index[new_edge] = e_index++;
                        vertex_map[new_edge] = *v;
                    }
        }

        dynamic_properties dp;
        for (typeof(properties.begin()) iter = properties.begin(); iter != properties.end(); ++iter)
        {
            if (iter->second->key() != typeid(graph_property_tag))
            {
                if (iter->second->key() == typeid(typename graph_traits<Graph>::vertex_descriptor))
                    dp.insert(iter->first, auto_ptr<dynamic_property_map>(new dynamic_property_map_wrap<vertex_map_t>(vertex_map, *iter->second)));
                else 
                    dp.insert(iter->first, auto_ptr<dynamic_property_map>(new dynamic_property_map_wrap<edge_map_t>(edge_map, *iter->second)));
            }
            else
            {
                dp.insert(iter->first, auto_ptr<dynamic_property_map>(iter->second));
            }
        }

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
                dp.property("vertex_id", line_vertex_index);
                write_graphviz(stream, line_graph, dp, string("vertex_id"));
            }
            else
            {
                write_graphml(stream, line_graph, line_vertex_index, dp, true);
            }
            stream.reset();
        }
        catch (ios_base::failure &e)
        {
            throw GraphException("error writing to file '" + file + "':" + e.what());
        }

        for (typeof(dp.begin()) iter = dp.begin(); iter != dp.end(); ++iter)
            if (iter->second->key() == typeid(graph_property_tag))
                iter->second = 0;
    }

    template <class DescriptorMap>
    class dynamic_property_map_wrap: public dynamic_property_map
    {
    public:
        dynamic_property_map_wrap(DescriptorMap& edge_map, dynamic_property_map& dm): _descriptor_map(edge_map), _dm(dm) {}
        any get(const any& key)       
        {
            typename property_traits<DescriptorMap>::key_type k = any_cast<typename property_traits<DescriptorMap>::key_type>(key);
            return _dm.get(_descriptor_map[k]);
        }

        string get_string(const any& key)
        {
            typename property_traits<DescriptorMap>::key_type k = any_cast<typename property_traits<DescriptorMap>::key_type>(key);
            return _dm.get_string(_descriptor_map[k]);
        }

        void put(const any& key, const any& value)
        {
            return _dm.put(_descriptor_map[any_cast<typename property_traits<DescriptorMap>::key_type>(key)], value);
        }

        const type_info& key() const
        {
            return typeid(typename property_traits<DescriptorMap>::key_type);
        }

        const type_info& value() const
        {
            return _dm.value();
        }
    private:
        DescriptorMap& _descriptor_map;
        dynamic_property_map& _dm;
    };

};

void GraphInterface::GetLineGraph(string out_file, string format)
{
    bool directed = _directed;
    _directed = false;
    check_filter(*this, bind<void>(get_line_graph(), _1, var(_edge_index), var(_properties), out_file, format), reverse_check(), always_undirected());
    _directed = directed;
}
