// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2013 Tiago de Paula Peixoto <tiago@skewed.de>
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

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_properties.hh"
#include "graph_selectors.hh"

#include <boost/mpl/for_each.hpp>
#include <boost/mpl/contains.hpp>
#include <boost/python/extract.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

template <class Descriptor, class GraphTgt, class GraphSrc, class IndexMap>
struct copy_property_dispatch
{
    copy_property_dispatch(const Descriptor& src_d, const Descriptor& tgt_d,
                           const GraphTgt& tgt, const GraphSrc& src,
                           boost::any& prop_src, boost::any& prop_tgt,
                           IndexMap& index_map, bool& found)
        : src_d(src_d), tgt_d(tgt_d), tgt(tgt), src(src), prop_src(prop_src),
          prop_tgt(prop_tgt), index_map(index_map), found(found) {}


    const Descriptor& src_d;
    const Descriptor& tgt_d;
    const GraphTgt& tgt;
    const GraphSrc& src;
    boost::any& prop_src;
    boost::any& prop_tgt;
    IndexMap& index_map;
    bool& found;

    template <class PropertyMap>
    void operator()(PropertyMap) const
    {
        PropertyMap* psrc = any_cast<PropertyMap>(&prop_src);
        if (psrc == NULL)
            return;
        if (prop_tgt.empty())
            prop_tgt = PropertyMap(index_map);
        PropertyMap* ptgt = any_cast<PropertyMap>(&prop_tgt);
        if (ptgt == NULL)
            return;
        (*ptgt)[tgt_d] = (*psrc)[src_d];
        found = true;
    }
};

template <class PropertyMaps, class Descriptor, class GraphTgt, class GraphSrc,
          class IndexMap>
void copy_property(const Descriptor& src_d, const Descriptor& tgt_d,
                   boost::any& prop_src, boost::any& prop_tgt,
                   const GraphTgt& tgt, const GraphSrc& src,
                   IndexMap& index_map)
{
    bool found = false;
    mpl::for_each<PropertyMaps>(copy_property_dispatch<Descriptor, GraphTgt, GraphSrc, IndexMap>
                                (src_d, tgt_d, tgt, src, prop_src, prop_tgt, index_map,
                                 found));
    if (!found)
        throw ValueException("Cannot find property map type.");
}

struct do_graph_copy
{
    template <class GraphTgt, class GraphSrc, class TgtVertexIndexMap,
              class SrcVertexIndexMap,  class TgtEdgeIndexMap,
              class SrcEdgeIndexMap>
    void operator()(const GraphSrc& src, GraphTgt& tgt,
                    TgtVertexIndexMap src_vertex_index,
                    SrcVertexIndexMap tgt_vertex_index,
                    TgtEdgeIndexMap src_edge_index,
                    SrcEdgeIndexMap tgt_edge_index,
                    vector<pair<reference_wrapper<boost::any>,reference_wrapper<boost::any> > >& vprops,
                    vector<pair<reference_wrapper<boost::any>,reference_wrapper<boost::any> > >& eprops) const
    {
        vector<size_t> index_map(num_vertices(src));
        typename graph_traits<GraphSrc>::vertex_iterator v, v_end;
        for (tie(v, v_end) = vertices(src); v != v_end; ++v)
        {
            if (src_vertex_index[*v] >= index_map.size())
                index_map.resize(src_vertex_index[*v]+1);
            typename graph_traits<GraphTgt>::vertex_descriptor new_v =
                add_vertex(tgt);
            index_map[src_vertex_index[*v]] = tgt_vertex_index[new_v];

            for (size_t i = 0; i < vprops.size(); ++i)
                copy_property<writable_vertex_properties>
                    (*v, new_v, vprops[i].first.get(), vprops[i].second.get(),
                     src, tgt, tgt_vertex_index);
        }

        size_t e_idx = 0;
        typename graph_traits<GraphSrc>::edge_iterator e, e_end;
        for (tie(e, e_end) = edges(src); e != e_end; ++e)
        {
            size_t s = index_map[src_vertex_index[source(*e, src)]];
            size_t t = index_map[src_vertex_index[target(*e, src)]];
            typedef typename graph_traits<GraphTgt>::edge_descriptor edge_t;
            edge_t new_e = add_edge(vertex(s,tgt), vertex(t,tgt), tgt).first;
            tgt_edge_index[new_e] = e_idx;
            ++e_idx;

            for (size_t i = 0; i < eprops.size(); ++i)
                copy_property<writable_edge_properties>
                    (edge_t(*e), new_e, eprops[i].first.get(), eprops[i].second.get(),
                     src, tgt, tgt_edge_index);
        }
    }
};

// copy constructor
GraphInterface::GraphInterface(const GraphInterface& gi, bool keep_ref,
                               python::object ovprops, python::object oeprops)
    :_state(keep_ref ? gi._state : shared_ptr<state_t>(new state_t())),
     _vertex_index(get(vertex_index, _state->_mg)),
     _edge_index(get(edge_index_t(), _state->_mg)),
     _reversed(gi._reversed),
     _directed(gi._directed),
     _vertex_filter_map(_vertex_index),
     _vertex_filter_invert(false),
     _vertex_filter_active(false),
     _edge_filter_map(_edge_index),
     _edge_filter_invert(false),
     _edge_filter_active(false)
{
    if (keep_ref)
        return;

    vector<pair<reference_wrapper<boost::any>,reference_wrapper<boost::any> > > vprops;
    for (int i = 0; i < python::len(ovprops); ++i)
    {
        vprops.push_back(make_pair(ref(python::extract<boost::any&>(ovprops[i][0])()),
                                   ref(python::extract<boost::any&>(ovprops[i][1])())));
    }
    vector<pair<reference_wrapper<boost::any>,reference_wrapper<boost::any> > > eprops;
    for (int i = 0; i < python::len(oeprops); ++i)
    {
        eprops.push_back(make_pair(ref(python::extract<boost::any&>(oeprops[i][0])()),
                                   ref(python::extract<boost::any&>(oeprops[i][1])())));
    }

    run_action<>()
        (const_cast<GraphInterface&>(gi),
         bind<void>(do_graph_copy(), _1, ref(_state->_mg),
                    gi._vertex_index, _vertex_index,
                    gi._edge_index, _edge_index, ref(vprops),
                    ref(eprops)))();

    _state->_nedges = num_edges(_state->_mg);
    _state->_max_edge_index = _state->_nedges > 0 ? num_edges(_state->_mg) - 1 : 0;

    // filters will be copied in python
}
