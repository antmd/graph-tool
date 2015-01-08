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
    boost::mpl::for_each<PropertyMaps>(copy_property_dispatch<Descriptor, GraphTgt, GraphSrc, IndexMap>
                                       (src_d, tgt_d, tgt, src, prop_src, prop_tgt, index_map,
                                        found));
    if (!found)
        throw ValueException("Cannot find property map type.");
}

template <class GraphSrc, class GraphTgt, class IndexMap, class SrcIndexMap,
          class TgtIndexMap>
struct copy_vertex_property_dispatch
{
    copy_vertex_property_dispatch(const GraphSrc& src, const GraphTgt& tgt,
                                  boost::any& prop_src, boost::any& prop_tgt,
                                  IndexMap& index_map,
                                  SrcIndexMap& src_vertex_index,
                                  TgtIndexMap& tgt_vertex_index,
                                  bool& found)
        : src(src), tgt(tgt), prop_src(prop_src),
          prop_tgt(prop_tgt), index_map(index_map),
          src_vertex_index(src_vertex_index),
          tgt_vertex_index(tgt_vertex_index), found(found) {}


    const GraphSrc& src;
    const GraphTgt& tgt;
    boost::any& prop_src;
    boost::any& prop_tgt;
    IndexMap& index_map;
    SrcIndexMap& src_vertex_index;
    TgtIndexMap& tgt_vertex_index;
    bool& found;

    template <class PropertyMap>
    void operator()(PropertyMap) const
    {
        PropertyMap* psrc = any_cast<PropertyMap>(&prop_src);
        if (psrc == NULL)
            return;
        if (prop_tgt.empty())
            prop_tgt = PropertyMap(tgt_vertex_index);
        PropertyMap* ptgt = any_cast<PropertyMap>(&prop_tgt);
        if (ptgt == NULL)
            return;
        found = true;

        typename PropertyMap::unchecked_t p_src =
            psrc->get_unchecked(num_vertices(src));
        typename PropertyMap::unchecked_t p_tgt =
            ptgt->get_unchecked(num_vertices(tgt));

        int i, N = num_vertices(src);
        #pragma omp parallel for default(shared) private(i) schedule(runtime) if (N > 100)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<GraphSrc>::vertex_descriptor v = vertex(i, src);
            if (v == graph_traits<GraphSrc>::null_vertex())
                continue;
            typename graph_traits<GraphTgt>::vertex_descriptor new_v =
                vertex(index_map[i], tgt);
            p_tgt[new_v] = p_src[v];
        }
    }
};

template <class PropertyMaps, class GraphSrc, class GraphTgt,
          class IndexMap, class SrcIndexMap, class TgtIndexMap>
void copy_vertex_property(boost::any& prop_src, boost::any& prop_tgt,
                          const GraphSrc& src, const GraphTgt& tgt,
                          IndexMap& index_map, SrcIndexMap& src_vertex_index,
                          TgtIndexMap& tgt_vertex_index)
{
    bool found = false;
    boost::mpl::for_each<PropertyMaps>(copy_vertex_property_dispatch<GraphSrc, GraphTgt,
                                                                     IndexMap, SrcIndexMap,
                                                                     TgtIndexMap>
            (src, tgt, prop_src, prop_tgt, index_map, src_vertex_index,
             tgt_vertex_index, found));
    if (!found)
        throw ValueException("Cannot find property map type.");
}


template <class GraphSrc, class GraphTgt, class IndexMap, class SrcIndexMap>
struct copy_edge_property_dispatch
{
    copy_edge_property_dispatch(const GraphSrc& src, const GraphTgt& tgt,
                                boost::any& prop_src, boost::any& prop_tgt,
                                IndexMap& index_map,
                                SrcIndexMap& src_edge_index,
                                size_t max_src_edge_index,
                                bool& found)
        : src(src), tgt(tgt), prop_src(prop_src),
          prop_tgt(prop_tgt), index_map(index_map),
          src_edge_index(src_edge_index),
          max_src_edge_index(max_src_edge_index), found(found) {}


    const GraphSrc& src;
    const GraphTgt& tgt;
    boost::any& prop_src;
    boost::any& prop_tgt;
    IndexMap& index_map;
    SrcIndexMap& src_edge_index;
    size_t max_src_edge_index;
    bool& found;

    template <class PropertyMap>
    void operator()(PropertyMap) const
    {
        PropertyMap* psrc = any_cast<PropertyMap>(&prop_src);
        if (psrc == NULL)
            return;
        if (prop_tgt.empty())
            prop_tgt = PropertyMap(get(edge_index_t(), tgt));
        PropertyMap* ptgt = any_cast<PropertyMap>(&prop_tgt);
        if (ptgt == NULL)
            return;
        found = true;

        typename PropertyMap::unchecked_t p_src =
            psrc->get_unchecked(max_src_edge_index + 1);
        typename PropertyMap::unchecked_t p_tgt =
            ptgt->get_unchecked(num_edges(tgt));

        int i, N = num_vertices(src);
        #pragma omp parallel for default(shared) private(i) schedule(runtime) if (N > 100)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<GraphSrc>::vertex_descriptor v = vertex(i, src);
            if (v == graph_traits<GraphSrc>::null_vertex())
                continue;

            typename graph_traits<GraphSrc>::out_edge_iterator e, e_end;
            for (tie(e, e_end) = out_edges(v, src); e != e_end; ++e)
            {
                typename graph_traits<GraphSrc>::vertex_descriptor s = source(*e, src);
                typename graph_traits<GraphSrc>::vertex_descriptor t = target(*e, src);
                if (!is_directed::apply<GraphSrc>::type::value && s > t)
                    continue;
                size_t ei = src_edge_index[*e];
                typename graph_traits<GraphTgt>::edge_descriptor new_e = index_map[ei];
                p_tgt[new_e] = p_src[*e];
            }
        }
    }
};

template <class PropertyMaps, class GraphSrc, class GraphTgt,
          class IndexMap, class SrcIndexMap>
void copy_edge_property(boost::any& prop_src, boost::any& prop_tgt,
                        const GraphSrc& src, const GraphTgt& tgt,
                        IndexMap& index_map, SrcIndexMap& src_vertex_index,
                        size_t max_src_edge_index)
{
    bool found = false;
    boost::mpl::for_each<PropertyMaps>(copy_edge_property_dispatch<GraphSrc, GraphTgt,
                                                                   IndexMap, SrcIndexMap>
            (src, tgt, prop_src, prop_tgt, index_map, src_vertex_index,
             max_src_edge_index, found));
    if (!found)
        throw ValueException("Cannot find property map type.");
}

struct do_graph_copy
{

    do_graph_copy(size_t max_src_edge_index)
        : max_src_edge_index(max_src_edge_index) {}
    size_t max_src_edge_index;

    template <class GraphTgt, class GraphSrc, class TgtVertexIndexMap,
              class SrcVertexIndexMap,  class TgtEdgeIndexMap,
              class SrcEdgeIndexMap, class OrderMap>
    void operator()(const GraphSrc& src, GraphTgt& tgt,
                    TgtVertexIndexMap src_vertex_index,
                    SrcVertexIndexMap tgt_vertex_index,
                    TgtEdgeIndexMap,
                    SrcEdgeIndexMap src_edge_index,
                    OrderMap vertex_order,
                    vector<pair<std::reference_wrapper<boost::any>,std::reference_wrapper<boost::any>>>& vprops,
                    vector<pair<std::reference_wrapper<boost::any>,std::reference_wrapper<boost::any>>>& eprops) const
    {
        vector<size_t> index_map(num_vertices(src));
        typename graph_traits<GraphSrc>::vertex_iterator v, v_end;
        for (tie(v, v_end) = vertices(src); v != v_end; ++v)
        {
            if (src_vertex_index[*v] >= index_map.size())
                index_map.resize(src_vertex_index[*v]+1);
            typename graph_traits<GraphTgt>::vertex_descriptor new_v =
                get(vertex_order, *v);
            while (new_v >= num_vertices(tgt))
                add_vertex(tgt);
            index_map[src_vertex_index[*v]] = tgt_vertex_index[new_v];
        }

        for (size_t i = 0; i < vprops.size(); ++i)
            copy_vertex_property<writable_vertex_properties>
                    (vprops[i].first.get(), vprops[i].second.get(),
                     src, tgt, index_map, src_vertex_index, tgt_vertex_index);


        vector<typename graph_traits<GraphTgt>::edge_descriptor> edge_map(num_edges(src));
        typename graph_traits<GraphSrc>::edge_iterator e, e_end;
        for (tie(e, e_end) = edges(src); e != e_end; ++e)
        {
            size_t s = index_map[src_vertex_index[source(*e, src)]];
            size_t t = index_map[src_vertex_index[target(*e, src)]];
            typedef typename graph_traits<GraphTgt>::edge_descriptor edge_t;
            edge_t new_e = add_edge(vertex(s,tgt), vertex(t,tgt), tgt).first;

            size_t ei = src_edge_index[*e];
            if (ei >= edge_map.size())
                edge_map.resize(ei + 1);

            edge_map[ei] = new_e;
        }

        for (size_t i = 0; i < eprops.size(); ++i)
            copy_edge_property<writable_edge_properties>
                (eprops[i].first.get(), eprops[i].second.get(),
                 src, tgt, edge_map, src_edge_index, max_src_edge_index);

    }
};

// copy constructor
GraphInterface::GraphInterface(const GraphInterface& gi, bool keep_ref,
                               python::object ovprops, python::object oeprops,
                               python::object vorder)
    :_mg(keep_ref ? gi._mg : std::shared_ptr<multigraph_t>(new multigraph_t())),
     _vertex_index(get(vertex_index, *_mg)),
     _edge_index(get(edge_index_t(), *_mg)),
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

    if (vorder == python::object())
    {
        // simple copying
        *_mg = *gi._mg;
        return;
    }

    vector<pair<std::reference_wrapper<boost::any>,std::reference_wrapper<boost::any>>> vprops;
    for (int i = 0; i < python::len(ovprops); ++i)
    {
        vprops.push_back(make_pair(std::ref(python::extract<boost::any&>(ovprops[i][0])()),
                                   std::ref(python::extract<boost::any&>(ovprops[i][1])())));
    }
    vector<pair<std::reference_wrapper<boost::any>,std::reference_wrapper<boost::any>>> eprops;
    for (int i = 0; i < python::len(oeprops); ++i)
    {
        eprops.push_back(make_pair(std::ref(python::extract<boost::any&>(oeprops[i][0])()),
                                   std::ref(python::extract<boost::any&>(oeprops[i][1])())));
    }

    boost::any avorder = python::extract<boost::any>(vorder)();
    run_action<>()
        (const_cast<GraphInterface&>(gi),
         std::bind(do_graph_copy(gi._mg->get_last_index()),
                   std::placeholders::_1, std::ref(*_mg),
                   gi._vertex_index, _vertex_index, gi._edge_index,
                   _edge_index, std::placeholders::_2, std::ref(vprops),
                   std::ref(eprops)),
         vertex_scalar_properties())(avorder);
    // filters will be copied in python
}
