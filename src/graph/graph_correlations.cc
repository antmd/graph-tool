// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006  Tiago de Paula Peixoto <tiago@forked.de>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
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

    template <class Graph, class Hist>
    void operator()(Graph &g, Hist &hist) const
    {
        typename graph_traits<Graph>::edge_iterator e, e_begin, e_end;
        tie(e_begin, e_end) = edges(g);
        for (e = e_begin; e != e_end; ++e)
        {
            typename Hist::key_type key;
            key.first = _deg1(source(*e,g),g);
            key.second = _deg2(target(*e,g),g);
            hist[key]++;
            if(is_convertible<typename graph_traits<Graph>::directed_category, undirected_tag>::value)
            {
        	swap(key.first, key.second);
        	hist[key]++;
            }
        }
    }
    DegreeSelector1& _deg1;
    DegreeSelector2& _deg2;
};

//==============================================================================
// GetDegreeCorrelationHistogram(deg1, deg2)
// retrieves the degree correlation histogram                
//==============================================================================

template <class SecondDegreeSelectors>
struct choose_vertex_correlation_histogram
{
    choose_vertex_correlation_histogram(const GraphInterface &g, GraphInterface::deg_t deg1, 
        				GraphInterface::deg_t deg2, GraphInterface::hist2d_t &hist)
        : _g(g), _hist(hist) 
    {
        tie(_deg1, _deg_name1) = get_degree_type(deg1);
        tie(_deg2, _deg_name2) = get_degree_type(deg2);
    }

    template <class DegreeSelector1>
    struct check_second_degree
    {
        check_second_degree(choose_vertex_correlation_histogram<SecondDegreeSelectors> &parent):_parent(parent) {}
        template <class DegreeSelector2>
        void operator()(DegreeSelector2)
        {
            if (mpl::at<degree_selector_index,DegreeSelector2>::type::value == _parent._deg2)
            {
        	DegreeSelector1 deg1(_parent._deg_name1, _parent._g);
        	DegreeSelector2 deg2(_parent._deg_name2, _parent._g);
        	check_filter(_parent._g, bind<void>(get_correlation_histogram<DegreeSelector1,DegreeSelector2>(deg1,deg2), _1, var(_parent._hist)), 
        		     reverse_check(),directed_check());
            }
        }
        choose_vertex_correlation_histogram<SecondDegreeSelectors> _parent;
    };

    template <class DegreeSelector>
    void operator()(DegreeSelector)
    {
        if (mpl::at<degree_selector_index,DegreeSelector>::type::value == _deg1)
            mpl::for_each<SecondDegreeSelectors>(check_second_degree<DegreeSelector>(*this));
    }
    const GraphInterface &_g;
    GraphInterface::hist2d_t &_hist;
    GraphInterface::degree_t _deg1;
    string _deg_name1;
    GraphInterface::degree_t _deg2;
    string _deg_name2;
    
};

GraphInterface::hist2d_t 
GraphInterface::GetVertexCorrelationHistogram(GraphInterface::deg_t deg1, GraphInterface::deg_t deg2) const
{
    hist2d_t hist;
    typedef mpl::vector<in_degreeS, out_degreeS, total_degreeS, scalarS> degree_selectors;
    try 
    {
        mpl::for_each<degree_selectors>(choose_vertex_correlation_histogram<degree_selectors>(*this, deg1, deg2, hist));
    }
    catch (dynamic_get_failure &e)
    {
        throw GraphException("error getting scalar property: " + string(e.what()));
    }
    return hist;
}

//==============================================================================
// get_edge_correlation_histogram
// retrieves the generalized vertex-edge-vertex correlation histogram                
//==============================================================================

template <class DegreeSelector1, class DegreeSelector2>
struct get_edge_correlation_histogram
{
    get_edge_correlation_histogram(DegreeSelector1& deg1, DegreeSelector2& deg2, scalarS& scalar)
        : _edge_scalar(scalar), _deg1(deg1), _deg2(deg2) {}

    template <class Graph, class Hist>
    void operator()(Graph &g, Hist &hist) const
    {
        typename graph_traits<Graph>::edge_iterator e, e_begin, e_end;
        tie(e_begin, e_end) = edges(g);
        for (e = e_begin; e != e_end; ++e)
        {
            typename Hist::key_type key;
            get<0>(key) = _deg1(source(*e,g),g);
            get<1>(key) = _edge_scalar(*e, g);
            get<2>(key) = _deg2(target(*e,g),g);
            hist[key]++;
            if(is_convertible<typename graph_traits<Graph>::directed_category, undirected_tag>::value)
            {
        	swap(get<0>(key), get<2>(key));
        	hist[key]++;
            }
        }
    }
    scalarS& _edge_scalar;
    DegreeSelector1& _deg1;
    DegreeSelector2& _deg2;
};

//==============================================================================
// GetEdgeVertexCorrelationHistogram(deg1, scalar, deg2)
// retrieves the degree-edge-degree correlation histogram 
//==============================================================================

template <class SecondDegreeSelectors>
struct choose_edge_vertex_correlation_histogram
{
    choose_edge_vertex_correlation_histogram(const GraphInterface& g, GraphInterface::deg_t deg1,  scalarS& edge_scalar, 
        				     GraphInterface::deg_t deg2, GraphInterface::hist3d_t& hist)
        : _g(g), _edge_scalar(edge_scalar), _hist(hist) 
    {
        tie(_deg1, _deg_name1) = get_degree_type(deg1);
        tie(_deg2, _deg_name2) = get_degree_type(deg2);
    }

    template <class DegreeSelector1>
    struct check_second_degree
    {
        check_second_degree(choose_edge_vertex_correlation_histogram<SecondDegreeSelectors> &parent):_parent(parent) {}
        template <class DegreeSelector2>
        void operator()(DegreeSelector2)
        {
            if (mpl::at<degree_selector_index,DegreeSelector2>::type::value == _parent._deg2)
            {
        	DegreeSelector1 deg1(_parent._deg_name1, _parent._g);
        	DegreeSelector2 deg2(_parent._deg_name2, _parent._g);
        	check_filter(_parent._g, bind<void>(get_edge_correlation_histogram<DegreeSelector1,DegreeSelector2>(deg1,deg2,_parent._edge_scalar),
        					    _1, var(_parent._hist)), 
        		     reverse_check(),directed_check());
            }
        }
        choose_edge_vertex_correlation_histogram<SecondDegreeSelectors> _parent;
    };

    template <class DegreeSelector>
    void operator()(DegreeSelector)
    {
        if (mpl::at<degree_selector_index,DegreeSelector>::type::value == _deg1)
            mpl::for_each<SecondDegreeSelectors>(check_second_degree<DegreeSelector>(*this));
    }
    const GraphInterface& _g;
    scalarS& _edge_scalar;
    GraphInterface::hist3d_t& _hist;
    GraphInterface::degree_t _deg1;
    string _deg_name1;
    GraphInterface::degree_t _deg2;
    string _deg_name2;    
};

GraphInterface::hist3d_t 
GraphInterface::GetEdgeVertexCorrelationHistogram(GraphInterface::deg_t deg1, string edge_scalar, GraphInterface::deg_t deg2) const
{
    hist3d_t hist;

    scalarS scalar(edge_scalar, *this);
    try
    {
        typedef mpl::vector<in_degreeS, out_degreeS, total_degreeS, scalarS> degree_selectors;
        mpl::for_each<degree_selectors>(choose_edge_vertex_correlation_histogram<degree_selectors>(*this, deg1, scalar, deg2, hist));
    }
    catch (dynamic_get_failure& e)
    {
        throw GraphException("error getting scalar property: " + string(e.what()));
    }
    return hist;
}



