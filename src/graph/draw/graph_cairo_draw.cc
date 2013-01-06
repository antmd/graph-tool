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

#ifdef HAVE_CAIROMM

#include "graph_filtering.hh"

#include <boost/python.hpp>
#include <boost/utility/enable_if.hpp>

#include "graph_selectors.hh"
#include "graph_properties.hh"

#include <iostream>

#if (GCC_VERSION >= 40400)
#   include <tr1/unordered_map>
#else
#   include <boost/tr1/unordered_map.hpp>
#endif

#include <cairommconfig.h>
#include <cairomm/context.h>
#include <cairomm/surface.h>
#include <pycairo/pycairo.h>

#include <boost/mpl/map/map50.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;


enum vertex_attr_t {
    VERTEX_SHAPE = 100,
    VERTEX_COLOR,
    VERTEX_FILL_COLOR,
    VERTEX_SIZE,
    VERTEX_ASPECT,
    VERTEX_ANCHOR,
    VERTEX_PENWIDTH,
    VERTEX_HALO,
    VERTEX_HALO_COLOR,
    VERTEX_TEXT,
    VERTEX_TEXT_COLOR,
    VERTEX_TEXT_POSITION,
    VERTEX_FONT_FAMILY,
    VERTEX_FONT_SLANT,
    VERTEX_FONT_WEIGHT,
    VERTEX_FONT_SIZE,
    VERTEX_SURFACE,
    VERTEX_PIE_FRACTIONS,
    VERTEX_PIE_COLORS
};

enum edge_attr_t {
    EDGE_COLOR = 200,
    EDGE_PENWIDTH,
    EDGE_START_MARKER,
    EDGE_MID_MARKER,
    EDGE_END_MARKER,
    EDGE_MARKER_SIZE,
    EDGE_CONTROL_POINTS,
    EDGE_DASH_STYLE,
    EDGE_SLOPPY
};

enum vertex_shape_t {
    SHAPE_CIRCLE = 300,
    SHAPE_TRIANGLE,
    SHAPE_SQUARE,
    SHAPE_PENTAGON,
    SHAPE_HEXAGON,
    SHAPE_HEPTAGON,
    SHAPE_OCTAGON,
    SHAPE_DOUBLE_CIRCLE,
    SHAPE_DOUBLE_TRIANGLE,
    SHAPE_DOUBLE_SQUARE,
    SHAPE_DOUBLE_PENTAGON,
    SHAPE_DOUBLE_HEXAGON,
    SHAPE_DOUBLE_HEPTAGON,
    SHAPE_DOUBLE_OCTAGON,
    SHAPE_PIE
};

enum edge_marker_t {
    MARKER_SHAPE_NONE = 400,
    MARKER_SHAPE_ARROW,
    MARKER_SHAPE_CIRCLE,
    MARKER_SHAPE_SQUARE,
    MARKER_SHAPE_DIAMOND,
    MARKER_SHAPE_BAR
};

typedef pair<double, double> pos_t;
typedef tuple<double, double, double, double> color_t;
typedef tr1::unordered_map<int, boost::any> attrs_t;

typedef mpl::map28<
    mpl::pair<mpl::int_<VERTEX_SHAPE>, vertex_shape_t>,
    mpl::pair<mpl::int_<VERTEX_COLOR>, color_t>,
    mpl::pair<mpl::int_<VERTEX_FILL_COLOR>, color_t>,
    mpl::pair<mpl::int_<VERTEX_SIZE>, double>,
    mpl::pair<mpl::int_<VERTEX_ASPECT>, double>,
    mpl::pair<mpl::int_<VERTEX_ANCHOR>, int32_t>,
    mpl::pair<mpl::int_<VERTEX_PENWIDTH>, double>,
    mpl::pair<mpl::int_<VERTEX_HALO>, uint8_t>,
    mpl::pair<mpl::int_<VERTEX_HALO_COLOR>, color_t>,
    mpl::pair<mpl::int_<VERTEX_TEXT>, string>,
    mpl::pair<mpl::int_<VERTEX_TEXT_COLOR>, color_t>,
    mpl::pair<mpl::int_<VERTEX_TEXT_POSITION>, double>,
    mpl::pair<mpl::int_<VERTEX_FONT_FAMILY>, string>,
    mpl::pair<mpl::int_<VERTEX_FONT_SLANT>, int32_t>,
    mpl::pair<mpl::int_<VERTEX_FONT_WEIGHT>, int32_t>,
    mpl::pair<mpl::int_<VERTEX_FONT_SIZE>, double>,
    mpl::pair<mpl::int_<VERTEX_SURFACE>, python::object>,
    mpl::pair<mpl::int_<VERTEX_PIE_FRACTIONS>, vector<double> >,
    mpl::pair<mpl::int_<VERTEX_PIE_COLORS>, vector<color_t> >,
    mpl::pair<mpl::int_<EDGE_COLOR>, color_t>,
    mpl::pair<mpl::int_<EDGE_PENWIDTH>, double>,
    mpl::pair<mpl::int_<EDGE_START_MARKER>, edge_marker_t>,
    mpl::pair<mpl::int_<EDGE_MID_MARKER>, edge_marker_t>,
    mpl::pair<mpl::int_<EDGE_END_MARKER>, edge_marker_t>,
    mpl::pair<mpl::int_<EDGE_MARKER_SIZE>, double>,
    mpl::pair<mpl::int_<EDGE_CONTROL_POINTS>, vector<double> >,
    mpl::pair<mpl::int_<EDGE_DASH_STYLE>, vector<double> >,
    mpl::pair<mpl::int_<EDGE_SLOPPY>, uint8_t> >
        attr_types;

namespace std
{
ostream& operator<<(ostream& out, const color_t& c)
{
    out << get<0>(c) << " " << get<1>(c) << " " << get<2>(c) << " " << get<3>(c);
    return out;
}

istream& operator>>(istream& in, color_t& c)
{
    in >> get<0>(c) >> get<1>(c) >> get<2>(c) >>  get<3>(c);
    return in;
}
}

istream& operator>>(istream& in, vertex_shape_t& c)
{
    int tmp;
    in >> tmp;
    c = vertex_shape_t(tmp);
    return in;
}

istream& operator>>(istream& in, edge_marker_t& c)
{
    int tmp;
    in >> tmp;
    c = edge_marker_t(tmp);
    return in;
}

template <class T1, class T2>
struct specific_convert;

template <class Type1, class Type2>
struct Converter
{
    Type1 operator()(const Type2& v) const
    {
        return do_convert(v, is_convertible<Type2,Type1>());
    }

    Type1 do_convert(const Type2& v, mpl::bool_<true>) const
    {
        return Type1(v);
    }

    Type1 do_convert(const Type2& v, mpl::bool_<false>) const
    {
        return specific_convert<Type1,Type2>()(v);
    }

    template <class T1, class T2, class Enable = void>
    struct specific_convert
    {
        T1 operator()(const T2& v) const
        {
            return dispatch(v, typename is_convertible<T2, T1>::type());
        }

        T1 dispatch(const T2& v, mpl::true_) const
        {
            return T1(v);
        }

        T1 dispatch(const T2& v, mpl::false_) const
        {
            return lexical_cast<T1>(v);
        }
    };

    template <class T1> // noop
    struct specific_convert<T1, T1>
    {
        T1 operator()(const T1& v) const
        {
            return v;
        }
    };

    // specific specializations
    // string
    template <class T1>
    struct specific_convert<T1, string>
    {
        T1 operator()(const string& v) const
        {
            //uint8_t is not char, it is bool!
            if (is_same<T1, uint8_t>::value)
                return convert<T1,int>()(lexical_cast<int>(v));
            else
                return lexical_cast<T1>(v);
        }
    };

    // python::object
    template <class T1>
    struct specific_convert<T1,python::object>
    {
        T1 operator()(const python::object& v) const
        {
            python::extract<T1> x(v);
            if (x.check())
                return x();
            else
                throw bad_lexical_cast();
        }
    };

    template <class T2>
    struct specific_convert<string, T2,
                            typename enable_if
                                <typename mpl::not_<
                                     boost::is_same<T2,python::object> >::type>::type>
    {
        string operator()(const T2& v) const
        {
            //uint8_t is not char, it is bool!
            if (is_same<T2, uint8_t>::value)
                return convert<string, int>()(lexical_cast<int>(v));
            else
                return lexical_cast<string>(v);
        }
    };

    // vectors
    template <class T1, class T2>
    struct specific_convert<vector<T1>, vector<T2> >
    {
        vector<T1> operator()(const vector<T2>& v) const
        {
            vector<T1> v2(v.size());
            convert<T1,T2> c;
            for (size_t i = 0; i < v.size(); ++i)
                v2[i] = c(v[i]);
            return v2;
        }
    };

    // color
    template <class T2>
    struct specific_convert<color_t, vector<T2> >
    {
        specific_convert<double, T2> c;
        color_t operator()(const vector<T2>& cv) const
        {
            if (cv.size() < 4)
                throw bad_lexical_cast();
            return make_tuple(c(cv[0]), c(cv[1]), c(cv[2]), c(cv[3]));
        }
    };


    // vertex_shape_t
    template <class T2>
    struct specific_convert<vertex_shape_t, T2,
                            typename enable_if
                            <typename mpl::and_<typename mpl::not_<
                                                    boost::is_same<T2,string> >::type,
                                                typename mpl::not_<
                                                    boost::is_same<T2,python::object> >::type>::type>::type>
    {
        specific_convert<int, T2> c;
        vertex_shape_t operator()(const T2& v) const
        {
            return static_cast<vertex_shape_t>(c(v));
        }
    };

    // edge_marker_t
    template <class T2>
    struct specific_convert<edge_marker_t, T2,
                            typename enable_if
                            <typename mpl::and_<typename mpl::not_<
                                                    boost::is_same<T2,string> >::type,
                                                typename mpl::not_<
                                                    boost::is_same<T2,python::object> >::type>::type>::type>
    {
        specific_convert<int, T2> c;
        edge_marker_t operator()(const T2& v) const
        {
            return static_cast<edge_marker_t>(c(v));
        }
    };
};


template <class Descriptor>
class AttrDict
{
public:
    AttrDict(Descriptor descriptor, attrs_t& attrs, attrs_t& defaults)
        : _descriptor(descriptor), _attrs(attrs), _defaults(defaults) {}

    template <class Value>
    Value get(int k)
    {
        typeof(_attrs.begin()) iter = _attrs.find(k);
        if (iter != _attrs.end())
        {
            typedef DynamicPropertyMapWrap<Value, Descriptor, Converter> pmap_t;
            pmap_t pmap(any_cast<pmap_t>(iter->second));
            return pmap.get(_descriptor);
        }
        try
        {
            return any_cast<Value>(_defaults[k]);
        }
        catch (bad_any_cast&)
        {
            throw ValueException("Error getting attribute " + lexical_cast<string>(k) +
                                 ", wanted: " +
                                 python::detail::gcc_demangle(typeid(Value).name()) +
                                 ", got: " +
                                 python::detail::gcc_demangle(_defaults[k].type().name()));
        }
    }

private:
    Descriptor _descriptor;
    attrs_t& _attrs;
    attrs_t& _defaults;
};

void draw_polygon(size_t N, double radius, Cairo::Context& cr)
{
    cr.save();
    cr.rotate(M_PI * (1. / 2 - 1. / N));
    cr.move_to(radius, 0);
    for (size_t i = 0; i < N; ++i)
    {
        double angle = (2 * M_PI * (i + 1)) / N;
        cr.line_to(radius * cos(angle), radius * sin(angle));
    }
    cr.close_path();
    cr.restore();
}

double get_polygon_anchor(size_t N, double radius, double angle)
{
    double theta = angle - M_PI * (1. / 2 - 1. / N);
    if (N % 2 == 0)
        theta += M_PI / N;
    if (theta > 2 * M_PI)
        theta -= 2 * M_PI;
    if (theta < 2 * M_PI)
        theta += 2 * M_PI;
    theta = fmod(theta, 2. * M_PI / N);
    if (theta > M_PI / N)
        theta -= 2. * M_PI / N;
    return radius * cos(M_PI / N) / cos(theta);
}


void draw_pie(double radius, const vector<double>& f,
              const vector<color_t>& colors, Cairo::Context& cr)
{
    if (colors.empty())
        throw ValueException("No pie colors!");
    double s = 0;
    for (size_t i = 0; i < f.size(); ++i)
        s += f[i];
    double last = 0;
    double pos = 0;
    cr.save();
    cr.begin_new_path();
    for (size_t i = 0; i < f.size(); ++i)
    {
        pos += f[i];
        double angle = (2 * pos * M_PI) / s;
        cr.move_to(0, 0);
        cr.arc(0, 0, radius, last, angle);
        last = angle;
        size_t j = i % colors.size();
        cr.set_source_rgba(get<0>(colors[j]),
                           get<1>(colors[j]),
                           get<2>(colors[j]),
                           get<3>(colors[j]));
        cr.fill();
    }
    cr.restore();
}


void move_radially(pos_t& pos, const pos_t& origin, double dr)
{
    double angle = atan2(pos.second - origin.second,
                         pos.first - origin.first);
    if (angle < 0)
        angle += 2 * M_PI;
    pos.first += dr * cos(angle);
    pos.second += dr * sin(angle);
}

double get_user_dist(Cairo::Context& cr)
{
    double x = 1./sqrt(2), y = 1./sqrt(2);
    cr.device_to_user_distance(x, y);
    return sqrt(x * x + y * y);
}

void get_surface_size(Cairo::RefPtr<Cairo::Surface> sfc,
                      double& width, double& height)
{
    Cairo::RefPtr<Cairo::Context> cr = Cairo::Context::create(sfc);
    double x1, x2, y1, y2;

    cr->get_clip_extents(x1, y1, x2, y2);

    width = x2 - x1;
    height = y2 - y1;
}



template <class Descriptor>
class VertexShape
{
public:
    VertexShape(pos_t& pos, AttrDict<Descriptor> attrs)
        : _pos(pos), _attrs(attrs)

    {}

    double get_size(Cairo::Context& cr)
    {
        double size = _attrs.template get<double>(VERTEX_SIZE);
        size *= get_user_dist(cr);

        string text = _attrs.template get<string>(VERTEX_TEXT);
        if (text != "")
        {
            double text_pos = _attrs.template get<double>(VERTEX_TEXT_POSITION);
            if (text_pos == -1)
            {
                cr.select_font_face(_attrs.template get<string>(VERTEX_FONT_FAMILY),
                                    static_cast<Cairo::FontSlant>(_attrs.template get<int32_t>(VERTEX_FONT_SLANT)),
                                    static_cast<Cairo::FontWeight>(_attrs.template get<int32_t>(VERTEX_FONT_WEIGHT)));
                cr.set_font_size(_attrs.template get<double>(VERTEX_FONT_SIZE) *
                                 get_user_dist(cr));
                Cairo::TextExtents extents;
                cr.get_text_extents(text, extents);
                double s = max(extents.width, extents.height) * 1.4;
                vertex_shape_t shape = _attrs.template get<vertex_shape_t>(VERTEX_SHAPE);
                if (shape >= SHAPE_DOUBLE_CIRCLE)
                {
                    s /= 0.7;
                    double pw = _attrs.template get<double>(VERTEX_PENWIDTH);
                    pw *= get_user_dist(cr);
                    s += pw;
                }
                size = max(size, s);
            }
        }
        return size;
    }

    pos_t get_anchor(const pos_t& origin, Cairo::Context& cr)
    {
        int anchor_type =_attrs.template get<int32_t>(VERTEX_ANCHOR);
        if (anchor_type == 0)
            return _pos;

        double angle = atan2(_pos.second - origin.second,
                             _pos.first - origin.first);
        if (angle < 0)
            angle += 2 * M_PI;
        double r = get_size(cr) / 2;
        double dr = r;

        double pw = _attrs.template get<double>(VERTEX_PENWIDTH);
        pw *= get_user_dist(cr);
        r += pw / 2.5;

        pos_t anchor;
        size_t nsides = 0;
        vertex_shape_t shape = _attrs.template get<vertex_shape_t>(VERTEX_SHAPE);
        switch (shape)
        {
        case SHAPE_TRIANGLE:
        case SHAPE_SQUARE:
        case SHAPE_PENTAGON:
        case SHAPE_HEXAGON:
        case SHAPE_HEPTAGON:
        case SHAPE_OCTAGON:
        case SHAPE_DOUBLE_TRIANGLE:
        case SHAPE_DOUBLE_SQUARE:
        case SHAPE_DOUBLE_PENTAGON:
        case SHAPE_DOUBLE_HEXAGON:
        case SHAPE_DOUBLE_HEPTAGON:
        case SHAPE_DOUBLE_OCTAGON:
            nsides = shape - SHAPE_TRIANGLE + 3;
            if (nsides > 8)
                nsides -= 7;
            dr = get_polygon_anchor(nsides, r, angle);
            break;
        case SHAPE_CIRCLE:
        case SHAPE_DOUBLE_CIRCLE:
        case SHAPE_PIE:
            dr = r;
            break;
        default:
            throw ValueException("Invalid vertex shape: " +
                                 lexical_cast<string>(int(_attrs.template get<vertex_shape_t>(VERTEX_SHAPE))));
        }

        double aspect = _attrs.template get<double>(VERTEX_ASPECT);

        anchor = _pos;
        anchor.first -= dr * cos(angle) * aspect;
        anchor.second -= dr * sin(angle);

        return anchor;
    }

    pos_t get_pos()
    {
        return _pos;
    }

    void draw(Cairo::Context& cr, bool outline=false)
    {
        color_t color, fillcolor;
        double size, pw, aspect;
        size = get_size(cr);
        aspect = _attrs.template get<double>(VERTEX_ASPECT);

        string text = _attrs.template get<string>(VERTEX_TEXT);
        double text_pos = 0;
        if (text != "" && !outline)
        {
            cr.select_font_face(_attrs.template get<string>(VERTEX_FONT_FAMILY),
                                static_cast<Cairo::FontSlant>(_attrs.template get<int32_t>(VERTEX_FONT_SLANT)),
                                static_cast<Cairo::FontWeight>(_attrs.template get<int32_t>(VERTEX_FONT_WEIGHT)));
            cr.set_font_size(_attrs.template get<double>(VERTEX_FONT_SIZE) *
                             get_user_dist(cr));
            text_pos = _attrs.template get<double>(VERTEX_TEXT_POSITION);
        }

        if (!outline)
            cr.save();
        cr.translate(_pos.first, _pos.second);

        if (_attrs.template get<uint8_t>(VERTEX_HALO) && !outline)
        {
            color_t c = _attrs.template get<color_t>(VERTEX_HALO_COLOR);
            cr.set_source_rgba(get<0>(c), get<1>(c), get<2>(c), get<3>(c));
            cr.save();
            cr.scale(aspect, 1.0);
            cr.arc(0, 0, size, 0, 2 * M_PI);
            cr.restore();
            cr.fill();
        }

        python::object osrc = _attrs.template get<python::object>(VERTEX_SURFACE);
        if (osrc == python::object())
        {
            pw =_attrs.template get<double>(VERTEX_PENWIDTH);
            pw *= get_user_dist(cr);
            cr.set_line_width(pw);

            color = _attrs.template get<color_t>(VERTEX_COLOR);
            cr.set_source_rgba(get<0>(color), get<1>(color), get<2>(color),
                               get<3>(color));

            size_t nsides = 0;
            vertex_shape_t shape = _attrs.template get<vertex_shape_t>(VERTEX_SHAPE);
            switch (shape)
            {
            case SHAPE_CIRCLE:
            case SHAPE_DOUBLE_CIRCLE:
                cr.save();
                cr.scale(aspect, 1.0);
                cr.arc(0, 0, size / 2., 0, 2 * M_PI);
                cr.close_path();
                cr.restore();
                if (shape == SHAPE_DOUBLE_CIRCLE && !outline)
                {
                    cr.stroke();
                    cr.save();
                    cr.scale(aspect, 1.0);
                    cr.arc(0, 0, min(size / 2 - 2 * pw,
                                     size * 0.8 / 2),
                           0, 2 * M_PI);
                    cr.restore();
                }
                break;
            case SHAPE_PIE:
                {
                    if (!outline)
                    {
                        vector<double> f = _attrs.template get<vector<double> >(VERTEX_PIE_FRACTIONS);
                        vector<color_t> fcolors = _attrs.template get<vector<color_t> >(VERTEX_PIE_COLORS);
                        draw_pie(size / 2 + pw / 2, f, fcolors, cr);
                    }
                    else
                    {
                        cr.save();
                        cr.arc(0, 0, size / 2., 0, 2 * M_PI);
                        cr.close_path();
                        cr.restore();
                    }
                }
                break;
            case SHAPE_TRIANGLE:
            case SHAPE_SQUARE:
            case SHAPE_PENTAGON:
            case SHAPE_HEXAGON:
            case SHAPE_HEPTAGON:
            case SHAPE_OCTAGON:
            case SHAPE_DOUBLE_TRIANGLE:
            case SHAPE_DOUBLE_SQUARE:
            case SHAPE_DOUBLE_PENTAGON:
            case SHAPE_DOUBLE_HEXAGON:
            case SHAPE_DOUBLE_HEPTAGON:
            case SHAPE_DOUBLE_OCTAGON:
                nsides = shape - SHAPE_TRIANGLE + 3;
                if (nsides > 8)
                    nsides -= 7;
                cr.save();
                cr.scale(aspect, 1.0);
                draw_polygon(nsides, size / 2, cr);
                cr.restore();
                if (shape >= SHAPE_DOUBLE_TRIANGLE && !outline)
                {
                    cr.stroke();
                    cr.save();
                    cr.scale(aspect, 1.0);
                    draw_polygon(nsides, min(size / 2 - 2 * pw,
                                             size * 0.8 / 2), cr);
                    cr.restore();
                }
                break;
            default:
                throw ValueException("Invalid vertex shape: " +
                                     lexical_cast<string>(int(_attrs.template get<vertex_shape_t>(VERTEX_SHAPE))));
            }

            if (!outline && shape != SHAPE_PIE)
            {
                fillcolor = _attrs.template get<color_t>(VERTEX_FILL_COLOR);
                cr.set_source_rgba(get<0>(fillcolor), get<1>(fillcolor),
                                   get<2>(fillcolor), get<3>(fillcolor));
                cr.fill_preserve();

                cr.set_source_rgba(get<0>(color), get<1>(color), get<2>(color),
                                       get<3>(color));
                cr.stroke();
            }
        }
        else
        {
            if (!outline)
            {
                double swidth, sheight;
                PycairoSurface* src = (PycairoSurface*) osrc.ptr();
                Cairo::RefPtr<Cairo::Surface> surface(new Cairo::Surface(src->surface));
                get_surface_size(surface, swidth, sheight);
                Cairo::RefPtr<Cairo::SurfacePattern> pat(Cairo::SurfacePattern::create(surface));
                //pat->set_extend(Cairo::EXTEND_REPEAT);

                double r = size / sqrt(2);
                double scale = r / max(swidth / aspect, sheight);

                Cairo::Matrix m = Cairo::identity_matrix();
                m.translate(swidth / 2, sheight / 2);
                m.scale(1. / scale, 1. / scale);
                pat->set_matrix(m);

                cr.set_source(pat);
                cr.rectangle(-r * aspect / 2, -r / 2, r * aspect, r);
                cr.fill();
            }
        }

        if (text != "" && !outline)
        {
            cr.save();
            Cairo::TextExtents extents;
            cr.get_text_extents(text, extents);

            if (text_pos < 0)
            {
                cr.translate(-extents.width/2 - extents.x_bearing,
                             extents.height / 2);
            }
            else
            {
                pos_t origin;
                origin.first = _pos.first + size * cos(text_pos);
                origin.second = _pos.second + size * sin(text_pos);
                pos_t anchor = get_anchor(origin, cr);
                double angle = atan2(_pos.second - anchor.second,
                                     _pos.first - anchor.first);
                anchor.first = size * 1.2 * cos(angle);
                anchor.second = size * 1.2 * sin(angle);
                if (anchor.first < 0)
                    anchor.first -= extents.width;
                if (anchor.second > 0)
                    anchor.second += extents.height;
                cr.translate(anchor.first, anchor.second);
            }
            color = _attrs.template get<color_t>(VERTEX_TEXT_COLOR);
            cr.set_source_rgba(get<0>(color), get<1>(color), get<2>(color),
                               get<3>(color));
            cr.show_text(text);
            cr.begin_new_path();
            cr.restore();
        }

        if (!outline)
            cr.restore();
        else
            cr.translate(-_pos.first, -_pos.second);

    }

    template <class, class>
    friend class EdgeShape;

private:
    pos_t _pos;
protected:
    AttrDict<Descriptor> _attrs;
};

template <class Descriptor, class VertexShape>
class EdgeShape
{
public:
    EdgeShape(VertexShape& s, VertexShape& t, AttrDict<Descriptor> attrs)
        : _s(s), _t(t), _attrs(attrs) {}

    void draw(Cairo::Context& cr)
    {
        pos_t pos_begin, pos_end;

        vector<double> controls =
            _attrs.template get<vector<double> >(EDGE_CONTROL_POINTS);

        edge_marker_t start_marker = _attrs.template get<edge_marker_t>(EDGE_START_MARKER);
        edge_marker_t end_marker = _attrs.template get<edge_marker_t>(EDGE_END_MARKER);
        bool sloppy = _attrs.template get<uint8_t>(EDGE_SLOPPY);

        pos_begin = _s.get_anchor(_t.get_pos(), cr);
        pos_end = _t.get_anchor(_s.get_pos(), cr);

        cr.save();

        if (controls.size() >= 4)
        {
            double angle = 0;
            double len = 0;
            if (pos_begin != pos_end)
            {
                angle = atan2(pos_end.second - pos_begin.second,
                              pos_end.first - pos_begin.first);
                len = sqrt(pow(pos_end.first - pos_begin.first, 2) +
                           pow(pos_end.second - pos_begin.second, 2));

                cr.save();
                cr.translate(pos_begin.first, pos_begin.second);
                cr.rotate(angle);
                cr.scale(len, 1.);

                for (size_t i = 0; i < controls.size(); i += 2)
                    cr.user_to_device(controls[i], controls[i + 1]);
                cr.restore();

                for (size_t i = 0; i < controls.size(); i += 2)
                    cr.device_to_user(controls[i], controls[i + 1]);
            }
            else
            {
                pos_begin = _s.get_pos();
                len = M_PI * _s.get_size(cr);
                cr.save();
                cr.translate(pos_begin.first, pos_begin.second);
                cr.scale(len / sqrt(2), len / sqrt(2));
                for (size_t i = 0; i < controls.size(); i += 2)
                    cr.user_to_device(controls[i], controls[i + 1]);
                cr.restore();

                for (size_t i = 0; i < controls.size(); i += 2)
                    cr.device_to_user(controls[i], controls[i + 1]);
            }

            size_t N = controls.size();
            pos_begin = _s.get_anchor(make_pair(controls[0],
                                                controls[1]), cr);
            pos_end = _t.get_anchor(make_pair(controls[N - 2],
                                              controls[N - 1]), cr);
        }


        pos_t pos_end_c = pos_end, pos_begin_c = pos_begin;
        double marker_size = _attrs.template get<double>(EDGE_MARKER_SIZE);
        marker_size *= get_user_dist(cr);

        if (start_marker == MARKER_SHAPE_NONE && !sloppy)
            pos_begin = _s.get_pos();
        else if ((start_marker != MARKER_SHAPE_NONE && start_marker != MARKER_SHAPE_BAR))
            move_radially(pos_begin_c, _s.get_pos(), marker_size / 2);
        if (end_marker == MARKER_SHAPE_NONE && !sloppy)
            pos_end = _t.get_pos();
        else if ((end_marker != MARKER_SHAPE_NONE && end_marker != MARKER_SHAPE_BAR))
            move_radially(pos_end_c, _t.get_pos(), marker_size / 2);

        color_t color = _attrs.template get<color_t>(EDGE_COLOR);
        double pw;
        pw = _attrs.template get<double>(EDGE_PENWIDTH);
        pw *= get_user_dist(cr);
        cr.set_line_width(pw);

        double a = get<3>(color);
        a *= get<3>(_s._attrs.template get<color_t>(VERTEX_COLOR));
        a *= get<3>(_s._attrs.template get<color_t>(VERTEX_FILL_COLOR));
        a *= get<3>(_t._attrs.template get<color_t>(VERTEX_COLOR));
        a *= get<3>(_t._attrs.template get<color_t>(VERTEX_FILL_COLOR));

        if (!sloppy && a < 1)
        {
            // set the clip region to the correct size for better push/pop_group
            // performance
            draw_edge_markers(pos_begin, pos_end, controls, marker_size, cr);
            draw_edge_line(pos_begin_c, pos_end_c, controls, cr);
            double sx1, sy1, sx2, sy2;
            cr.get_stroke_extents(sx1, sy1, sx2, sy2);
            cr.begin_new_path();
            cr.rectangle(sx1, sy1, sx2 - sx1, sy2 - sy1);
            _s.draw(cr, true);
            _t.draw(cr, true);
            cr.set_fill_rule(Cairo::FILL_RULE_EVEN_ODD);
            cr.clip();

            // seamlessly blend in separate surface
            cr.push_group();
            draw_edge_markers(pos_begin, pos_end, controls, marker_size, cr);
            cr.set_source_rgba(get<0>(color), get<1>(color), get<2>(color), 1);
            cr.fill();
            draw_edge_line(pos_begin_c, pos_end_c, controls, cr);
            cr.set_line_width(pw);
            cr.stroke();
            vector<double> empty;
            cr.set_dash(empty, 0);
            cr.pop_group_to_source();
            cr.reset_clip();
            cr.paint_with_alpha(get<3>(color));
        }
        else
        {
            cr.set_source_rgba(get<0>(color), get<1>(color), get<2>(color), get<3>(color));
            draw_edge_markers(pos_begin, pos_end, controls, marker_size, cr);
            cr.fill();
            draw_edge_line(pos_begin_c, pos_end_c, controls, cr);
            cr.stroke();
        }

        cr.restore();
    }

    void draw_edge_markers(pos_t& pos_begin, pos_t& pos_end,
                           vector<double>& controls,
                           double marker_size,
                           Cairo::Context& cr)
    {
        double len = sqrt(pow(pos_end.first - pos_begin.first, 2) +
                          pow(pos_end.second - pos_begin.second, 2));
        double angle_b, angle_e;
        if (controls.size() >= 4)
        {
            size_t N = controls.size();
            angle_b = atan2(controls[1] - pos_begin.second,
                            controls[0] - pos_begin.first);
            angle_e = atan2(pos_end.second - controls[N - 1],
                            pos_end.first - controls[N - 2]);
        }
        else
        {
            angle_b = angle_e = atan2(pos_end.second - pos_begin.second,
                                      pos_end.first - pos_begin.first);
        }

        cr.save();
        cr.translate(pos_end.first, pos_end.second);
        cr.rotate(angle_e);
        draw_marker(EDGE_END_MARKER, marker_size, cr);
        cr.restore();

        cr.save();
        cr.translate(pos_begin.first, pos_begin.second);
        cr.rotate(angle_b);
        cr.translate(marker_size, 0);
        draw_marker(EDGE_START_MARKER, marker_size, cr);
        cr.restore();

        if (angle_e == angle_b)
        {
            cr.save();
            cr.translate(pos_end.first, pos_end.second);
            cr.rotate(angle_b);
            cr.translate(-len/2. + marker_size / 2, 0);
            draw_marker(EDGE_MID_MARKER, marker_size, cr);
            cr.restore();
        }
    }

    void draw_edge_line(pos_t& pos_begin_c, pos_t& pos_end_c,
                        vector<double>& controls, Cairo::Context& cr)
    {
        cr.move_to(pos_begin_c.first, pos_begin_c.second);
        vector<double> dashes = _attrs.template get<vector<double> >(EDGE_DASH_STYLE);

        if (dashes.size() > 2)
        {
            double offset = dashes.back();
            dashes.pop_back();
            cr.set_dash(dashes, offset);
        }

        if (controls.size() >= 4)
        {
            size_t step = (controls.size() > 4) ? 6 : 4;
            for (size_t i = 0; i < controls.size(); i += step)
                if (i + 5 >= controls.size())
                    cr.curve_to(controls[i], controls[i + 1],
                                controls[i + 2], controls[i + 3],
                                pos_end_c.first, pos_end_c.second);
                else
                    cr.curve_to(controls[i], controls[i + 1],
                                controls[i + 2], controls[i + 3],
                                controls[i + 4], controls[i + 5]);
        }
        else
        {
            cr.line_to(pos_end_c.first, pos_end_c.second);
        }
    }

    void draw_marker(edge_attr_t attr, double size, Cairo::Context& cr)
    {
        edge_marker_t marker = _attrs.template get<edge_marker_t>(attr);
        switch (marker)
        {
        case MARKER_SHAPE_ARROW:
            {
                double angle = M_PI / 7.;
                double y = tan(angle) * size;
                double x = -size * 0.6;

                cr.move_to(0, 0);
                cr.line_to(-size, y);
                cr.line_to(x, 0);
                cr.line_to(-size, -y);
                cr.line_to(0, 0);
                cr.close_path();
            }
            break;
        case MARKER_SHAPE_CIRCLE:
            cr.arc(-size / 2, 0, size / 2., 0, 2 * M_PI);
            break;
        case MARKER_SHAPE_SQUARE:
            cr.save();
            cr.translate(-size / 2, 0);
            draw_polygon(4, size / 2, cr);
            cr.restore();
            break;
        case MARKER_SHAPE_DIAMOND:
            cr.save();
            cr.translate(-size / 2, 0);
            cr.rotate(M_PI / 4.);
            cr.scale(sqrt(2), 1);
            draw_polygon(4, size / 2, cr);
            cr.restore();
            break;
        case MARKER_SHAPE_BAR:
            {
                double x = 0, w = size / 4;
                if (attr == EDGE_START_MARKER)
                    x = -size + w;
                cr.move_to(x, 0);
                cr.line_to(x, -size / 2);
                cr.line_to(x - w, -size / 2);
                cr.line_to(x - w, size / 2);
                cr.line_to(x, size / 2);
                cr.close_path();
            }
            break;
        case MARKER_SHAPE_NONE:
            break;
        default:
            throw ValueException("Invalid edge marker: " +
                                 lexical_cast<string>(int(marker)));
        }
    }

private:
    VertexShape _s;
    VertexShape _t;
    AttrDict<Descriptor> _attrs;
};

template <class Graph, class VertexIterator, class PosMap>
void draw_vertices(Graph& g, pair<VertexIterator,VertexIterator> v_range,
                   PosMap pos_map, attrs_t& attrs, attrs_t& defaults,
                   Cairo::Context& cr)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    for(VertexIterator v = v_range.first; v != v_range.second; ++v)
    {
        pos_t pos;
        pos.first = pos_map[*v][0];
        pos.second = pos_map[*v][1];
        VertexShape<vertex_t> vs(pos, AttrDict<vertex_t>(*v, attrs, defaults));
        vs.draw(cr);
    }
}

template <class Graph, class EdgeIterator, class PosMap>
void draw_edges(Graph& g, pair<EdgeIterator, EdgeIterator> e_range,
                PosMap pos_map, attrs_t& eattrs, attrs_t& edefaults,
                attrs_t& vattrs, attrs_t& vdefaults, Cairo::Context& cr)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    typedef typename graph_traits<Graph>::edge_descriptor edge_t;
    for(EdgeIterator e = e_range.first; e != e_range.second; ++e)
    {
        vertex_t s, t;
        s = source(*e, g);
        t = target(*e, g);

        pos_t spos, tpos;
        spos.first = pos_map[s][0];
        spos.second = pos_map[s][1];
        tpos.first = pos_map[t][0];
        tpos.second = pos_map[t][1];
        VertexShape<vertex_t> ss(spos, AttrDict<vertex_t>(s, vattrs, vdefaults));
        VertexShape<vertex_t> ts(tpos, AttrDict<vertex_t>(t, vattrs, vdefaults));

        EdgeShape<edge_t,VertexShape<vertex_t> > es(ss, ts,
                                                    AttrDict<edge_t>(*e, eattrs,
                                                                     edefaults));
        es.draw(cr);
    }
}

struct no_order {};

template <class Iterator>
struct ordered_range
{
    typedef typename iterator_traits<Iterator>::value_type val_t;

    ordered_range(const pair<Iterator, Iterator>& range)
        : _range(range)
    {
    }

    template <class Order>
    struct val_cmp
    {
        val_cmp(Order order): _order(order) {}
        bool operator()(const val_t& a, const val_t& b)
        {
            return get(_order, a) < get(_order, b);
        }
        Order _order;
    };

    template <class Order>
    pair<typename vector<val_t>::iterator, typename vector<val_t>::iterator>
    get_range(Order order)
    {
        if (_ordered.empty())
        {
            for (Iterator iter = _range.first; iter != _range.second; ++iter)
                _ordered.push_back(*iter);
            sort(_ordered.begin(), _ordered.end(), val_cmp<Order>(order));
        }
        return make_pair(_ordered.begin(), _ordered.end());
    }

    pair<Iterator, Iterator> get_range(no_order)
    {
        return _range;
    }

    pair<Iterator, Iterator> _range;
    vector<val_t> _ordered;
};


struct do_cairo_draw_edges
{
    template <class Graph, class PosMap, class EdgeOrder>
    void operator()(Graph& g, PosMap pos, EdgeOrder edge_order,
                    attrs_t& vattrs, attrs_t& eattrs, attrs_t& vdefaults,
                    attrs_t& edefaults, Cairo::Context& cr) const
    {
        ordered_range<typename graph_traits<Graph>::edge_iterator>
            edge_range(edges(g));
        draw_edges(g, edge_range.get_range(edge_order), pos, eattrs,
                   edefaults, vattrs, vdefaults, cr);
    }
};

struct do_cairo_draw_vertices
{
    template <class Graph, class PosMap, class VertexOrder>
    void operator()(Graph& g, PosMap pos, VertexOrder vertex_order,
                    attrs_t& vattrs, attrs_t& eattrs, attrs_t& vdefaults,
                    attrs_t& edefaults, Cairo::Context& cr) const
    {
        ordered_range<typename graph_traits<Graph>::vertex_iterator>
            vertex_range(vertices(g));
        draw_vertices(g, vertex_range.get_range(vertex_order), pos, vattrs,
                      vdefaults, cr);
    }
};

template <class Descriptor, class PropMaps>
struct get_pmap
{
    get_pmap(boost::any& opmap, boost::any& pmap, int type)
        : _opmap(opmap), _pmap(pmap), _type(type) {}
    boost::any& _opmap;
    boost::any& _pmap;
    int _type;

    template <class ValueType>
    void operator()(ValueType) const
    {
        typedef typename ValueType::second val_t;

        typedef DynamicPropertyMapWrap<val_t, Descriptor, Converter> pmap_t;

        if (_type == ValueType::first::value)
            _pmap = pmap_t(_opmap, PropMaps());
    }
};


template <class Descriptor, class PropMaps>
void populate_attrs(python::dict vattrs, attrs_t& attrs)
{
    python::list items = vattrs.items();
    for (int i = 0; i < python::len(items); ++i)
    {
        boost::any oattr = python::extract<boost::any>(items[i][1]);
        boost::any pmap;
        int type = python::extract<int>(items[i][0]);
        mpl::for_each<attr_types>(get_pmap<Descriptor,PropMaps>(oattr, pmap,
                                                                type));
        attrs[type] = pmap;
    }
}


struct get_dval
{
    get_dval(python::object& odval, boost::any& dval, int type)
        : _odval(odval), _dval(dval), _type(type) {}
    python::object& _odval;
    boost::any& _dval;
    int _type;

    template <class ValueType>
    void operator()(ValueType) const
    {
        typedef typename ValueType::second val_t;
        if (_type == int(ValueType::first::value))
            _dval = python::extract<val_t>(_odval)();
    }
};

void populate_defaults(python::dict odefaults, attrs_t& defaults)
{
    python::list items = odefaults.items();
    for (int i = 0; i < python::len(items); ++i)
    {
        python::object odval = items[i][1];
        boost::any dval;
        int type = python::extract<int>(items[i][0]);
        mpl::for_each<attr_types>(get_dval(odval, dval, type));
        if (dval.empty())
            throw ValueException("Invalid attribute type.");
        defaults[type] = dval;
    }
}

struct populate_edge_attrs
{
    template <class Graph>
    void operator()(Graph& g, python::dict oeattrs, attrs_t& eattrs,
                    python::dict oedefaults, attrs_t& edefaults) const
    {
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;
        populate_attrs<edge_t, edge_properties>(oeattrs, eattrs);
        populate_defaults(oedefaults, edefaults);
    }
};


void cairo_draw(GraphInterface& gi,
                boost::any pos,
                boost::any vorder,
                boost::any eorder,
                bool nodesfirst,
                python::dict ovattrs,
                python::dict oeattrs,
                python::dict ovdefaults,
                python::dict oedefaults,
                python::object ocr)
{
    attrs_t vattrs, eattrs, vdefaults, edefaults;
    typedef graph_traits<GraphInterface::multigraph_t>::vertex_descriptor vertex_t;

    populate_attrs<vertex_t, vertex_properties>(ovattrs, vattrs);
    populate_defaults(ovdefaults, vdefaults);

    run_action<graph_tool::detail::always_directed>()
        (gi, bind<void>(populate_edge_attrs(), _1, oeattrs,
                        ref(eattrs), oedefaults, ref(edefaults)))();

    typedef mpl::push_back<vertex_scalar_properties, no_order>::type
        vorder_t;
    typedef mpl::push_back<edge_scalar_properties, no_order>::type
        eorder_t;
    if (vorder.empty())
        vorder = no_order();
    if (eorder.empty())
        eorder = no_order();

    Cairo::Context cr(PycairoContext_GET(ocr.ptr()));
    if (nodesfirst)
        run_action<graph_tool::detail::always_directed>()
            (gi, bind<void>(do_cairo_draw_vertices(), _1, _2, _3,
                            ref(vattrs), ref(eattrs), ref(vdefaults),
                            ref(edefaults), ref(cr)),
             vertex_scalar_vector_properties(),
             vorder_t())(pos, vorder);
    run_action<graph_tool::detail::always_directed>()
        (gi, bind<void>(do_cairo_draw_edges(), _1, _2, _3,
                        ref(vattrs), ref(eattrs), ref(vdefaults),
                        ref(edefaults), ref(cr)),
         vertex_scalar_vector_properties(),
         eorder_t())(pos, eorder);
    if (!nodesfirst)
        run_action<graph_tool::detail::always_directed>()
            (gi, bind<void>(do_cairo_draw_vertices(), _1, _2, _3,
                            ref(vattrs), ref(eattrs), ref(vdefaults),
                            ref(edefaults), ref(cr)),
             vertex_scalar_vector_properties(),
             vorder_t())(pos, vorder);
}

struct do_apply_transforms
{
    template <class Graph, class PosMap>
    void operator()(Graph& g, PosMap pos, Cairo::Matrix& m) const
    {
        typename graph_traits<Graph>::vertex_iterator v, v_end;
        for(tie(v, v_end) = vertices(g); v != v_end; ++v)
        {
            pos[*v].resize(2);
            double x = pos[*v][0], y = pos[*v][1];
            m.transform_point(x, y);
            pos[*v][0] = x;
            pos[*v][1] = y;
        }
    }
};

void apply_transforms(GraphInterface& gi, boost::any pos, double xx, double yx,
                      double xy, double yy, double x0, double y0)
{
    Cairo::Matrix m(xx, yx, xy, yy, x0, y0);
    run_action<graph_tool::detail::always_directed>()
        (gi, bind<void>(do_apply_transforms(), _1, _2, ref(m)),
         vertex_scalar_vector_properties())(pos);
}

struct do_put_parallel_splines
{
    template <class Graph, class PosMap, class LabelMap, class SplinesMap>
    void operator()(Graph& g, PosMap pos, LabelMap l, SplinesMap spline,
                    double loop_angle, double parallel_distance) const
    {
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        typedef typename property_traits<SplinesMap>::key_type skey_t;

        pair<double, double> cm;
        cm.first = cm.second = 0;
        size_t n = 0;
        typename graph_traits<Graph>::vertex_iterator v, v_end;
        for(tie(v, v_end) = vertices(g); v != v_end; ++v)
        {
            cm.first += get(pos, *v)[0];
            cm.second += get(pos, *v)[1];
            ++n;
        }
        cm.first /= n;
        cm.second /= n;

        typename graph_traits<Graph>::edge_iterator e, e_end;
        for(tie(e, e_end) = edges(g); e != e_end; ++e)
        {
            if (target(*e, g) == source(*e, g))
            {
                vertex_t s = source(*e, g);
                typename graph_traits<Graph>::out_edge_iterator eo, eo_end;
                vector<edge_t> pes;
                for (tie(eo, eo_end) = out_edges(s, g); eo != eo_end; ++eo)
                {
                    if (target(*eo, g) == source(*eo, g))
                        pes.push_back(*eo);
                }

                pair<double, double> dist;
                dist.first = get(pos, s)[0] - cm.first;
                dist.second = get(pos, s)[1] - cm.second;
                double theta = atan2(dist.second, dist.first) - M_PI / 2;
                if (!isnan(loop_angle))
                    theta = loop_angle;
                typename property_traits<SplinesMap>::value_type sp(4), sp2(4);
                for (size_t j = 0; j < pes.size(); ++j)
                {
                    sp[0] = -0.5 * (j + 1);
                    sp[1] = 0.75 * (j + 1);
                    sp[2] = 0.5 * (j + 1);
                    sp[3] = 0.75 * (j + 1);

                    sp2[0] = sp[0] * cos(theta) - sp[1] * sin(theta);
                    sp2[1] = sp[0] * sin(theta) + sp[1] * cos(theta);
                    sp2[2] = sp[2] * cos(theta) - sp[3] * sin(theta);
                    sp2[3] = sp[2] * sin(theta) + sp[3] * cos(theta);

                    put(spline, skey_t(pes[j]), sp2);
                }
            }
            else if (get(l, *e) == 1)
            {
                vector<pair<edge_t, bool> > pes;
                vertex_t s = source(*e, g);
                typename graph_traits<Graph>::out_edge_iterator eo, eo_end;
                for (tie(eo, eo_end) = out_edges(s, g); eo != eo_end; ++eo)
                {
                    if (target(*eo, g) == target(*e, g))
                        pes.push_back(make_pair(*eo, true));
                }

                typename graph_traits<Graph>::in_edge_iterator ei, ei_end;
                for (tie(ei, ei_end) = in_edges(s, g); ei != ei_end; ++ei)
                {
                    if (source(*ei, g) == target(*e, g))
                        pes.push_back(make_pair(*ei, false));
                }

                typename property_traits<SplinesMap>::value_type sp(4, 0);
                double n = (pes.size() - 1.) / 2.;
                for (size_t j = 0; j < pes.size(); ++j)
                {
                    typedef typename property_traits<SplinesMap>::value_type::value_type val_t;
                    double one = pes[j].second ? 1 : -1;
                    sp[0] = val_t(0.3);
                    sp[1] = val_t(one * (j - n) * parallel_distance / n);
                    sp[2] = val_t(0.7);
                    sp[3] = val_t(one * (j - n) * parallel_distance / n);
                    put(spline, skey_t(pes[j].first), sp);
                }
            }
        }
    }
};

void put_parallel_splines(GraphInterface& gi, boost::any opos,
                          boost::any ol, boost::any splines,
                          double loop_angle, double parallel_distance)
{
    DynamicPropertyMapWrap<vector<double>, GraphInterface::vertex_t>
        pos(opos, vertex_scalar_vector_properties());
    DynamicPropertyMapWrap<int, GraphInterface::edge_t>
        l(ol, edge_scalar_properties());

    run_action<graph_tool::detail::always_directed>()
        (gi, bind<void>(do_put_parallel_splines(), _1, pos, l, _2, loop_angle,
                        parallel_distance),
         edge_scalar_vector_properties())(splines);
}


using namespace boost::python;

struct color_from_list
{
    color_from_list()
    {
        converter::registry::push_back
            (&convertible, &construct,
             boost::python::type_id<color_t>());
    }

    static void* convertible(PyObject* obj_ptr)
    {
        handle<> x(borrowed(obj_ptr));
        object o(x);
        size_t N = len(o);
        if (N < 4)
            return 0;
        return obj_ptr;
    }

    static void construct(PyObject* obj_ptr,
                          converter::rvalue_from_python_stage1_data* data)
    {
        handle<> x(borrowed(obj_ptr));
        object o(x);
        color_t c;
        assert(len(o) >= 4);
        get<0>(c) = extract<double>(o[0]);
        get<1>(c) = extract<double>(o[1]);
        get<2>(c) = extract<double>(o[2]);
        get<3>(c) = extract<double>(o[3]);
        void* storage =
            ( (boost::python::converter::rvalue_from_python_storage
               <color_t >*) data)->storage.bytes;
        new (storage) color_t(c);
        data->convertible = storage;
    }
};

struct color_vector_from_list
{
    color_vector_from_list()
    {
        converter::registry::push_back
            (&convertible, &construct,
             boost::python::type_id<vector<color_t> >());
    }

    static void* convertible(PyObject* obj_ptr)
    {
        handle<> x(borrowed(obj_ptr));
        object o(x);
        size_t N = len(o);
        if (N < 4)
            return 0;
        return obj_ptr;
    }

    static void construct(PyObject* obj_ptr,
                          converter::rvalue_from_python_stage1_data* data)
    {
        handle<> x(borrowed(obj_ptr));
        object o(x);
        vector<color_t> c;
        assert(len(o) >= 4);
        for (int i = 0; i < len(o) / 4; ++i)
        {
            c.push_back(color_t());
            get<0>(c[i]) = extract<double>(o[0 + 4 * i]);
            get<1>(c[i]) = extract<double>(o[1 + 4 * i]);
            get<2>(c[i]) = extract<double>(o[2 + 4 * i]);
            get<3>(c[i]) = extract<double>(o[3 + 4 * i]);
        }

        void* storage =
            ( (boost::python::converter::rvalue_from_python_storage
               <color_t >*) data)->storage.bytes;
        new (storage) vector<color_t>(c);
        data->convertible = storage;
    }
};


template <class Enum>
struct enum_from_int
{
    enum_from_int()
    {
        converter::registry::push_back
            (&convertible, &construct,
             boost::python::type_id<Enum>());
    }

    static void* convertible(PyObject* obj_ptr)
    {
        handle<> x(borrowed(obj_ptr));
        object o(x);
        extract<int> check(o);
        if (!check.check())
            return 0;
        return obj_ptr;
    }

    static void construct(PyObject* obj_ptr,
                          converter::rvalue_from_python_stage1_data* data)
    {
        handle<> x(borrowed(obj_ptr));
        object o(x);
        Enum val = Enum(extract<int>(o)());
        void* storage =
            ( (boost::python::converter::rvalue_from_python_storage
               <Enum>*) data)->storage.bytes;
        new (storage) Enum(val);
        data->convertible = storage;
    }
};


BOOST_PYTHON_MODULE(libgraph_tool_draw)
{
    def("cairo_draw", &cairo_draw);
    def("put_parallel_splines", &put_parallel_splines);
    def("apply_transforms", &apply_transforms);

    enum_<vertex_attr_t>("vertex_attrs")
        .value("shape", VERTEX_SHAPE)
        .value("color", VERTEX_COLOR)
        .value("fill_color", VERTEX_FILL_COLOR)
        .value("size", VERTEX_SIZE)
        .value("aspect", VERTEX_ASPECT)
        .value("anchor", VERTEX_ANCHOR)
        .value("pen_width", VERTEX_PENWIDTH)
        .value("halo", VERTEX_HALO)
        .value("halo_color", VERTEX_HALO_COLOR)
        .value("text", VERTEX_TEXT)
        .value("text_color", VERTEX_TEXT_COLOR)
        .value("text_position", VERTEX_TEXT_POSITION)
        .value("font_family", VERTEX_FONT_FAMILY)
        .value("font_slant", VERTEX_FONT_SLANT)
        .value("font_weight", VERTEX_FONT_WEIGHT)
        .value("font_size", VERTEX_FONT_SIZE)
        .value("surface", VERTEX_SURFACE)
        .value("pie_fractions", VERTEX_PIE_FRACTIONS)
        .value("pie_colors", VERTEX_PIE_COLORS);

    enum_<edge_attr_t>("edge_attrs")
        .value("color", EDGE_COLOR)
        .value("pen_width", EDGE_PENWIDTH)
        .value("start_marker", EDGE_START_MARKER)
        .value("mid_marker", EDGE_MID_MARKER)
        .value("end_marker", EDGE_END_MARKER)
        .value("marker_size", EDGE_MARKER_SIZE)
        .value("control_points", EDGE_CONTROL_POINTS)
        .value("dash_style", EDGE_DASH_STYLE)
        .value("sloppy", EDGE_SLOPPY);

    enum_<vertex_shape_t>("vertex_shape")
        .value("circle", SHAPE_CIRCLE)
        .value("triangle", SHAPE_TRIANGLE)
        .value("square", SHAPE_SQUARE)
        .value("pentagon", SHAPE_PENTAGON)
        .value("hexagon", SHAPE_HEXAGON)
        .value("heptagon", SHAPE_HEPTAGON)
        .value("octagon", SHAPE_OCTAGON)
        .value("double_circle", SHAPE_DOUBLE_CIRCLE)
        .value("double_triangle", SHAPE_DOUBLE_TRIANGLE)
        .value("double_square", SHAPE_DOUBLE_SQUARE)
        .value("double_pentagon", SHAPE_DOUBLE_PENTAGON)
        .value("double_hexagon", SHAPE_DOUBLE_HEXAGON)
        .value("double_heptagon", SHAPE_DOUBLE_HEPTAGON)
        .value("double_octagon", SHAPE_DOUBLE_OCTAGON)
        .value("pie", SHAPE_PIE);

    enum_<edge_marker_t>("edge_marker")
        .value("none", MARKER_SHAPE_NONE)
        .value("arrow", MARKER_SHAPE_ARROW)
        .value("circle", MARKER_SHAPE_CIRCLE)
        .value("square", MARKER_SHAPE_SQUARE)
        .value("diamond", MARKER_SHAPE_DIAMOND)
        .value("bar", MARKER_SHAPE_BAR);

    color_from_list();
    color_vector_from_list();
    enum_from_int<vertex_attr_t>();
    enum_from_int<edge_attr_t>();
    enum_from_int<vertex_shape_t>();
    enum_from_int<edge_marker_t>();
}

#else

#include <boost/python.hpp>
BOOST_PYTHON_MODULE(libgraph_tool_draw)
{
}

#endif // HAVE_CAIROMM
