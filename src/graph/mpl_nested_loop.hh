// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007-2010 Tiago de Paula Peixoto <tiago@forked.de>
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

#ifndef NESTED_FOR_LOOP_HH
#define NESTED_FOR_LOOP_HH

#include <boost/mpl/for_each.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/empty.hpp>
#include <boost/any.hpp>

namespace boost
{
namespace mpl
{
// The following is a implementation of a nested for_each loop, which runs a
// given Action functor for each combination of its arguments, given by the type
// ranges, as such:
//
//     struct foo
//     {
//         template<class T1, class T2, class T3>
//         void operator()(T1, T2, T3) const
//         {
//             ...
//         }
//     };
//
//     ...
//
//     typedef mpl::vector<int,float,long> r1;
//     typedef mpl::vector<string,double> r2;
//     typedef mpl::vector<size_t,char> r3;
//
//     nested_for_each<r1,r2,r3>(foo());
//
// The code above will run foo::operator(T1, T2, T3) for all combinations given
// by r1, r2 and r3. This provides a more general compile-time to run-time
// meta-programming than the more simple mpl::for_each().
//
// Implementation notes: This supports up to 5 nested type ranges. In principle
// this could be generalized, and a lot of the boiler plate code in the
// implementation could be eliminated, by making more thorough use of MPL's
// meta-functions. My attempts to do this, however, tried to make use of
// lambda's bind function, which requires the arity of the function being bound
// to be known in advance. I tried to implement a bind_last() function, which
// always binds the last argument of a given functor, independently of its exact
// arity, but it made use of the typeof operator a lot, and at a given point
// caused GCC (4.2.2) to ICE... Thus, I decided to write the inner loops by
// hand, and wait before GCC implements the auto and typedecl functionality of
// c++0x, making the life of C++ template meta-programmers less horrible.
//
// In any case, "five nested levels should be enough for everybody".

template <class TR1, class TR2, class TR3, class TR4, class TR5>
struct nested_for_each
{
    template <class Action>
    void operator()(Action a) const
    {
        for_each<TR1>(inner_loop1<Action>(a));
    }

    template <class Action>
    struct inner_loop1
    {
        inner_loop1(Action a): _a(a) {}

        template <class T1>
        void operator()(T1) const
        {
            typedef typename mpl::if_<mpl::empty<TR3>,
                                      TR2,
                                      mpl::vector<> >::type eval_range;
            if(mpl::empty<TR3>::type::value)
                for_each<eval_range>(eval_action2<Action,T1>(_a));
            else
                for_each<TR2>(inner_loop2<Action,T1>(_a));
        }

        Action _a;
    };

    template<class Action, class T1>
    struct inner_loop2
    {
        inner_loop2(Action a): _a(a) {}

        template <class T2>
        void operator()(T2) const
        {
            typedef typename mpl::if_<mpl::empty<TR4>,
                                      TR3,
                                      mpl::vector<> >::type eval_range;
            if(mpl::empty<TR4>::type::value)
                for_each<eval_range>(eval_action3<Action,T1,T2>(_a));
            else
                for_each<TR3>(inner_loop3<Action, T1, T2>(_a));
        }

        Action _a;
    };

    template<class Action, class T1, class T2>
    struct inner_loop3
    {
        inner_loop3(Action a): _a(a) {}

        template <class T3>
        void operator()(T3) const
        {
            typedef typename mpl::if_<mpl::empty<TR5>,
                                      TR4,
                                      mpl::vector<> >::type eval_range;
            if(mpl::empty<TR5>::type::value)
                for_each<eval_range>(eval_action4<Action,T1,T2,T3>(_a));
            else
                for_each<TR4>(inner_loop4<Action, T1, T2, T3>(_a));
        }

        Action _a;
    };

    template<class Action, class T1, class T2, class T3>
    struct inner_loop4
    {
        inner_loop4(Action a): _a(a) {}

        template <class T4>
        void operator()(T4) const
        {
            for_each<TR5>(eval_action5<Action, T1, T2, T3, T4>(_a));
        }

        Action _a;
    };

    template<class Action, class T1>
    struct eval_action2
    {
        eval_action2(Action a): _a(a) {}

        template <class T2>
        void operator()(T2) const
        {
            _a(T1(), T2());
        }

        Action _a;
    };

    template<class Action, class T1, class T2>
    struct eval_action3
    {
        eval_action3(Action a): _a(a) {}

        template <class T3>
        void operator()(T3) const
        {
            _a(T1(), T2(), T3());
        }

        Action _a;
    };

    template<class Action, class T1, class T2, class T3>
    struct eval_action4
    {
        eval_action4(Action a): _a(a) {}

        template <class T4>
        void operator()(T4) const
        {
            _a(T1(), T2(), T3(), T4());
        }

        Action _a;
    };

    template<class Action, class T1, class T2, class T3, class T4>
    struct eval_action5
    {
        eval_action5(Action a): _a(a) {}

        template <class T5>
        void operator()(T5) const
        {
            _a(T1(), T2(), T3(), T4(), T5());
        }

        Action _a;
        };
};

template <class TR1, class TR2, class TR3 = mpl::vector<>,
          class TR4 = mpl::vector<>, class TR5 = mpl::vector<> >
struct nested_for_each;

// The functor below wraps another functor Action, but only calls it for the
// correct argument types, determined at runtime. Together with
// nested_for_each() above, it provides general means of selecting static
// polymorphic implementation of algorithms at run-time, as such:
//
//    struct foo
//    {
//        template <class T1, class T2>
//        void operator()(T1 x, T2 y) const
//        {
//            ... do something with x and y ...
//        }x
//    };
//
//    typedef mpl::vector<double,std::string,int> types;
//
//    any x = double(1.0)               // determined at run-time
//    any y = std::string("user input") // determined at run-time
//
//    bool found;
//    mpl::nested_for_each<types,types>(mpl::select_types(foo, x, y));
//
// The last line will call foo::operator()(double, std::string) passing the
// values of x and y, and the found variable will contain the value 'true'.

template <class Action> struct selected_types; // forward decl.

template <class Action>
selected_types<Action>
select_types(Action a, bool& found, any a1 = any(), any a2 = any(),
             any a3 = any(), any a4 = any(), any a5 = any())
{
    return selected_types<Action>(a, found, a1, a2, a3, a4, a5);
}

template <class Action>
struct selected_types
{
    selected_types(Action a, bool& found, any a1, any a2, any a3, any a4,
                   any a5)
        : _a(a), _found(found), _a1(a1), _a2(a2), _a3(a3), _a4(a4), _a5(a5) {}

    template <class T1>
    void operator()(T1) const
    {
        const T1* a1 = any_cast<T1>(&_a1);
        if (a1 != 0)
        {
            _a(*a1);
            _found = true;
        }
    }

    template <class T1, class T2>
    void operator()(T1, T2) const
    {
        const T1* a1 = any_cast<T1>(&_a1);
        const T2* a2 = any_cast<T2>(&_a2);
        if (a1 != 0 && a2 != 0)
        {
            _a(*a1, *a2);
            _found = true;
        }
    }

    template <class T1, class T2, class T3>
    void operator()(T1, T2, T3) const
    {
        const T1* a1 = any_cast<T1>(&_a1);
        const T2* a2 = any_cast<T2>(&_a2);
        const T3* a3 = any_cast<T3>(&_a3);

        if (a1 != 0 && a2 != 0 && a3 != 0)
        {
            _a(*a1, *a2, *a3);
            _found = true;
        }
    }

    template <class T1, class T2, class T3, class T4>
    void operator()(T1, T2, T3, T4) const
    {
        const T1* a1 = any_cast<T1>(&_a1);
        const T2* a2 = any_cast<T2>(&_a2);
        const T3* a3 = any_cast<T3>(&_a3);
        const T4* a4 = any_cast<T4>(&_a4);

        if (a1 != 0 && a2 != 0 && a3 != 0 && a4 != 0)
        {
            _a(*a1, *a2, *a3, *a4);
            _found = true;
        }

    }

    template <class T1, class T2, class T3, class T4, class T5>
    void operator()(T1, T2, T3, T4, T5) const
    {
        const T1* a1 = any_cast<T1>(&_a1);
        const T2* a2 = any_cast<T2>(&_a2);
        const T3* a3 = any_cast<T3>(&_a3);
        const T4* a4 = any_cast<T4>(&_a4);
        const T5* a5 = any_cast<T5>(&_a5);

        if (a1 != 0 && a2 != 0 && a3 != 0 && a4 != 0 && a5 != 0)
        {
            _a(*a1, *a2, *a3, *a4, *a5);
            _found = true;
        }
    }

    Action _a;
    bool& _found;
    any _a1, _a2, _a3, _a4, _a5;
};

} // mpl namespace
} // boost namespace

#endif //NESTED_FOR_LOOP_HH
