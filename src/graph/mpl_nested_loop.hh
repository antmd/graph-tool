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

template <class Action, class... Args>
struct dispatch
{
    dispatch(Action a): _a(a) {}

    template <class T>
    dispatch<Action, Args..., T> join() const
    {
        return dispatch<Action, Args..., T>(_a);
    }

    template <class T>
    void operator()(T) const
    {
        _a(Args()..., T());
    }

    void operator()() const
    {
        _a(Args()...);
    }

    Action _a;
};

template <class TR1, class... TRS, class Action>
void nested_for_each_imp(Action a);

template <class Action, class... TRS>
struct inner_loop
{
    inner_loop(Action a): _a(a) {}

    template <class T>
    void operator()(T) const
    {
        nested_for_each_imp<TRS...>(_a.template join<T>());
    }

    Action _a;
};

template <class TR1, class... TRS, class Action>
void nested_for_each_imp(Action a)
{
    for_each<TR1>(inner_loop<Action, TRS...>(a));
}

template <class Action>
void nested_for_each_imp(Action a)
{
    a();
}

template <class... TRS, class Action>
void nested_for_each(Action a)
{
    auto b = dispatch<Action>(a);
    nested_for_each_imp<TRS...>(b);
}


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

template <class Action, class... Args>
selected_types<Action>
select_types(Action a, bool& found, Args... args)
{
    return selected_types<Action>(a, found, args...);
}

template <class Action>
struct selected_types
{
    template <class... Args>
    selected_types(Action a, bool& found, Args&&... args)
        : _a(a), _found(found)
    {
        _args = {args...};
    }


    template <class... Args, class T, class... Ts>
    void dispatch(unsigned int i, std::tuple<T, Ts...>, Args&&... args) const
    {
        assert(i < _args.size());
        T* a = const_cast<T*>(any_cast<T>(&_args[i]));
        if (a != 0)
            dispatch(i + 1, std::tuple<Ts...>(), std::forward<Args>(args)..., *a);
    }

    template <class... Args>
    void dispatch(unsigned int i, std::tuple<>, Args&&... args) const
    {
        _a(std::forward<Args>(args)...);
        _found = true;
    }

    template <class... Ts>
    void operator()(Ts&&... ts) const
    {
        dispatch(0, std::tuple<Ts...>());
    }

    Action _a;
    bool& _found;
    std::vector<any> _args;
};

} // mpl namespace
} // boost namespace

#endif //NESTED_FOR_LOOP_HH
