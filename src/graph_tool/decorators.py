#! /usr/bin/env python
#
# Copyright (C) 2007 Tiago de Paula Peixoto <tiago@forked.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Some useful decorators
"""

__author__="Tiago de Paula Peixoto <tiago@forked.de>"
__copyright__="Copyright 2008 Tiago de Paula Peixoto"
__license__="GPL version 3 or above"

import inspect, functools
import libgraph_tool_core as libcore

################################################################################
# Decorators
# Some useful function decorators which will be used
################################################################################

def _wraps(func):
    """This decorator works like the functools.wraps meta-decorator, but
    also preserves the function's argument signature. This uses eval, and is
    thus a bit of a hack, but there no better way I know of to do this."""
    def decorate(f):
        argspec = inspect.getargspec(func)
        args_call = inspect.formatargspec(argspec[0])
        argspec = inspect.formatargspec(argspec[0], defaults=argspec[3])
        argspec = argspec.lstrip("(").rstrip(")")
        wrap = eval("lambda %s: f%s" % (argspec, args_call), locals())
        return functools.wraps(func)(wrap)
    return decorate

def _attrs(**kwds):
    """Decorator which adds arbitrary attributes to methods"""
    def decorate(f):
        for k in kwds:
            setattr(f, k, kwds[k])
        return f
    return decorate

def _limit_args(allowed_vals):
    """Decorator which will limit the values of given arguments to a specified
    list of allowed values, and raise TypeError exception if the given value
    doesn't belong. 'allowed_vals' is a dict containing the allowed value list
    for each limited function argument."""
    def decorate(func):
        @_wraps(func)
        def wrap(*args, **kwargs):
            arg_names = inspect.getargspec(func)[0]
            arguments = zip(arg_names, args)
            arguments += [(k, kwargs[k]) for k in kwargs.keys()]
            for a in arguments:
                if allowed_vals.has_key(a[0]):
                    if a[1] not in allowed_vals[a[0]]:
                        raise TypeError("value for '%s' must be one of: %s" % \
                                         (a[0], ", ".join(allowed_vals[a[0]])))
            return func(*args, **kwargs)
        return wrap
    return decorate

def _require(arg_name, *allowed_types):
    """Decorator that lets you annotate function definitions with argument type
    requirements. These type requirements are automatically checked by the
    system at function invocation time."""
    def make_wrapper(f):
        if hasattr(f, "wrapped_args"):
            wrapped_args = getattr(f, "wrapped_args")
        else:
            code = f.func_code
            wrapped_args = list(code.co_varnames[:code.co_argcount])

        try:
            arg_index = wrapped_args.index(arg_name)
        except ValueError:
            raise NameError, arg_name

        @_wraps(f)
        def wrapper(*args, **kwargs):
            if len(args) > arg_index:
                arg = args[arg_index]
                if not isinstance(arg, allowed_types):
                    type_list = " or ".join(str(allowed_type) \
                                            for allowed_type in allowed_types)
                    raise TypeError, "Expected '%s' to be %s; was %s." % \
                          (arg_name, type_list, type(arg))
            else:
                if arg_name in kwargs:
                    arg = kwargs[arg_name]
                    if not isinstance(arg, allowed_types):
                        type_list = " or ".join(str(allowed_type) \
                                                for allowed_type in \
                                                allowed_types)
                        raise TypeError, "Expected '%s' to be %s; was %s." %\
                              (arg_name, type_list, type(arg))

            return f(*args, **kwargs)
        wrapper.wrapped_args = wrapped_args
        return wrapper
    return make_wrapper
