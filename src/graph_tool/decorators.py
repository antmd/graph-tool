#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
#
# Copyright (C) 2006-2015 Tiago de Paula Peixoto <tiago@skewed.de>
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

from __future__ import division, absolute_import, print_function

__author__ = "Tiago de Paula Peixoto <tiago@skewed.de>"
__copyright__ = "Copyright 2006-2015 Tiago de Paula Peixoto"
__license__ = "GPL version 3 or above"

import inspect
import functools
import sys

################################################################################
# Decorators
# Some useful function decorators which will be used
################################################################################

# exec statement in python 2.7 and exec() function in 3.2 are mutually exclusive
if sys.hexversion > 0x03000000:
    def exec_function(source, filename, global_map):
        exec(compile(source, filename, "exec"), global_map)
else:
    eval(compile("""\
def exec_function(source, filename, global_map):
    exec compile(source, filename, "exec") in global_map
""","<exec_function>", "exec"))


def _wraps(func):
    """This decorator works like the functools.wraps meta-decorator, but
    also preserves the function's argument signature. This uses eval, and is
    thus a bit of a hack, but there no better way I know of to do this."""
    def decorate(f):
        argspec = inspect.getargspec(func)
        ___wrap_defaults = defaults = argspec[-1]
        if defaults is not None:
            def_string = ["___wrap_defaults[%d]" % d for
                          d in range(len(defaults))]
            def_names = argspec[0][-len(defaults):]
        else:
            def_string = None
            def_names = None
        args_call = inspect.formatargspec(argspec[0], defaults=def_names)
        argspec = inspect.formatargspec(argspec[0], defaults=def_string)
        argspec = argspec.lstrip("(").rstrip(")")
        wf = "def %s(%s):\n    return f%s\n" % \
            (func.__name__, argspec, args_call)
        if def_string is not None:
            for d in def_string:
                wf = wf.replace("'%s'" % d, "%s" % d)
            for d in def_names:
                wf = wf.replace("'%s'" % d, "%s" % d)
        exec_function(wf, __file__, locals())
        return functools.wraps(func)(locals()[func.__name__])
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
            arguments = list(zip(arg_names, args))
            arguments += [(k, kwargs[k]) for k in list(kwargs.keys())]
            for a in arguments:
                if a[0] in allowed_vals:
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
            wrapped_args = f.wrapped_args
        else:
            code = f.__code__
            wrapped_args = list(code.co_varnames[:code.co_argcount])

        try:
            arg_index = wrapped_args.index(arg_name)
        except ValueError:
            raise NameError(arg_name)

        @_wraps(f)
        def wrapper(*args, **kwargs):
            if len(args) > arg_index:
                arg = args[arg_index]
                if not isinstance(arg, allowed_types):
                    type_list = " or ".join(str(allowed_type) \
                                            for allowed_type in allowed_types)
                    raise TypeError("Expected '%s' to be %s; was %s." % \
                                    (arg_name, type_list, type(arg)))
            else:
                if arg_name in kwargs:
                    arg = kwargs[arg_name]
                    if not isinstance(arg, allowed_types):
                        type_list = " or ".join(str(allowed_type) \
                                                for allowed_type in \
                                                allowed_types)
                        raise TypeError("Expected '%s' to be %s; was %s." %\
                                        (arg_name, type_list, type(arg)))

            return f(*args, **kwargs)
        wrapper.wrapped_args = wrapped_args
        return wrapper
    return make_wrapper
