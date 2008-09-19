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

import sys, string, hashlib, os.path
from .. import core
from .. import libgraph_tool_core

try:
    import scipy.weave
except ImportError:
    raise libgraph_tool_core.raise_error\
          ("You need to have scipy installed to use 'run_action'.")

prefix = None
for d in [d + "/graph_tool" for d in sys.path]:
    if os.path.exists(d):
        prefix = d
        break

inc_prefix = prefix + "/include"
cxxflags = libgraph_tool_core.mod_info().cxxflags + " -I%s" % inc_prefix

# this is the code template which defines the action function object
support_template = open(prefix + "/run_action/run_action_template.hh").read()

# property map types
props = ""
for d in ["vertex", "edge", "graph"]:
    for t in core.value_types():
        props += (("typedef typename %s_prop_t::template as<%s >::type" + \
                   " %sprop_%s_t;\n") % \
                  (d,t,d[0],t.replace(" ","_").replace("::","_")\
                   .replace("<","_").replace(">","_"))).\
                   replace("__","_")

def inline(g, code, arg_names=[], local_dict=None,
           global_dict=None, force=0, compiler="gcc", verbose=0,
           auto_downcast=1, support_code="", libraries=[],
           library_dirs=[], extra_compile_args=[],
           runtime_library_dirs=[], extra_objects=[],
           extra_link_args=[], mask_ret=[]):
    """Compile (if necessary) and run the C++ action specified by 'code',
    using weave."""

    # we need to have different template names for each actions, to avoid
    # strange RTTI issues. We'll therefore append an md5 hash of the code (plus
    # grah_tool version string) to each action's name
    code_hash = hashlib.md5(code + core.__version__).hexdigest()

    # each term on the expansion will properly unwrap a tuple pointer value
    # to a reference with the appropriate name and type
    exp_term = """typename boost::remove_pointer<typename tr1::tuple_element<%d,Args>::type>::type& %s =
                          *tr1::get<%d>(_args);"""
    arg_expansion = "\n".join([ exp_term % (i,arg_names[i],i) for i in \
                                xrange(0, len(arg_names))])

    # handle returned values
    return_vals = ""
    for arg in arg_names:
        if arg not in mask_ret:
            return_vals += 'return_vals["%s"] = %s;\n' % (arg, arg)

    support_template = string.Template(globals()["support_template"])

    support_code += support_template.substitute(code_hash=code_hash,
                                                property_map_types=props,
                                                arg_expansion=arg_expansion,
                                                code=code,
                                                return_vals = return_vals)

    # insert a hash value of the support_code into the code below, to force
    # recompilation when support_code (and module version) changes
    support_hash = hashlib.md5(support_code + core.__version__).hexdigest()

    # the actual inline code will just call g.RunAction() on the underlying
    # GraphInterface instance. The inline arguments will be packed into a
    # tuple of pointers.
    code = string.Template(r"""
    python::object pg(python::handle<>
                        (python::borrowed((PyObject*)(self___graph))));
    GraphInterface& g = python::extract<GraphInterface&>(pg);
    RunAction(g, make_action(tr1::make_tuple(${args}), return_val));
    // support code hash: ${support_hash}
    """).substitute(args=", ".join(["&%s" %a for a in arg_names]),
                    code_hash=code_hash, support_hash=support_hash)

    # we need to get the locals and globals of the _calling_ function. Thus, we
    # need to go deeper into the call stack
    call_frame = sys._getframe(1)
    if local_dict is None:
        local_dict = call_frame.f_locals
    if global_dict is None:
        global_dict = call_frame.f_globals
    local_dict["self___graph"] = g._Graph__graph # the graph interface

    # RTLD_GLOBAL needs to be set in dlopen() if we want typeinfo and
    # friends to work properly across DSO boundaries. See
    # http://gcc.gnu.org/faq.html#dso
    orig_dlopen_flags = sys.getdlopenflags()
    sys.setdlopenflags(core.RTLD_NOW|core.RTLD_GLOBAL)

    # call weave and pass all the updated kw arguments
    ret_vals = \
             scipy.weave.inline(code, ["self___graph"] + arg_names, force=force,
                                local_dict=local_dict, global_dict=global_dict,
                                compiler=compiler, verbose=verbose,
                                auto_downcast=auto_downcast,
                                support_code=support_code,
                                libraries=libraries,
                                library_dirs=sys.path + library_dirs,
                                extra_compile_args=[cxxflags] + \
                                        extra_compile_args,
                                runtime_library_dirs=runtime_library_dirs,
                                extra_objects=extra_objects,
                                extra_link_args=extra_link_args)
    # check if exception was thrown
    if ret_vals["__exception_thrown"]:
        libgraph_tool_core.raise_error(ret_vals["__exception_error"])
    else:
        del ret_vals["__exception_thrown"]
        del ret_vals["__exception_error"]
    sys.setdlopenflags(orig_dlopen_flags) # reset dlopen to normal case to
                                          # avoid unnecessary symbol collision
    return ret_vals
