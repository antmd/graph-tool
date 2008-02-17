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


## FIXME

    @_handle_exceptions
    @_lazy_load
    def run_action(self, code, arg_names=[], local_dict=None,
                   global_dict=None, force=0, compiler="gcc", verbose=0,
                   auto_downcast=1, support_code="", libraries=[],
                   library_dirs=[], extra_compile_args=[],
                   runtime_library_dirs=[], extra_objects=[],
                   extra_link_args=[]):
        """Compile (if necessary) and run the C++ action specified by 'code',
        using weave."""
        try:
            import scipy.weave
        except ImportError:
            raise GraphError(self, "You need to have scipy installed to use" + \
                             " 'run_action'.")

        prefix_dir = libgraph_tool_core.mod_info().install_prefix
        python_dir = libgraph_tool_core.mod_info().python_dir
        python_dir = string.Template(python_dir).substitute(prefix=prefix_dir)
        cxxflags = libgraph_tool_core.mod_info().cxxflags

        # this is the code template which defines the action functor
        support_template = r"""
        #include <map>
        #include <set>
        #include <list>
        #include <tr1/unordered_set>
        #include <tr1/unordered_map>
        #include <boost/lambda/lambda.hpp>
        #include <boost/lambda/bind.hpp>
        #include <boost/tuple/tuple.hpp>
        #include <boost/type_traits.hpp>
        #include "${include_prefix}/graph.hh"
        #include "${include_prefix}/graph_filtering.hh"
        #include "${include_prefix}/graph_properties.hh"

        using namespace boost;
        using namespace boost::tuples;
        using namespace std;
        using namespace graph_tool;

        template <class IndexMap>
        struct prop_bind_t
        {
            template <class Value>
            struct as
            {
                typedef vector_property_map<Value,IndexMap> type;
            };
        };

        struct action_${code_hash}
        {
            template <class Graph, class VertexIndex, class EdgeIndex,
                      class Args>
            void operator()(Graph& g, VertexIndex vertex_index,
                            EdgeIndex edge_index,
                            dynamic_properties& properties,
                            const Args& args) const
            {
                // convenience typedefs
                typedef typename graph_traits<Graph>::vertex_descriptor
                    vertex_t;
                typedef typename graph_traits<Graph>::vertex_iterator
                    vertex_iter_t;
                typedef typename graph_traits<Graph>::edge_descriptor edge_t;
                typedef typename graph_traits<Graph>::edge_iterator edge_iter_t;
                typedef typename graph_traits<Graph>::out_edge_iterator
                    out_edge_iter_t;
                typedef typename in_edge_iteratorS<Graph>::type in_edge_iter_t;

                typedef prop_bind_t<VertexIndex> vertex_prop_t;
                typedef prop_bind_t<EdgeIndex> edge_prop_t;

                typedef typename vertex_prop_t::template as<bool>::type
                    vprop_bool_t;
                typedef typename vertex_prop_t::template as<int>::type
                    vprop_int_t;
                typedef typename vertex_prop_t::template as<long>::type
                    vprop_long_t;
                typedef typename vertex_prop_t::template as<long long>::type
                    vprop_long_long_t;
                typedef typename vertex_prop_t::template as<size_t>::type
                    vprop_size_t_t;
                typedef typename vertex_prop_t::template as<double>::type
                    vprop_double_t;
                typedef typename vertex_prop_t::template as<float>::type
                    vprop_float_t;
                typedef typename vertex_prop_t::template as<string>::type
                    vprop_string_t;

                typedef typename edge_prop_t::template as<bool>::type
                    eprop_bool_t;
                typedef typename edge_prop_t::template as<long>::type
                    eprop_long_t;
                typedef typename edge_prop_t::template as<int>::type
                    eprop_int_t;
                typedef typename edge_prop_t::template as<long long>::type
                    eprop_long_long_t;
                typedef typename edge_prop_t::template as<size_t>::type
                    eprop_size_t_t;
                typedef typename edge_prop_t::template as<double>::type
                    eprop_double_t;
                typedef typename edge_prop_t::template as<float>::type
                    eprop_float_t;
                typedef typename edge_prop_t::template as<string>::type
                    eprop_string_t;

                // the arguments will be expanded below
                ${arg_expansion}

                // the actual code
                ${code}
            }
        };

        """
        # we need to have different template names for each actions, to avoid
        # strange RTTI issues. We'll therefore add an md5 hash of the code to
        # each action's name
        import hashlib
        code_hash = hashlib.md5(code).hexdigest()

        # each term on the expansion will properly unwrap a tuple pointer value
        # to a reference with the appropriate name and type
        exp_term = """typename boost::remove_pointer<typename element<%d,Args>
                                                     ::type>::type& %s =
                          *get<%d>(args);"""
        arg_expansion = "\n".join([ exp_term % (i,arg_names[i],i) for i in \
                                    xrange(0, len(arg_names))])
        support_template = string.Template(support_template)
        inc_prefix = python_dir + "/graph_tool/include"
        support_code = support_template.substitute(code_hash=code_hash,
                                                   arg_expansion=arg_expansion,
                                                   code=code,
                                                   include_prefix = inc_prefix)\
                                                   + support_code

        # insert a hash value of the support_code into the code below, to force
        # recompilation when support_code (and module version) changes
        support_hash = hashlib.md5(support_code + __version__).hexdigest()

        # the actual inline code will just call g.RunAction() on the underlying
        # GraphInterface instance. The inline arguments will be packed into a
        # tuple of pointers.
        code = string.Template(r"""
        python::object pg(python::handle<>
                            (python::borrowed((PyObject*)(self___graph))));
        GraphInterface& g = python::extract<GraphInterface&>(pg);
        g.RunAction(action_${code_hash}(), make_tuple(${args}));
        // support code hash: ${support_hash}
        """).substitute(args=", ".join(["&%s" %a for a in arg_names]),
                        code_hash=code_hash, support_hash=support_hash)

        # we need to get the locals and globals of the _calling_ function. We
        # need to go deeper into the call stack due to all the function
        # decorators being used.
        call_frame = sys._getframe(5)
        if local_dict is None:
            local_dict = {}
            local_dict.update(call_frame.f_locals)
        if global_dict is None:
            global_dict = call_frame.f_globals
        local_dict["self___graph"] = self.__graph # the graph interface

        # RTLD_GLOBAL needs to be set in dlopen() if we want typeinfo and
        # friends to work properly across DSO boundaries. See
        # http://gcc.gnu.org/faq.html#dso
        sys.setdlopenflags(RTLD_NOW|RTLD_GLOBAL)

        # call weave and pass all the updated kw arguments
        scipy.weave.inline(code, ["self___graph"] + arg_names, force=force,
                           local_dict=local_dict, global_dict=global_dict,
                           compiler=compiler, verbose=verbose,
                           auto_downcast=auto_downcast,
                           support_code=support_code,
                           libraries=["graph_tool"] + libraries,
                           library_dirs=sys.path + library_dirs,
                           extra_compile_args=[cxxflags] + \
                                            extra_compile_args,
                           runtime_library_dirs=runtime_library_dirs,
                           extra_objects=extra_objects,
                           extra_link_args=["-L" + python_dir + "/graph_tool/",
                                            "-Wl,-E"] + extra_link_args)

        sys.setdlopenflags(_orig_dlopen_flags) # reset dlopen to normal case to
                                               # avoid unnecessary symbol
                                               # collision
