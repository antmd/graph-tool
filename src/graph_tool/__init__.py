#! /usr/bin/env python
# graph_tool.py -- a general graph manipulation python module
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

__author__="Tiago de Paula Peixoto <tiago@forked.de>"
__copyright__="Copyright 2007 Tiago de Paula Peixoto"
__license__="GPL version 3 or above"
__URL__="http://graph-tool.forked.de"

__all__ = ["Graph", "GraphError"]

import sys
# RTLD_GLOBAL needs to be set in dlopen() if we want typeinfo and friends to
# work properly across DSO boundaries. See http://gcc.gnu.org/faq.html#dso

# The "except" is because the dl module raises a system error on ia64 and x86_64
# systems because "int" and addresses are different sizes.
try:
    from dl import RTLD_NOW, RTLD_GLOBAL
except ImportError:
    RTLD_NOW = 2
    RTLD_GLOBAL = 256
_orig_dlopen_flags = sys.getdlopenflags()

sys.setdlopenflags(RTLD_NOW|RTLD_GLOBAL)
import libgraph_tool
sys.setdlopenflags(_orig_dlopen_flags) # reset it to normal case to avoid
                                       # unnecessary symbol collision
__version__ = libgraph_tool.mod_info().version

import os, os.path, re, struct, fcntl, termios, gzip, bz2, string,\
       textwrap, time, signal, traceback, shutil, time, math, inspect, \
       functools, types

################################################################################
# Utility functions
################################################################################

def _open_file(name, mode="w"):
    """Open the file correctly, according to its file name"""
    if name == "-":
        if mode == "r":
            return sys.stdin
        else:
            return sys.stdout
    name = os.path.expanduser(name)
    if name.endswith("bz2"):
        return bz2.BZ2File(name, mode)
    if name.endswith("gz"):
        return gzip.GzipFile(name, mode)
    return file(name, mode)

def _eval_expr(expr, vars=dict(), edit_expr=False):
    """Evaluats a given python expression into a callable function, which
    returns the evaluation of the last statement. If exp is not a string, return
    it intact."""

    if expr.__class__ != str:
        return expr
    else:
        import random
        vars.update(dict([(x, eval("random."+x)) for x in dir(random) \
                          if not x.startswith("_")]))
        import math
        vars.update(dict([(x, eval("math."+x)) for x in dir(math) \
                          if not x.startswith("_")]))

        # Load some specific functions form scipy if available (and only on
        # demand)
        try:
            if 'gamaincc' in expr:
                from scipy.special import gammaincc
                vars['gamaincc'] = gammaincc
            if 'fsolve' in expr:
                from scipy.optimize import fsolve
                vars['fsolve'] = fsolve
        except ImportError:
            pass

        vars['inv_poisson'] = inv_poisson
        vars['inv_exponential'] = inv_exponential
        vars['inv_power_law'] = inv_power_law

        # Evaluate all statements
        r = re.compile(r"((('[^']*')|(\"[^\"']*\"))|([^;]+)(?=;|$))")
        statements = [x[0].strip() for x in r.findall(expr)]
        file_re = re.compile(r"^import:(.*)")
        for s in statements[:-1]:
            m = file_re.match(s)
            if m != None:
                # Import code from file
                file_name = m.group(1)
                exec _open_file(file_name.strip(), mode="r") in vars
            else:
                exec s in vars

        # Last statement must be returned as a function
        last = statements[-1].strip()
        if last.startswith("lambda ") or last.startswith("lambda:") or \
               last in vars:
            return eval(last, vars)
        else:
            if edit_expr == False:
                return eval("lambda: " + last, vars)
            else:
                if edit_expr != "g":
                    return eval("lambda %s, g: " % edit_expr + last, vars)
                else:
                    return eval("lambda g: " + last, vars)

def _get_mean(hist):
    """Get the mean, and the standard deviation fo the mean, for a given
    histogram."""
    avg, dev, count = 0.0, 0.0, 0.0
    try:
        for k, v in hist.iteritems():
            avg += k*v
            count += v
        avg /= count
        for k, v in hist.iteritems():
            dev += (k - avg)**2
        dev = math.sqrt(dev/(count**2))
    except ZeroDivisionError:
        avg = dev = float("nan") # nans are ok, since graph can be empty
    return (avg, dev)

def _degree(name):
    """Retrieve the degree type from string"""
    deg = name
    if name == "in-degree" or name == "in":
        deg = libgraph_tool.Degree.In
    if name == "out-degree" or name == "out":
        deg = libgraph_tool.Degree.Out
    if name == "total-degree" or name == "total":
        deg = libgraph_tool.Degree.Total
    return deg

def _parse_range(range):
    """Parse a range in the form 'lower[*] upper[*]' (where an '*' indicates
    a closed boundary) or 'comp val', where 'comp' is [==|=|!=|>=|<=|>|<].
    """
    try:
        rcomp = re.compile(r"(>|<|>=|<=|==|=|!=)\s*([^=<>!]+)")
        rrange = re.compile(r"([^\*]+)(\*?)\s+([^\*]+)(\*?)")
        inverted = False
        if rcomp.match(range):
            comp = rcomp.match(range).group(1)
            thres = rcomp.match(range).group(2)
            thres = float(thres)
            if "<" in comp:
                ran = (float('-Inf'), thres)
                if "=" in comp:
                    inc = (False, True)
                else:
                    inc = (False, False)
            elif ">" in comp:
                ran = (thres, float('Inf'))
                if "=" in comp:
                    inc = (True, False)
                else:
                    inc = (False, False)
            elif comp == "==" or comp == "=" or comp == "!=":
                ran = (thres, thres)
                inc = (True, True)
                if comp == "!=":
                    inverted = True
        elif rrange.match(range):
            m = rrange.match(range)
            r1 = float(m.group(1))
            i1 = "*" in m.group(2)
            r2 = float(m.group(3))
            i2 = "*" in m.group(4)
            ran = (r1, r2)
            inc = (i1, i2)
        else:
            raise ValueError("invalid value for range: " + range)
    except (ValueError, TypeError), e:
        raise ValueError("invalid value for range: %s: %s " % (range, str(e)))
    return (ran, inc, inverted)

# Some other predefined functions
# TODO: put them in a more convenient place... an util submodule, perhaps?
def inv_poisson(p, m):
    """returns the inverse of the poisson distribution, with average m"""
    from scipy.optimize import fsolve
    from scipy.special import gammaincc
    return int(round(fsolve(lambda k, l: gammaincc(k, l)-p, m, (m))-0.5+1e-15))
    # FIXME: the 1e-15 constant hack is ugly

def inv_exponential(p, m):
    """returns the inverse of the discrete exponential distibution, with average
    m"""
    return int(round(log(1-p)/log(float(m)/(m+1))-0.5+1e-15))
    # FIXME: the 1e-15 constant hack is ugly

def inv_power_law(p, b):
    """returns the inverse of the discrete power-law distibution, with exponent
    b"""
    return int(round((1-p)**(-1/(b-1)) - 1))

################################################################################
# Decorators
# Some useful function decorators which will be used below
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

def _handle_exceptions(func):
    """Decorator which will catch and properly propagate exceptions raised
    by GraphInterface in all the methods"""
    @_wraps(func)
    def wrap(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except (IOError, RuntimeError), e:
            raise GraphError(args[0], str(e))
    return wrap

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

def _lazy_load(func):
    """Decorator which will call the 'load' method before executing the
    decorated function, thus implementing 'lazy loading'."""
    @_wraps(func)
    def wrap(*args, **kwargs):
        self = args[0]
        if self.lazy_filename != None and self.lazy_format != None:
            self.load(self.lazy_filename, self.lazy_format)
            self.lazy_filename = None
            self.lazy_format = None
        return func(*args, **kwargs)
    return wrap

################################################################################
# Graph class
# The main graph interface
################################################################################

class GraphError(Exception):
    """Exception raised by the Graph class"""
    def __init__(self, graph, message):
        self.graph = graph
        self.message = message
    def __str__(self):
        return self.message

class Graph(object):
    """The main graph type which encapsulates the GraphInterface"""

    def __init__(self):
        self.__graph = libgraph_tool.GraphInterface()
        self.lazy_filename = None
        self.lazy_format = None

    # Graph access

    @_handle_exceptions
    @_lazy_load
    def vertices(self):
        "Return iterator over the vertices"
        return self.__graph.Vertices()

    @_handle_exceptions
    @_lazy_load
    def edges(self):
        "Return iterator over the edges"
        return self.__graph.Edges()

    @_handle_exceptions
    @_lazy_load
    def add_vertex(self):
        "Add a new vertex to the graph, and return it"
        return self.__graph.AddVertex()

    @_handle_exceptions
    @_lazy_load
    def remove_vertex(self, vertex):
        "Remove a vertex from the graph"
        return self.__graph.RemoveVertex(vertex)

    @_handle_exceptions
    @_lazy_load
    def add_edge(self, source, target):
        "Add a new edge from 'source' to 'target' to the graph, and return it"
        return self.__graph.AddEdge(source, target)

    @_handle_exceptions
    @_lazy_load
    def remove_edge(self, edge):
        "Remove a edge from the graph"
        return self.__graph.RemoveEdge(edge)

    @_handle_exceptions
    @_lazy_load
    def __get_vertex_properties(self):
        return self.__graph.GetVertexProperties()
    vertex_properties = property(__get_vertex_properties,
                                 doc="Dictionary of vertex properties")

    @_handle_exceptions
    @_lazy_load
    def __get_edge_properties(self):
        return self.__graph.GetEdgeProperties()
    edge_properties = property(__get_edge_properties,
                                 doc="Dictionary of edge properties")

    @_handle_exceptions
    @_lazy_load
    def __get_graph_properties(self):
        class GraphPropertyDict(dict):
            """Wrapper for the dict of graph properties, which sets the value on
            the property map when changed in the dict."""
            def __init__(self, g, old):
                dict.__init__(self)
                dict.update(self,old)
                self.g = g
            def __setitem__(self, key, val):
                g = self.g._Graph__graph
                try:
                    g.GetGraphProperties()[key][g] = val
                except KeyError:
                    raise
                except:
                    raise ValueError("wrong value for type '%s': '%s'" % \
                                     (g.GetGraphProperties()[key].get_type(),
                                      str(val)))
                dict.__setitem__(self, key, val)
            def __delitem__(self, key):
                g = self.g._Graph__graph
                g.RemoveGraphProperty(key)
                dict.__delitem__(self, key)
            def get_type(self, key):
                """Return the type of the internal property"""
                return self.g._Graph__graph.GetGraphProperties()[key].get_type()
        valdict = dict((k, v[self.__graph]) \
                       for k,v in self.__graph.GetGraphProperties().iteritems())
        return GraphPropertyDict(self, valdict)

    @_handle_exceptions
    @_lazy_load
    def __set_graph_properties(self, val):
        props = self.__graph.GetGraphProperties()
        for k,v in val.iteritems():
            if not props.has_key(k):
                if v.__class__ == str:
                    typename = "string"
                elif v.__class__ == float:
                    typename = "double"
                elif v.__class__ == int:
                    typename = "int"
                elif v.__class__ == bool:
                    typename = "boolean"
                else:
                    raise TypeError("invalid value type for graph " + \
                                    "property '%s'" % k)
                self.add_graph_property(k, typename)
                props = self.__graph.GetGraphProperties()
            props[k][self.__graph] = v

    graph_properties = property(__get_graph_properties, __set_graph_properties,
                                doc="Dictionary of graph properties")

    @_handle_exceptions
    @_lazy_load
    def add_vertex_property(self, name, type):
        """Add a new (uninitialized) vertex property map of type 'type', and
        return it"""
        self.edit_vertex_property(name, type)
        return self.vertex_properties[name]

    @_handle_exceptions
    @_lazy_load
    def add_edge_property(self, name, type):
        """Add a new (uninitialized) edge property map of type 'type', and
        return it"""
        self.edit_edge_property(name, type)
        return self.edge_properties[name]

    @_handle_exceptions
    @_lazy_load
    def add_graph_property(self, name, type):
        """Add a new (uninitialized) graph property map of type 'type'"""
        self.edit_graph_property(name, type)

    # Basic options (File i/o and such)
    __groups = ["Basic Options"]

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    def load(self, filename, format="auto"):
        """Load graph from 'filename'. The format is guessed from the file name,
        or can be specified by 'format', which can be either 'xml' or 'dot'."""

        if format == 'auto':
            self.__graph.ReadFromFile(filename)
        else:
            self.__graph.ReadFromFile(filename, format)
        self.lazy_filename = filename
        self.lazy_format = format

    @_handle_exceptions
    def lazy_load(self, filename, format="auto"):
        """Calls load() with the same arguments when the graph is first accessed."""

        self.lazy_filename = filename
        self.lazy_format = format

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def save(self, filename, format="auto"):
        """Save graph to file. The format is guessed from the 'file' name, or
        can be specified by 'format', which can be either 'xml' or 'dot'."""

        if format == 'auto':
            self.__graph.WriteToFile(filename)
        else:
            self.__graph.WriteToFile(filename, format)

    # Graph generation
    __groups.append("Graph generation")

    @_attrs(opt_group=__groups[-1], first_subopt="N")
    @_handle_exceptions
    @_lazy_load
    def correlated_configurational_model \
            (self, N = 10000, pjk = "lambda j, k: 1.0",
             pjk_ceil = "lambda j, k: pjk(j, k)", pjk_m = 1.0,
             inv_pjk_ceil = "lambda p,r:(inv_poisson(p, 2), inv_poisson(p, 2))",
             corr = "lambda jl, kl, j, k: 1.0",
             corr_ceil = "lambda jl, kl, j, k: corr(jl,kl,j,k)", corr_m = 1.0,
             inv_corr_ceil = "lambda p, r, a, b: inv_pjk_ceil(p, r)",
             undirected = False, progress = False, seed = int(time.time())):
        """Generate graph using the configurational model with arbitrary degree
        correlations. See documentation for details."""

        f = self.__graph.GenerateCorrelatedConfigurationalModel
        exp_vars = dict()
        exp_vars["pjk"] = _eval_expr(pjk, exp_vars)
        exp_vars["pjk_ceil"] = _eval_expr(pjk_ceil, exp_vars)
        exp_vars["inv_pjk_ceil"] = _eval_expr(inv_pjk_ceil, exp_vars)
        exp_vars["corr"] = _eval_expr(corr, exp_vars)
        exp_vars["corr_ceil"] = _eval_expr(corr_ceil, exp_vars)
        exp_vars["inv_corr_ceil"] = _eval_expr(inv_corr_ceil, exp_vars)
        f(N, exp_vars["pjk"], exp_vars["pjk_ceil"], exp_vars["inv_pjk_ceil"],
          pjk_m, exp_vars["corr"], exp_vars["corr_ceil"],
          exp_vars["inv_corr_ceil"], corr_m, undirected, seed, progress)

    # Graph filtering
    __groups.append("Filtering")

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def directed(self):
        """Treat graph as directed (default)."""
        self.__graph.SetDirected(True)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def undirected(self):
        """Treat graph as undirected."""
        self.__graph.SetDirected(False)

    @_handle_exceptions
    @_lazy_load
    def set_directed(self, is_directed):
        """Set the directedness of the graph"""
        self.__graph.SetDirected(is_directed)

    @_handle_exceptions
    @_lazy_load
    def is_directed(self):
        """Get the directedness of the graph"""
        return self.__graph.GetDirected()

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def reverse(self):
        """Reverse the direction of the edges."""
        self.__graph.SetReversed(not self.__graph.GetReversed())

    @_handle_exceptions
    @_lazy_load
    def set_reverse(self, is_reversed):
        """Reverse the direction of the edges, if 'reversed' is True, or
        maintain the original direction otherwise."""
        self.__graph.SetReversed(is_reversed)

    @_handle_exceptions
    @_lazy_load
    def get_reverse(self):
        """Return 'True' if the edges are reversed, and 'False' otherwise."""
        return self.__graph.GetReversed()

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def exclude_vertex_range(self, property, range):
        """Choose vertex property and range to exclude."""
        (ran, inc, inverted) = _parse_range(range)
        self.__graph.SetVertexFilterProperty(property)
        self.__graph.SetVertexFilterRange(ran, inc, not inverted)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def keep_vertex_range(self, property, range):
        """Choose vertex property and range to keep (and exclude the rest)."""
        (ran, inc, inverted) = _parse_range(range)
        self.__graph.SetVertexFilterProperty(property)
        self.__graph.SetVertexFilterRange(ran, inc, inverted)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def reset_vertex_filter(self):
        """Remove edge filter."""
        self.__graph.SetVertexFilterProperty('')

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def exclude_edge_range(self, property, range):
        """Choose edge property and range to exclude."""
        (ran, inc, inverted) = _parse_range(range)
        self.__graph.SetEdgeFilterProperty(property)
        self.__graph.SetEdgeFilterRange(ran, inc, not inverted)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def keep_edge_range(self, property, range):
        """Choose edge property and range to keep (and exclude the rest)."""
        (ran, inc, inverted) = _parse_range(range)
        self.__graph.SetEdgeFilterProperty(property)
        self.__graph.SetEdgeFilterRange(ran, inc, inverted)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    def reset_edge_filter(self):
        """Remove edge filter"""
        self.__graph.SetEdgeFilterProperty('')

    # Graph modification
    __groups.append("Graph Modification")

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def edit_vertex_property(self, property, type="double", expr=None):
        """Edit the selected vertex property"""
        self.__graph.EditVertexProperty(property, type,
                                        _eval_expr(expr,
                                                   edit_expr="v"), self)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def edit_edge_property(self, property, type="double", expr=None):
        """Edit the selected edge property"""
        self.__graph.EditEdgeProperty(property, type,
                                      _eval_expr(expr,
                                                 edit_expr="e"), self)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def edit_graph_property(self, property, type="double", expr=None):
        """Edit the selected graph property"""
        self.__graph.EditGraphProperty(property, type,
                                       _eval_expr(expr,
                                                  edit_expr="g"), self)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def remove_vertex_property(self, property):
        """Remove the selected vertex property"""
        self.__graph.RemoveVertexProperty(property)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def remove_edge_property(self, property):
        """Remove the selected edge property"""
        self.__graph.RemoveEdgeProperty(property)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def remove_graph_property(self, property):
        """Remove the selected graph property"""
        self.__graph.RemoveGraphProperty(property)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def insert_vertex_index_property(self, property):
        """Insert vertex index as property"""
        self.__graph.InsertVertexIndexProperty(property)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def insert_edge_index_property(self, property):
        """Insert edge index as property"""
        self.__graph.InsertEdgeIndexProperty(property)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def list_properties(self):
        """List all properties"""
        w = max([len(x) for x in self.graph_properties.keys() + \
                 self.vertex_properties.keys() + \
                 self.edge_properties.keys()]) + 4
        w = w if w > 14 else 14
        for k,v in self.graph_properties.iteritems():
            print "%%-%ds (graph)   (type: %%s, val: %%s)" % w % \
                  (k, self.graph_properties.get_type(k), str(v))
        for k,v in self.vertex_properties.iteritems():
            print "%%-%ds (vertex)  (type: %%s)" % w % (k, v.get_type())
        for k,v in self.edge_properties.iteritems():
            print "%%-%ds (edge)    (type: %%s)" % w % (k, v.get_type())

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def purge_vertices(self):
        """Remove all vertices of the graph which are currently being filtered
        out, and return to the unfiltered state"""
        self.__graph.PurgeVertices()
        self.__graph.SetVertexFilterProperty('')

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def purge_edges(self):
        """Remove all edges of the graph which are currently being filtered out,
        and return to the unfiltered state"""
        self.__graph.PurgeEdges()
        self.__graph.SetEdgeFilterProperty('')

    @_attrs(opt_group=__groups[-1], fist_subopt="parallel_edges")
    @_handle_exceptions
    @_lazy_load
    @_limit_args({"strategy":["correlated", "uncorrelated"]})
    def random_rewire(self, strategy="uncorrelated", parallel_edges=False,
                      self_loops=False, seed=int(time.time())):
        """Randomly rewire the edges of the graph"""
        self.__graph.RandomRewire(strategy, self_loops, parallel_edges, seed)

    # Basic graph statistics
    __groups.append("Basic Statistics")

    @_attrs(opt_group=__groups[-1], has_output=True)
    @_handle_exceptions
    @_lazy_load
    def number_of_vertices(self):
        """Get the number of vertices."""
        return self.__graph.GetNumberOfVertices()

    @_attrs(opt_group=__groups[-1], has_output=True)
    @_handle_exceptions
    @_lazy_load
    def number_of_edges(self):
        """Get the number of edges"""
        return self.__graph.GetNumberOfEdges()

    @_attrs(opt_group=__groups[-1], has_output=True)
    @_handle_exceptions
    @_lazy_load
    def vertex_histogram(self, degree):
        """Get the vertex degree/property histogram"""
        return self.__graph.GetVertexHistogram(_degree(degree))

    @_attrs(opt_group=__groups[-1], has_output=True)
    @_handle_exceptions
    @_lazy_load
    def edge_histogram(self, property):
        """Get the edge property histogram"""
        return self.__graph.GetEdgeHistogram(property)

    @_attrs(opt_group=__groups[-1], has_output=True)
    @_handle_exceptions
    @_lazy_load
    def average_vertex_property(self, degree):
        """Get the average of the vertex property"""
        return _get_mean(self.__graph.GetVertexHistogram(_degree(degree)))

    @_attrs(opt_group=__groups[-1], has_output=True)
    @_handle_exceptions
    @_lazy_load
    def average_edge_property(self, property):
        """Get the average of the edge property"""
        return _get_mean(self.__graph.GetEdgeHistogram(property))

    @_attrs(opt_group=__groups[-1], has_output=True)
    @_handle_exceptions
    @_lazy_load
    def combined_vertex_histogram(self, degree1, degree2):
        """Get the combined (degree1, degree2) histogram. Scalar properties are
        also accepted as degree1 or degree2"""
        return self.__graph.GetCombinedVertexHistogram(_degree(degree1),
                                                       _degree(degree2))

    @_attrs(opt_group=__groups[-1], has_output=True)
    @_handle_exceptions
    @_lazy_load
    def distance_histogram(self, weight=None):
        """Get the distance histogram"""
        if weight == None:
            weight = ""
        self.__graph.GetDistanceHistogram(weight)

    @_attrs(opt_group=__groups[-1], has_output=True)
    @_handle_exceptions
    @_lazy_load
    def average_distance(self, weight=None):
        """Get the averarge distance"""
        if weight == None:
            weight = ""
        return _get_mean(self.__graph.GetDistanceHistogram(weight))

    @_attrs(opt_group=__groups[-1], has_output=True)
    @_handle_exceptions
    @_lazy_load
    def average_harmonic_distance(self, weight=None):
        """Get the averarge harmonic distance"""
        if weight == None:
            weight = ""
        hist = self.__graph.GetDistanceHistogram(weight)
        avg, err = _get_mean(dict((1.0/k, v) for k, v in hist.iteritems()))
        return (1.0/avg, err/(avg**2))

    @_attrs(opt_group=__groups[-1], fist_subopt="samples", has_output=True)
    @_handle_exceptions
    @_lazy_load
    def sampled_distance_histogram(self, samples=1000, weight=None,
                                   seed=int(time.time())):
        """Get the sampled distance histogram"""
        if weight == None:
            weight = ""
        return self.__graph.GetSampledDistanceHistogram(weight, samples, seed)

    @_attrs(opt_group=__groups[-1], has_output=True, first_subopt="samples")
    @_handle_exceptions
    @_lazy_load
    def average_sampled_distance(self, samples=1000, weight=None,
                                 seed=int(time.time())):
        """Get the average sampled distance."""
        if weight == None:
            weight = ""
        return _get_mean(self.__graph.GetSampledDistanceHistogram(weight,
                                                                  samples,
                                                                  seed))

    @_attrs(opt_group=__groups[-1], has_output=True, first_subopt="samples")
    @_handle_exceptions
    @_lazy_load
    def average_sampled_harmonic_distance(self, samples=1000, weight=None,
                                 seed=int(time.time())):
        """Get the average sampled harmonic distance."""
        if weight == None:
            weight = ""
        hist = self.__graph.GetSampledDistanceHistogram(weight, samples, seed)
        avg, err = _get_mean(dict([(1.0/k, v) for k, v in hist.iteritems()]))
        (1.0/avg, err/(avg**2))

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def label_components(self, property):
        """Label components to property"""
        self.__graph.LabelComponents(property)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def label_parallel_edges(self, property):
        """Label parallel edges to property"""
        self.__graph.LabelParallelEdges(property)

    @_attrs(opt_group=__groups[-1], has_output=True)
    @_handle_exceptions
    @_lazy_load
    def reciprocity(self):
        """Get the edge reciprocity"""
        self.__graph.GetReciprocity()

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def minimum_spanning_tree(self, property,  weight=None):
        """Mark the minimum spanning tree edges in property"""
        if weight == None:
            weight = ""
        self.__graph.GetMinimumSpanningTree(weight, property)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def line_graph(self, file, format="xml"):
        """Save the corresponding line graph to file"""
        self.__graph.GetLineGraph(file, format)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def betweenness_centrality(self, vertex_betweeness, edge_betweeness,
                               weight=None):
        """Calculate and store the vertex and/or edge betweenness centrality"""
        if weight == None:
            weight = ""
        self.__graph.GetBetweenness(weight, edge_betweeness, vertex_betweeness)

    @_attrs(opt_group=__groups[-1], has_output=True)
    @_handle_exceptions
    @_lazy_load
    def central_point_dominance(self, vertex_betweeness):
        """Calculate central point dominance, given the 'vertex_betweeness'
        vertex property"""
        return self.__graph.GetCentralPointDominance(vertex_betweeness)

    __groups.append("Correlations")

    @_attrs(opt_group=__groups[-1], has_output=True)
    @_handle_exceptions
    @_lazy_load
    def average_combined_vertex_correlation(self, degree1, degree2):
        """Get the average of degree2 in function of degree1. Scalar properties
        are also accepted as degree1 or degree2"""
        f = self.__graph.GetAverageCombinedVertexCorrelation
        return f(_degree(degree1), _degree(degree2))

    @_attrs(opt_group=__groups[-1], has_output=True)
    @_handle_exceptions
    @_lazy_load
    def vertex_correlation_histogram(self, origin_degree, neighbour_degree,
                                     weight=None):
        """Get the degree correlation histogram. Scalar properties are also
        accepted in place of *_degree. An optional edge weight property can be
        passed by weight"""
        if weight == None:
            weight = ""
        f = self.__graph.GetVertexCorrelationHistogram
        return f(_degree(origin_degree), _degree(neighbour_degree), weight)

    @_attrs(opt_group=__groups[-1], has_output=True)
    @_handle_exceptions
    @_lazy_load
    def average_nearest_neighbours_correlation(self, origin_degree,
                                               neighbour_degree, weight=None):
        """Get the average nearest neighbours correlation. Scalar properties are
        also accepted in place of *_degree. An optional edge weight property can
        be passed by 'weight'"""
        if weight == None:
            weight = ""
        f = self.__graph.GetAverageNearestNeighboursCorrelation
        return f(_degree(origin_degree), _degree(neighbour_degree), weight)

    @_attrs(opt_group=__groups[-1], has_output=True)
    @_handle_exceptions
    @_lazy_load
    def edge_vertex_correlation_histogram(self, degree_source, edge_prop,
                                          degree_target):
        """Get the source degree vs. edge scalar vs. target degree correlation
        histogram. Scalar properties are also accepted in place of degree_*"""
        f = self.__graph.GetEdgeVertexCorrelationHistogram
        return f(_degree(degree_source), edge_prop, _degree(degree_target))

    @_attrs(opt_group=__groups[-1], has_output=True)
    @_handle_exceptions
    @_lazy_load
    def assortativity_coefficient(self, degree):
        """Get the assortativity coefficient. Scalar properties are also
        accepted in place of degree"""
        return self.__graph.GetAssortativityCoefficient(_degree(degree))

    @_attrs(opt_group=__groups[-1], has_output=True)
    @_handle_exceptions
    @_lazy_load
    def scalar_assortativity_coefficient(self, degree):
        """Get the scalar assortativity coefficient. Scalar properties are also
        accepted in place of degree"""
        return self.__graph.GetScalarAssortativityCoefficient(_degree(degree))

    __groups.append("Clustering")

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def local_clustering_coefficient(self, property):
        """Set the local clustering coefficient to vertex property"""
        self.__graph.SetLocalClusteringToProperty(property)

    @_attrs(opt_group=__groups[-1], has_output=True)
    @_handle_exceptions
    @_lazy_load
    def global_clustering_coefficient(self, file):
        """Get the global clustering coefficient"""
        return self.__graph.GetGlobalClustering()

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def extended_clustering_coefficient(self, prefix, max):
        """Set the extended clustering coefficients c1 to cmax to vertex
        properties prefix1 to prefixmax"""
        self.__graph.SetExtendedClusteringToProperty(prefix, max)

    __groups.append("Layout")

    @_attrs(opt_group=__groups[-1], first_subopt="type")
    @_handle_exceptions
    @_lazy_load
    @_limit_args({"type":["fg-all-pairs", "fg-grid", "kw"]})
    def spring_block_layout(self, property, type="fg-all-pairs", iter=1000,
                            weight=None, seed=int(time.time())):
        """Compute the spring block layout. The positions will be stored in
        property."""
        if weight == None:
            weight = ""
        self.__graph.ComputeGraphLayoutSpringBlock(property, weight, type, iter,
                                                   seed)

    @_attrs(opt_group=__groups[-1], first_subopt="type")
    @_handle_exceptions
    @_lazy_load
    @_limit_args({"type":["square", "circle", "heart"]})
    def gursoy_atun_layout(self, property, type="square", iter=1000,
                           weight=None, seed=int(time.time())):
        """Compute the Gursoy-Atun layout."""
        if weight == None:
            weight = ""
        self.__graph.ComputeGraphLayoutGursoy(property, weight,type, iter, seed)

    __groups.append("Community")

    @_attrs(opt_group=__groups[-1], first_subopt="g")
    @_handle_exceptions
    @_lazy_load
    @_limit_args({"corr":["correlated", "uncorrelated", "random"]})
    def community_structure(self, property, g=1.0, n=1000, tmin=0.01, tmax=1.0,
                            spins=0, corr="uncorrelated", weight=None,
                            seed=int(time.time()), verbose=False, history=""):
        """Calculate the community structure and assign it to property."""
        strat_map = {'correlated': libgraph_tool.CommCorr.Correlated,
                     'uncorrelated': libgraph_tool.CommCorr.Uncorrelated,
                     'random': libgraph_tool.CommCorr.ErdosReyni}
        if weight == None:
            weight = ""
        self.__graph.GetCommunityStructure(g, strat_map[correlated], n, tmin,
                                           tmax, spins, seed, verbose, history,
                                           weight, property)

    @_attrs(opt_group=__groups[-1], has_output=True)
    @_handle_exceptions
    @_lazy_load
    def modularity(self, property, weight=None):
        """Calculate the modularity, given a community partition specified by
        property"""
        if weight == None:
            weight = ""
        self.__graph.GetModularity(weight, property)

    @_attrs(opt_group=__groups[-1])
    @_handle_exceptions
    @_lazy_load
    def community_graph(self, property, size_property, file, format):
        """Obtain the graph of communities, given a community partition
        specified by property. The resulting graph will have a vertex property
        named SIZE-property, which will contain the size of the corresponding
        communities."""
        if format == 'auto':
            format = ''
        self.__graph.GetCommunityNetwork(property, size_property, file, format)

    __groups.append("Plugins")

    @_attrs(opt_group=__groups[-1], first_subopt="arg_names")
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

        prefix_dir = libgraph_tool.mod_info().install_prefix
        python_dir = libgraph_tool.mod_info().python_dir
        python_dir = string.Template(python_dir).substitute(prefix=prefix_dir)
        cxxflags = libgraph_tool.mod_info().cxxflags
        
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

