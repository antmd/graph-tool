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

import sys
# RTLD_GLOBAL needs to be set in dlopen() if we want typeinfo and friends to
# work properly across DSO boundaries. See http://gcc.gnu.org/faq.html#dso

# The "except" is because the dl module raises a system error on ia64 and x86_64
# systems because "int" and addresses are different sizes.
try:
    from dl import RTLD_LAZY, RTLD_NOW, RTLD_GLOBAL
except ImportError:
    RTLD_LAZY = 1
    RTLD_NOW = 2
    RTLD_GLOBAL = 256
_orig_dlopen_flags = sys.getdlopenflags()

sys.setdlopenflags(RTLD_LAZY|RTLD_GLOBAL)
import libgraph_tool_correlations
sys.setdlopenflags(_orig_dlopen_flags) # reset it to normal case to avoid
                                       # unnecessary symbol collision

from .. core import _degree, _prop
from numpy import *

__all__ = ["assortativity", "scalar_assortativity",
           "corr_hist", "combined_corr_hist", "avg_neighbour_corr"]

def assortativity(g, deg):
    return libgraph_tool_correlations.\
           assortativity_coefficient(g.underlying_graph(), _degree(g, deg))

def scalar_assortativity(g, deg):
    return libgraph_tool_correlations.\
           scalar_assortativity_coefficient(g.underlying_graph(),
                                            _degree(g, deg))

def corr_hist(g, deg1, deg2, bins=[[1],[1]], weight=None):
    ret = libgraph_tool_correlations.\
          vertex_correlation_histogram(g.underlying_graph(), _degree(g, deg1),
                                       _degree(g, deg2), _prop("e", g, weight),
                                       bins[0], bins[1])
    return [ret[0], [ret[1][0], ret[1][1]]]

def combined_corr_hist(g, deg1, deg2, bins=[[1],[1]]):
    ret = libgraph_tool_correlations.\
          vertex_combined_correlation_histogram(g.underlying_graph(),
                                                _degree(g, deg1),
                                                _degree(g, deg2),
                                                bins[0], bins[1])
    return [ret[0], [ret[1][0], ret[1][1]]]

def avg_neighbour_corr(g, deg1, deg2, bins=[[1],[1]], weight=None):
    ret = corr_hist(g, deg1, deg2, bins, weight)
    xbins = ret[1][0]
    ybins = ret[1][1]
    counts = ret[0]
    avg = empty((counts.shape[0]))
    dev = empty((counts.shape[0]))
    mask = empty((counts.shape[0]), dtype=dtype('bool'))
    n_masked = 0
    for i in xrange(0, len(ret[0])):
        N = counts[i,:].sum()
        if N > 0:
            avg[i] = average(ybins, weights=counts[i,:])
            dev[i] = sqrt(average((ybins-avg[i])**2,
                                  weights=counts[i,:]))/sqrt(N)
            mask[i] = False
        else:
            mask[i] = True
            n_masked += 1
    if n_masked > 0: # remove empty bins
        navg = empty(len(avg) - n_masked)
        ndev = empty(len(dev) - n_masked)
        nxbins = empty(len(xbins) - n_masked)
        cum = 0
        for i in xrange(0, len(avg)):
            if not mask[i]:
                navg[i-cum] = avg[i]
                ndev[i-cum] = dev[i]
                nxbins[i-cum] = xbins[i]
            else:
                cum += 1
        avg = navg
        dev = ndev
        xbins = nxbins
    return [avg, dev, xbins]
