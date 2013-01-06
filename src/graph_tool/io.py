#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# graph_tool -- a general graph manipulation python module
#
# Copyright (C) 2006-2013 Tiago de Paula Peixoto <tiago@skewed.de>
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

from __future__ import division, absolute_import, print_function

import pickle
import base64
import atexit
from io import BytesIO
from . import libgraph_tool_core

# IStream and OStream need to be tweaked a little to become a real file-like
# object...


def IStream_read(self, n=None):
    if n == None:
        data = ""
        new_data = None
        while new_data != "":
            new_data = self.Read(1)
            data += new_data
        return data
    else:
        return self.Read(n)


def IStream_readline(self, n=None):
    c = None
    line = ""
    while c != "" and c != "\n" and len(line) < n:
        c = self.Read(1)
        line += c
    return line


def OStream_write(self, s):
    data = s
    if not isinstance(data, str) and isinstance(data, bytes):
        data = data.decode('utf-8')
    self.Write(data, len(s))

libgraph_tool_core.IStream.read = IStream_read
libgraph_tool_core.IStream.readline = IStream_readline
libgraph_tool_core.OStream.write = OStream_write


# define and set the pickler/unpickler functions
def pickler(stream, obj):
    sstream = BytesIO()
    pickle.dump(obj, sstream, -1)
    stream.write(base64.b64encode(sstream.getvalue()))

def unpickler(stream):
    data = stream.read().encode('utf-8')
    sstream = BytesIO(base64.b64decode(data))
    return pickle.load(sstream)

libgraph_tool_core.set_pickler(pickler)
libgraph_tool_core.set_unpickler(unpickler)

def clean_picklers():
    libgraph_tool_core.set_pickler(None)
    libgraph_tool_core.set_unpickler(None)

atexit.register(clean_picklers)
