# ===========================================================================
#            http://autoconf-archive.cryp.to/ax_boost_graph.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_BOOST_GRAPH
#
# DESCRIPTION
#
#   This macro checks to see if the Boost.Graph library is installed. It
#   also attempts to guess the currect library name using several attempts.
#   It tries to build the library name using a user supplied name or suffix
#   and then just the raw library.
#
#   If the library is found, HAVE_BOOST_GRAPH is defined and
#   BOOST_GRAPH_LIB is set to the name of the library.
#
#   This macro calls AC_SUBST(BOOST_GRAPH_LIB).
#
#   In order to ensure that the Graph headers are specified on the include
#   path, this macro requires AX_GRAPH to be called.
#
# LAST MODIFICATION
#
#   2008-04-12
#
# COPYLEFT
#
#   Copyright (c) 2008 Michael Tindal
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 2 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Macro Archive. When you make and
#   distribute a modified version of the Autoconf Macro, you may extend this
#   special exception to the GPL to apply to your modified version as well.

AC_DEFUN([AX_BOOST_GRAPH],
[AC_CACHE_CHECK(whether the Boost::Graph library is available,
ac_cv_boost_graph,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 CPPFLAGS_SAVE=$CPPFLAGS
 if test "x$GRAPH_INCLUDE_DIR" != "x"; then
   CPPFLAGS="-I$GRAPH_INCLUDE_DIR $CPPFLAGS"
 fi
 AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[
 #include <boost/graph/adjacency_list.hpp>
 using namespace boost;
 adjacency_list<> g;
 ]],
 [[return 0;]]),
 ac_cv_boost_graph=yes, ac_cv_boost_graph=no)
 AC_LANG_RESTORE
 CPPFLAGS="$CPPFLAGS_SAVE"
])
if test "$ac_cv_boost_graph" = "yes"; then
  AC_DEFINE(HAVE_BOOST_GRAPH,,[define if the Boost::Graph library is available])
  ax_graph_lib=boost_graph
  AC_ARG_WITH([boost-graph],AS_HELP_STRING([--with-boost-graph],[specify the boost graph library or suffix to use]),
  [if test "x$with_boost_graph" != "xno"; then
     ax_graph_lib=$with_boost_graph
     ax_boost_graph_lib=boost_graph-$with_boost_graph
   fi])
  for ax_lib in $ax_graph_lib $ax_boost_graph_lib boost_graph boost_graph-mt boost_graph-mt-py2.5; do
    AC_CHECK_LIB($ax_lib, exit, [BOOST_GRAPH_LIB=$ax_lib break])
  done
  AC_SUBST(BOOST_GRAPH_LIB)
fi
])dnl
