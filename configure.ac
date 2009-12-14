dnl Process this file with autoconf to produce a configure script.

dnl graph-tool package version number
dnl An odd micro number indicates in-progress development.
dnl An even micro number indicates a released version.
m4_define(graph_tool_version_major, 2)
m4_define(graph_tool_version_minor, 2)
m4_define(graph_tool_version_micro, 1)

AC_INIT([graph-tool],
        [graph_tool_version_major().graph_tool_version_minor().graph_tool_version_micro()],
        [http://graph-tool.forked.de])

GRAPH_TOOL_VERSION_MAJOR=graph_tool_version_major()
GRAPH_TOOL_VERSION_MINOR=graph_tool_version_minor()
GRAPH_TOOL_VERSION_MICRO=graph_tool_version_micro()
AC_SUBST(GRAPH_TOOL_VERSION_MAJOR)
AC_SUBST(GRAPH_TOOL_VERSION_MINOR)
AC_SUBST(GRAPH_TOOL_VERSION_MICRO)

AC_CONFIG_SRCDIR(src)
AM_INIT_AUTOMAKE([foreign])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
AM_PROG_CC_C_O
AM_CONFIG_HEADER(config.h)
AC_CONFIG_MACRO_DIR([m4])

AC_ISC_POSIX
AC_SUBST(CXXFLAGS,$CXXFLAGS)
AC_GNU_SOURCE
AC_PROG_CXX
AC_LANG_CPLUSPLUS
AC_PROG_CC_STDC
AC_HEADER_STDC
AC_DISABLE_STATIC
AC_PROG_LIBTOOL
AC_PROG_INSTALL

dnl turn on warnings
[CXXFLAGS="${CXXFLAGS} -Wall"]

dnl Checks for options

dnl Turn debugging on or off
AC_MSG_CHECKING(whether to enable debug info...)

AC_ARG_ENABLE([debug], [AC_HELP_STRING([--enable-debug@<:@=full@:>@],
        [compile debug information [default=yes] ])],
        [AC_MSG_RESULT(yes)]
        [AC_DEFINE([DEBUG], 1, [compile debug info])]
        [CXXFLAGS=" ${CXXFLAGS} -ggdb3 "],
        [AC_MSG_RESULT(no)]
        )

dnl set template depth and hidden visibility

[CXXFLAGS="${CXXFLAGS} -ftemplate-depth-150"]

AC_MSG_CHECKING(whether to enable function inlining...)

AC_ARG_ENABLE([inline], [AC_HELP_STRING([--disable-inline],
        [disable inline functions [default=enabled] ])],
        if test $enableval = no; then
           [AC_MSG_RESULT(no)]
           [AC_DEFINE([NO_INLINE], 1, [disable function inlining])]
           [CXXFLAGS="${CXXFLAGS} -fno-inline -fno-default-inline"]
           NO_INLINE=yes
        else
           [AC_MSG_RESULT(yes)]
        fi
        ,
        [AC_MSG_RESULT(yes)]
        NO_INLINE=no
        )

AM_CONDITIONAL(NO_INLINE, test $NO_INLINE = yes)

AC_MSG_CHECKING(whether to enable optimizations...)

AC_ARG_ENABLE([optimization], [AC_HELP_STRING([--disable-optimization],
        [disable optimization [default=enabled] ])],
        [AC_MSG_RESULT(no)]
        ,
        [AC_MSG_RESULT(yes)]
        [CXXFLAGS="${CXXFLAGS} -O99"]
        )

AC_MSG_CHECKING(whether to enable symbol visibility manipulation...)

AC_ARG_ENABLE([visibility], [AC_HELP_STRING([--disable-visibility],
        [disable symbol visibility manipulation [default=enabled] ])],
        [AC_MSG_RESULT(no)]
        ,
        [AC_MSG_RESULT(yes)]
        [CXXFLAGS="${CXXFLAGS} -fvisibility=default -fvisibility-inlines-hidden"]
        )


AC_MSG_CHECKING(whether to enable graph filtering...)
AC_ARG_ENABLE([graph-filtering], [AC_HELP_STRING([--disable-graph-filtering],
        [disable graph filtering [default=enabled] ])],
        if test $enableval = no; then
           [AC_MSG_RESULT(no)]
           [AC_DEFINE([NO_GRAPH_FILTERING], 1, [disable graph filtering])]
           NO_GRAPH_FILTERING=yes
        else
           [AC_MSG_RESULT(yes)]
        fi
        ,
        [AC_MSG_RESULT(yes)]
        NO_GRAPH_FILTERING=no
        )

AC_MSG_CHECKING(whether to enable parallel algorithms with openmp...)
AC_ARG_ENABLE([openmp], [AC_HELP_STRING([--enable-openmp],
        [enable openmp [default=disabled] ])],
        if test $enableval = yes; then
           [AC_MSG_RESULT(yes)]
           [AC_DEFINE([USING_OPENMP], 1, [using openmp])]
           USING_OPENMP=yes
           [CXXFLAGS="${CXXFLAGS} -fopenmp"]
           [OPENMP_LDFLAGS=" -lgomp "]
        fi
        ,
        [AC_MSG_RESULT(no)]
        [AC_DEFINE([USING_OPENMP], 0, [using openmp])]
        USING_OPENMP=no
        [CXXFLAGS="${CXXFLAGS} -Wno-unknown-pragmas"]
        [OPENMP_LDFLAGS=""]
        )
AC_SUBST(OPENMP_LDFLAGS)

dnl Override prefix with default value if not set
if test "x${prefix}" = "xNONE"; then
  prefix=${ac_default_prefix}
fi

dnl Checks for programs.

dnl Checks for libraries.

dnl math library
AC_CHECK_LIB(m,main)

dnl bzip2 compression library
AC_CHECK_LIB(bz2,main)

dnl expat
AC_CHECK_LIB(expat,main)

dnl python
[TEMP_CPPFLAGS="${CPPFLAGS}"]
AM_PATH_PYTHON([2.5])
AC_PYTHON_DEVEL([2.5])
[CPPFLAGS="${CPPFLAGS} ${PYTHON_CPPFLAGS}"]
[LDFLAGS="${LDFLAGS} ${PYTHON_LDFLAGS}"]

dnl boost
AX_BOOST_BASE([1.40.0])
AX_BOOST_PYTHON
if test "$BOOST_PYTHON_LIB" = ""; then
   AC_MSG_ERROR([No usable boost::python found])
fi
AX_BOOST_IOSTREAMS
if test "$BOOST_IOSTREAMS_LIB" = ""; then
   AC_MSG_ERROR([No usable boost::iostreams found])
fi
AX_BOOST_REGEX
if test "$BOOST_REGEX_LIB" = ""; then
   AC_MSG_ERROR([No usable boost::regex found])
fi
AX_BOOST_GRAPH
if test "$BOOST_GRAPH_LIB" = ""; then
   AC_MSG_ERROR([No usable boost::graph found])
fi
[CPPFLAGS="${CPPFLAGS} ${BOOST_CPPFLAGS}"]

dnl CGAL
AC_CHECK_LIB(CGAL,main)
[CGAL_LIBADD="-lCGAL"]
AC_SUBST(CGAL_LIBADD)
[CGAL_FLAGS="-frounding-math"]
AC_SUBST(CGAL_FLAGS)

dnl Checks for header files.

dnl numpy
AC_PYTHON_MODULE(numpy, fatal)
NUMPY_DIR="${PYTHON_SITE_PKG}/numpy"
AC_ARG_WITH([numpy], [AC_HELP_STRING([--with-numpy=DIR],
                     [directory where numpy is installed
                      [default=${PYTHON_SITE_PKG}/numpy] ])],
        NUMPY_DIR=$withval)
[NUMPY_DIR="${NUMPY_DIR}/core/include/numpy"]
AC_CHECK_HEADER([${NUMPY_DIR}/arrayobject.h],
                [CPPFLAGS="${CPPFLAGS} -I${NUMPY_DIR}"],
                [AC_MSG_ERROR([Numpy extension header not found])],
                [#include "Python.h"])

dnl scipy
AC_PYTHON_MODULE(scipy, fatal)
[SCIPY_DIR="${PYTHON_SITE_PKG}/scipy"]
AC_ARG_WITH([scipy], [AC_HELP_STRING([--with-scipy=DIR],
                     [scipy install directory
                      [default=${PYTHON_SITE_PKG}/scipy] ])],
            [SCIPY_DIR=$withval])
AC_CHECK_HEADER([${SCIPY_DIR}/weave/scxx/object.h],
                [AC_DEFINE([HAVE_SCIPY], [1], [using scipy's weave])]
                [CPPFLAGS="${CPPFLAGS} -I${SCIPY_DIR}"],
                [AC_MSG_ERROR([Scipy weave/scxx header not found])])

dnl expat
AC_CHECK_HEADER([expat.h],
                [],
                [AC_MSG_ERROR([Expat header not found])])

dnl cgal
AC_CHECK_HEADER([CGAL/version.h],
                [],
                [AC_MSG_ERROR([CGAL headers not found])])

dnl Checks for typedefs, structures, and compiler characteristics.
dnl Checks for library functions.

dnl Checks for Additional stuff.

dnl Set PACKAGE SOURCE DIR in config.h.
packagesrcdir=`cd $srcdir && pwd`

dnl Set PACKAGE PREFIX
if test "x${prefix}" = "xNONE"; then
  packageprefix=${ac_default_prefix}
else
  packageprefix=${prefix}
fi

dnl Set PACKAGE DATA & DOC DIR
packagedatadir=share/${PACKAGE}
packagedocdir=doc/${PACKAGE}


dnl Subst PACKAGE_DATA_DIR.
NO_PREFIX_PACKAGE_DATA_DIR="${packagedatadir}"
AC_SUBST(NO_PREFIX_PACKAGE_DATA_DIR)
PACKAGE_DATA_DIR="${packageprefix}/${packagedatadir}"
AC_SUBST(PACKAGE_DATA_DIR)

dnl Subst PACKAGE_DOC_DIR.
NO_PREFIX_PACKAGE_DOC_DIR="${packagedocdir}"
AC_SUBST(NO_PREFIX_PACKAGE_DOC_DIR)
PACKAGE_DOC_DIR="${packageprefix}/${packagedocdir}"
AC_SUBST(PACKAGE_DOC_DIR)

dnl Subst PYTHON_DIR.
AC_DEFINE_UNQUOTED([INSTALL_PREFIX],"${prefix}", [python prefix])
AC_DEFINE_UNQUOTED([PYTHON_DIR], "${PYTHON_SITE_PKG}", [The directory name for the site-packages subdirectory of the standard Python install tree.])

AC_DEFINE_UNQUOTED([CXXFLAGS],"${CXXFLAGS}", [c++ compilation options])
AC_DEFINE_UNQUOTED([CPPFLAGS],"${CPPFLAGS}", [c++ preprocessor compilation options])
AC_DEFINE_UNQUOTED([LDFLAGS],"${LDFLAGS}", [linker options])

AC_DEFINE_UNQUOTED([PACKAGE_DATA_DIR], "${packageprefix}/${packagedatadir}", [package data dir])
AC_DEFINE_UNQUOTED([PACKAGE_DOC_DIR], "${packageprefix}/${packagedocdir}", [package doc dir])
AC_DEFINE_UNQUOTED([PACKAGE_SOURCE_DIR], "${packagesrcdir}", [package source dir])

AC_DEFINE_UNQUOTED([AUTHOR], "Tiago de Paula Peixoto <tiago@forked.de>", [program author(s)])

AC_DEFINE_UNQUOTED([COPYRIGHT],
    "Copyright (C) 2008 Tiago de Paula Peixoto\nThis is free software; see the source for copying conditions.  There is NO\nwarranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.",
	[copyright info])

# git information
AC_DEFINE([GIT_COMMIT], "esyscmd(git show | head -n 1 | sed 's/commit //' |  grep -o -e '.\{8\}' | head -n 1 |tr -d '\n')", [git HEAD commit hash])
AC_DEFINE([GIT_COMMIT_DATE], "esyscmd(git log -1 | head -n 3 | grep 'Date:' | sed s/'Date:   '// | tr -d '\n')", [git HEAD commit date])

# set some automake variables

# global CPP flags
AC_SUBST(CPPFLAGS)

# global CXX flags
AC_SUBST(CXXFLAGS)

# CPP flags for submodules
[MOD_CPPFLAGS="-I\$(srcdir) -I\$(srcdir)/.. -I\$(srcdir)/../../boost-workaround \
-DHAVE_CONFIG_H"]
AC_SUBST(MOD_CPPFLAGS)

# default directory for submodules
[MOD_DIR="${pythondir}/graph_tool"]
AC_SUBST(MOD_DIR)

# default LIBADD flags for submodules
[MOD_LIBADD="${PYTHON_LDFLAGS} ${BOOST_IOSTREAMS_LIB} -l${BOOST_PYTHON_LIB} \
${BOOST_REGEX_LIB} -l${BOOST_GRAPH_LIB} ${OPENMP_LDFLAGS} -lexpat"]
AC_SUBST(MOD_LIBADD)

# needed for typeinfo objects to work across DSO boundaries.
# see http://gcc.gnu.org/faq.html#dso
[MOD_LDFLAGS="-module  -avoid-version  -export-dynamic  -no-undefined\
-Wl,-E"]
AC_SUBST(MOD_LDFLAGS)

AC_OUTPUT([
Makefile
src/Makefile
src/graph/Makefile
src/graph/correlations/Makefile
src/graph/generation/Makefile
src/graph/stats/Makefile
src/graph/centrality/Makefile
src/graph/clustering/Makefile
src/graph/community/Makefile
src/graph/util/Makefile
src/graph/topology/Makefile
src/graph/flow/Makefile
src/graph/layout/Makefile
src/graph_tool/Makefile
])

# ugly, dirty, smelly hack to force libtool to link to libstc++ before anything
# else, to avoid strange exception-related segfaults in
# __cxa_allocate_exception() see:
# http://www.nabble.com/preventing-boost::python::throw_error_already_set-td20958723.html
# http://wiki.fifengine.de/Segfault_in_cxa_allocate_exception
AC_MSG_NOTICE([Modifying libtool to correctly link libstdc++ before anything else...])
libtool=`cat libtool | sed "s/-nostdlib/-nostdlib -lstdc++/g"`
echo "$libtool" > libtool
