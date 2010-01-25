# SYNOPSIS
#
#   AX_BOOST_REGEX
#
# DESCRIPTION
#
#   This macro checks to see if the Boost.Regex library is installed. It
#   also attempts to guess the currect library name using several attempts.
#   It tries to build the library name using a user supplied name or suffix
#   and then just the raw library.
#
#   If the library is found, HAVE_BOOST_REGEX is defined and
#   BOOST_REGEX_LIB is set to the name of the library.
#
#   This macro calls AC_SUBST(BOOST_REGEX_LIB).
#
#   In order to ensure that the Regex headers are specified on the include
#   path, this macro requires AX_REGEX to be called.
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

AC_DEFUN([AX_BOOST_REGEX],
[AC_CACHE_CHECK(whether the Boost::Regex library is available,
ac_cv_boost_regex,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 CPPFLAGS_SAVE=$CPPFLAGS
 if test "x$REGEX_INCLUDE_DIR" != "x"; then
   CPPFLAGS="-I$REGEX_INCLUDE_DIR $CPPFLAGS"
 fi
 AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[
#include <boost/regex.hpp>
]],
[[boost::regex r(); return 0;]]),
 ac_cv_boost_regex=yes, ac_cv_boost_regex=no)
 AC_LANG_RESTORE
 CPPFLAGS="$CPPFLAGS_SAVE"
])
if test "$ac_cv_boost_regex" = "yes"; then
  AC_DEFINE(HAVE_BOOST_REGEX,,[define if the Boost::Regex library is available])
  ax_regex_lib=boost_regex
  AC_ARG_WITH([boost-regex],AS_HELP_STRING([--with-boost-regex],[specify the boost regex library or suffix to use]),
  [if test "x$with_boost_regex" != "xno"; then
     ax_regex_lib=$with_boost_regex
     ax_boost_regex_lib=boost_regex-$with_boost_regex
   fi])
  for ax_lib in $ax_regex_lib $ax_boost_regex_lib boost_regex boost_regex-mt boost_regex-mt-py2.5; do
    AC_CHECK_LIB($ax_lib, exit, [BOOST_REGEX_LIB=$ax_lib break])
  done
  AC_SUBST(BOOST_REGEX_LIB)
fi
])dnl
