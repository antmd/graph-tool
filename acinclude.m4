dnl @synopsis AX_PYTHON
dnl
dnl This macro does a complete Python development environment check.
dnl
dnl It recurses through several python versions (from 2.1 to 2.4 in
dnl this version), looking for an executable. When it finds an
dnl executable, it looks to find the header files and library.
dnl
dnl It sets PYTHON_BIN to the name of the python executable,
dnl PYTHON_INCLUDE_DIR to the directory holding the header files, and
dnl PYTHON_LIB to the name of the Python library.
dnl
dnl This macro calls AC_SUBST on PYTHON_BIN (via AC_CHECK_PROG),
dnl PYTHON_INCLUDE_DIR and PYTHON_LIB.
dnl
dnl @category InstalledPackages
dnl @author Michael Tindal <mtindal@paradoxpoint.com>
dnl @version 2004-09-20
dnl @license GPLWithACException

AC_DEFUN([AX_PYTHON],
[AC_MSG_CHECKING(for python build information)
AC_MSG_RESULT([])
for python in python2.4 python2.3 python2.2 python2.1 python; do
AC_CHECK_PROGS(PYTHON_BIN, [$python])
ax_python_bin=$PYTHON_BIN
if test x$ax_python_bin != x; then
   AC_CHECK_LIB($ax_python_bin, main, ax_python_lib=$ax_python_bin, ax_python_lib=no)
   AC_CHECK_HEADER([$ax_python_bin/Python.h],
   [[ax_python_header=`locate $ax_python_bin/Python.h | sed -e s,/Python.h,,`]],
   ax_python_header=no)
   if test $ax_python_lib != no; then
     if test $ax_python_header != no; then
       break;
     fi
   fi
fi
done
if test x$ax_python_bin = x; then
   ax_python_bin=no
fi
if test x$ax_python_header = x; then
   ax_python_header=no
fi
if test x$ax_python_lib = x; then
   ax_python_lib=no
fi

AC_MSG_RESULT([  results of the Python check:])
AC_MSG_RESULT([    Binary:      $ax_python_bin])
AC_MSG_RESULT([    Library:     $ax_python_lib])
AC_MSG_RESULT([    Include Dir: $ax_python_header])

if test x$ax_python_header != xno; then
  PYTHON_INCLUDE_DIR=$ax_python_header
  AC_SUBST(PYTHON_INCLUDE_DIR)
fi
if test x$ax_python_lib != xno; then
  PYTHON_LIB=$ax_python_lib
  AC_SUBST(PYTHON_LIB)
fi
])dnl
AC_DEFUN([AM_GNU_GETTEXT_VERSION], [])
dnl @synopsis AX_BOOST([MINIMUM-VERSION], [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
dnl
dnl Test for the Boost C++ libraries of a particular version (or newer)
dnl
dnl If no path to the installed boost library is given the macro
dnl searchs under /usr, /usr/local, and /opt, and evaluates the
dnl $BOOST_ROOT environment variable. Further documentation is
dnl available at <http://randspringer.de/boost/index.html>.
dnl
dnl This macro calls:
dnl
dnl   AC_SUBST(BOOST_CPPFLAGS) / AC_SUBST(BOOST_LDFLAGS)
dnl   AC_SUBST(BOOST_FILESYSTEM_LIB)
dnl   AC_SUBST(BOOST_PROGRAM_OPTIONS_LIB)
dnl   AC_SUBST(BOOST_THREAD_LIB)
dnl   AC_SUBST(BOOST_IOSTREAMS_LIB)
dnl   AC_SUBST(BOOST_SERIALIZATION_LIB)
dnl   AC_SUBST(BOOST_WSERIALIZATION_LIB)
dnl   AC_SUBST(BOOST_SIGNALS_LIB)
dnl   AC_SUBST(BOOST_DATE_TIME_LIB)
dnl   AC_SUBST(BOOST_REGEX_LIB)
dnl   AC_SUBST(BOOST_UNIT_TEST_FRAMEWORK_LIB)
dnl
dnl And sets:
dnl
dnl   HAVE_BOOST
dnl   HAVE_BOOST_FILESYSTEM
dnl   HAVE_BOOST_PROGRAM_OPTIONS
dnl   HAVE_BOOST_THREAD
dnl   HAVE_BOOST_IOSTREAMS
dnl   HAVE_BOOST_SERIALIZATION
dnl   HAVE_BOOST_SIGNALS
dnl   HAVE_BOOST_DATE_TIME
dnl   HAVE_BOOST_REGEX
dnl   HAVE_BOOST_UNIT_TEST_FRAMEWORK
dnl
dnl @category InstalledPackages
dnl @category Cxx
dnl @author Thomas Porschberg <thomas@randspringer.de>
dnl @version 2006-06-15
dnl @license AllPermissive

AC_DEFUN([AX_BOOST],
[
    AC_ARG_WITH([boost],
                AS_HELP_STRING([--with-boost=DIR],
                [use boost (default is YES) specify the root directory for boost library (optional)]),
                [
                if test "$withval" = "no"; then
		            want_boost="no"
                elif test "$withval" = "yes"; then
                    want_boost="yes"
                    ac_boost_path=""
                else
			        want_boost="yes"
            		ac_boost_path="$withval"
		        fi
            	],
                [want_boost="yes"])

    AC_CANONICAL_BUILD
	if test "x$want_boost" = "xyes"; then
        AC_REQUIRE([AC_PROG_CC])
		boost_lib_version_req=ifelse([$1], ,1.20.0,$1)
		boost_lib_version_req_shorten=`expr $boost_lib_version_req : '\([[0-9]]*\.[[0-9]]*\)'`
		boost_lib_version_req_major=`expr $boost_lib_version_req : '\([[0-9]]*\)'`
		boost_lib_version_req_minor=`expr $boost_lib_version_req : '[[0-9]]*\.\([[0-9]]*\)'`
		boost_lib_version_req_sub_minor=`expr $boost_lib_version_req : '[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
		if test "x$boost_lib_version_req_sub_minor" = "x" ; then
			boost_lib_version_req_sub_minor="0"
    	fi
		WANT_BOOST_VERSION=`expr $boost_lib_version_req_major \* 100000 \+  $boost_lib_version_req_minor \* 100 \+ $boost_lib_version_req_sub_minor`
		AC_MSG_CHECKING(for boostlib >= $boost_lib_version_req)
		succeeded=no

		dnl first we check the system location for boost libraries
		dnl this location ist chosen if boost libraries are installed with the --layout=system option
		dnl or if you install boost with RPM
		if test "$ac_boost_path" != ""; then
			BOOST_LDFLAGS="-L$ac_boost_path/lib"
			BOOST_CPPFLAGS="-I$ac_boost_path/include"
		else
			for ac_boost_path_tmp in /usr /usr/local /opt ; do
				if test -d "$ac_boost_path_tmp/include/boost" && test -r "$ac_boost_path_tmp/include/boost"; then
					BOOST_LDFLAGS="-L$ac_boost_path_tmp/lib"
					BOOST_CPPFLAGS="-I$ac_boost_path_tmp/include"
					break;
				fi
			done
		fi

		CPPFLAGS_SAVED="$CPPFLAGS"
		CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
		export CPPFLAGS

		LDFLAGS_SAVED="$LDFLAGS"
		LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
		export LDFLAGS

	AC_LANG_PUSH(C++)
     	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
@%:@include <boost/version.hpp>
]],
       [[
#if BOOST_VERSION >= $WANT_BOOST_VERSION
// Everything is okay
#else
#  error Boost version is too old
#endif

		]])],
    	[
         AC_MSG_RESULT(yes)
		 succeeded=yes
		 found_system=yes
         ifelse([$2], , :, [$2])
       ],
       [
       ])
       AC_LANG_POP([C++])
		dnl if we found no boost with system layout we search for boost libraries
		dnl built and installed without the --layout=system option or for a staged(not installed) version
		if test "x$succeeded" != "xyes"; then
			_version=0
			if test "$ac_boost_path" != ""; then
                BOOST_LDFLAGS="-L$ac_boost_path/lib"
				if test -d "$ac_boost_path" && test -r "$ac_boost_path"; then
					for i in `ls -d $ac_boost_path/include/boost-* 2>/dev/null`; do
						_version_tmp=`echo $i | sed "s#$ac_boost_path##" | sed 's/\/include\/boost-//' | sed 's/_/./'`
						V_CHECK=`expr $_version_tmp \> $_version`
						if test "$V_CHECK" = "1" ; then
							_version=$_version_tmp
						fi
						VERSION_UNDERSCORE=`echo $_version | sed 's/\./_/'`
						BOOST_CPPFLAGS="-I$ac_boost_path/include/boost-$VERSION_UNDERSCORE"
					done
				fi
			else
				for ac_boost_path in /usr /usr/local /opt ; do
					if test -d "$ac_boost_path" && test -r "$ac_boost_path"; then
						for i in `ls -d $ac_boost_path/include/boost-* 2>/dev/null`; do
							_version_tmp=`echo $i | sed "s#$ac_boost_path##" | sed 's/\/include\/boost-//' | sed 's/_/./'`
							V_CHECK=`expr $_version_tmp \> $_version`
							if test "$V_CHECK" = "1" ; then
								_version=$_version_tmp
								best_path=$ac_boost_path
							fi
						done
					fi
				done

				VERSION_UNDERSCORE=`echo $_version | sed 's/\./_/'`
				BOOST_CPPFLAGS="-I$best_path/include/boost-$VERSION_UNDERSCORE"
				BOOST_LDFLAGS="-L$best_path/lib"

	    		if test "x$BOOST_ROOT" != "x"; then
                    if test -d "$BOOST_ROOT" && test -r "$BOOST_ROOT" && test -d "$BOOST_ROOT/stage/lib" && test -r "$BOOST_ROOT/stage/lib"; then
						version_dir=`expr //$BOOST_ROOT : '.*/\(.*\)'`
						stage_version=`echo $version_dir | sed 's/boost_//' | sed 's/_/./g'`
						stage_version_shorten=`expr $stage_version : '\([[0-9]]*\.[[0-9]]*\)'`
						V_CHECK=`expr $stage_version_shorten \>\= $_version`
						if test "$V_CHECK" = "1" ; then
							AC_MSG_NOTICE(We will use a staged boost library from $BOOST_ROOT)
							BOOST_CPPFLAGS="-I$BOOST_ROOT"
							BOOST_LDFLAGS="-L$BOOST_ROOT/stage/lib"
						fi
					fi
	    		fi
			fi

			CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
			export CPPFLAGS
			LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
			export LDFLAGS

            AC_LANG_PUSH(C++)
            AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
@%:@include <boost/version.hpp>
]],
       [[
#if BOOST_VERSION >= $WANT_BOOST_VERSION
// Everything is okay
#else
#  error Boost version is too old
#endif

		]])],
    	[
         AC_MSG_RESULT(yes ($_version))
		 succeeded=yes
         ifelse([$2], , :, [$2])
       ],
       [
         AC_MSG_RESULT(no ($_version))
         ifelse([$3], , :, [$3])
       ])
    	AC_LANG_POP([C++])
		fi

		if test "$succeeded" != "yes" ; then
			if test "$_version" = "0" ; then
				AC_MSG_ERROR([[We could not detect the boost libraries (version $boost_lib_version_req_shorten or higher). If you have a staged boost library (still not installed) please specify \$BOOST_ROOT in your environment and do not give a PATH to --with-boost option.  If you are sure you have boost installed, then check your version number looking in <boost/version.hpp>. See http://randspringer.de/boost for more documentation.]])
			else
				AC_MSG_ERROR('Your boost libraries seems to old (version $_version).  We need at least $boost_lib_version_shorten')
			fi
		else
			AC_SUBST(BOOST_CPPFLAGS)
			AC_SUBST(BOOST_LDFLAGS)
			AC_DEFINE(HAVE_BOOST,,[define if the Boost library is available])

			AC_CACHE_CHECK([whether the Boost::Filesystem library is available],
						   ax_cv_boost_filesystem,
						[AC_LANG_PUSH([C++])
			 AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[@%:@include <boost/filesystem/path.hpp>]],
                                   [[using namespace boost::filesystem;
                                   path my_path( "foo/bar/data.txt" );
                                   return 0;]]),
            				       ax_cv_boost_filesystem=yes, ax_cv_boost_filesystem=no)
                                   AC_LANG_POP([C++])
			])
			if test "$ax_cv_boost_filesystem" = "yes"; then
				AC_DEFINE(HAVE_BOOST_FILESYSTEM,,[define if the Boost::FILESYSTEM library is available])
				BN=boost_filesystem
				for ax_lib in $BN $BN-$CC $BN-$CC-mt $BN-$CC-mt-s $BN-$CC-s \
                              lib$BN lib$BN-$CC lib$BN-$CC-mt lib$BN-$CC-mt-s lib$BN-$CC-s \
                              $BN-mgw $BN-mgw $BN-mgw-mt $BN-mgw-mt-s $BN-mgw-s ; do
				    AC_CHECK_LIB($ax_lib, main,
                                 [BOOST_FILESYSTEM_LIB="-l$ax_lib" AC_SUBST(BOOST_FILESYSTEM_LIB) link_filesystem="yes" break],
                                 [link_filesystem="no"])
  				done
				if test "x$link_filesystem" = "xno"; then
					AC_MSG_NOTICE(Could not link against $ax_lib !)
				fi
			fi

			AC_CACHE_CHECK([whether the Boost::Program_Options library is available],
						   ax_cv_boost_program_options,
						   [AC_LANG_PUSH([C++])
			               AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[@%:@include <boost/program_options.hpp>]],
                                   [[boost::program_options::options_description generic("Generic options");
                                   return 0;]]),
                           ax_cv_boost_program_options=yes, ax_cv_boost_program_options=no)
                           AC_LANG_POP([C++])
			])
			if test "$ax_cv_boost_program_options" = yes; then
				AC_DEFINE(HAVE_BOOST_PROGRAM_OPTIONS,,[define if the Boost::PROGRAM_OPTIONS library is available])
				BN=boost_program_options
				for ax_lib in $BN $BN-$CC $BN-$CC-mt $BN-$CC-mt-s $BN-$CC-s \
                              lib$BN lib$BN-$CC lib$BN-$CC-mt lib$BN-$CC-mt-s lib$BN-$CC-s \
                              $BN-mgw $BN-mgw $BN-mgw-mt $BN-mgw-mt-s $BN-mgw-s ; do
				    AC_CHECK_LIB($ax_lib, main,
                                 [BOOST_PROGRAM_OPTIONS_LIB="-l$ax_lib" AC_SUBST(BOOST_PROGRAM_OPTIONS_LIB) link_program_options="yes" break],
                                 [link_program_options="no"])
  				done
				if test "x$link_program_options="no"" = "xno"; then
					AC_MSG_NOTICE(Could not link against $ax_lib !)
				fi
			fi

			AC_CACHE_CHECK(whether the Boost::Thread library is available,
						   ax_cv_boost_thread,
						[AC_LANG_PUSH([C++])
			 CXXFLAGS_SAVE=$CXXFLAGS

			 if test "x$build_os" = "xsolaris" ; then
  				 CXXFLAGS="-pthreads $CXXFLAGS"
			 elif test "x$build_os" = "xming32" ; then
				 CXXFLAGS="-mthreads $CXXFLAGS"
			 else
				CXXFLAGS="-pthread $CXXFLAGS"
			 fi
			 AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[@%:@include <boost/thread/thread.hpp>]],
                                   [[boost::thread_group thrds;
                                   return 0;]]),
                   ax_cv_boost_thread=yes, ax_cv_boost_thread=no)
			 CXXFLAGS=$CXXFLAGS_SAVE
             AC_LANG_POP([C++])
			])
			if test "x$ax_cv_boost_thread" = "xyes"; then
               if test "x$build_os" = "xsolaris" ; then
 				  BOOST_CPPFLAGS="-pthreads $BOOST_CPPFLAGS"
			   elif test "x$build_os" = "xming32" ; then
 				  BOOST_CPPFLAGS="-mthreads $BOOST_CPPFLAGS"
			   else
				  BOOST_CPPFLAGS="-pthread $BOOST_CPPFLAGS"
			   fi

				AC_SUBST(BOOST_CPPFLAGS)
				AC_DEFINE(HAVE_BOOST_THREAD,,[define if the Boost::THREAD library is available])
				BN=boost_thread
				for ax_lib in $BN $BN-$CC $BN-$CC-mt $BN-$CC-mt-s $BN-$CC-s \
                              lib$BN lib$BN-$CC lib$BN-$CC-mt lib$BN-$CC-mt-s lib$BN-$CC-s \
                              $BN-mgw $BN-mgw $BN-mgw-mt $BN-mgw-mt-s $BN-mgw-s ; do
				    AC_CHECK_LIB($ax_lib, main, [BOOST_THREAD_LIB="-l$ax_lib" AC_SUBST(BOOST_THREAD_LIB) link_thread="yes" break],
                                 [link_thread="no"])
  				done
				if test "x$link_thread" = "xno"; then
					AC_MSG_NOTICE(Could not link against $ax_lib !)
				fi
			fi

			AC_CACHE_CHECK(whether the Boost::IOStreams library is available,
						   ax_cv_boost_iostreams,
						[AC_LANG_PUSH([C++])
			 AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[@%:@include <boost/iostreams/filtering_stream.hpp>
												 @%:@include <boost/range/iterator_range.hpp>
												]],
                                   [[std::string  input = "Hello World!";
									 namespace io = boost::iostreams;
									 io::filtering_istream  in(boost::make_iterator_range(input));
									 return 0;
                                   ]]),
                   ax_cv_boost_iostreams=yes, ax_cv_boost_iostreams=no)
			 AC_LANG_POP([C++])
			])
			if test "x$ax_cv_boost_iostreams" = "xyes"; then
				AC_DEFINE(HAVE_BOOST_IOSTREAMS,,[define if the Boost::IOStreams library is available])
				BN=boost_iostreams
				for ax_lib in $BN $BN-$CC $BN-$CC-mt $BN-$CC-mt-s $BN-$CC-s \
                              lib$BN lib$BN-$CC lib$BN-$CC-mt lib$BN-$CC-mt-s lib$BN-$CC-s \
                              $BN-mgw $BN-mgw $BN-mgw-mt $BN-mgw-mt-s $BN-mgw-s ; do
				    AC_CHECK_LIB($ax_lib, main, [BOOST_IOSTREAMS_LIB="-l$ax_lib" AC_SUBST(BOOST_IOSTREAMS_LIB) link_thread="yes" break],
                                 [link_thread="no"])
  				done
				if test "x$link_thread" = "xno"; then
					AC_MSG_NOTICE(Could not link against $ax_lib !)
				fi
			fi

			AC_CACHE_CHECK(whether the Boost::Serialization library is available,
						   ax_cv_boost_serialization,
						[AC_LANG_PUSH([C++])
			 AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[@%:@include <fstream>
												 @%:@include <boost/archive/text_oarchive.hpp>
                                                 @%:@include <boost/archive/text_iarchive.hpp>
												]],
                                   [[std::ofstream ofs("filename");
									boost::archive::text_oarchive oa(ofs);
									 return 0;
                                   ]]),
                   ax_cv_boost_serialization=yes, ax_cv_boost_serialization=no)
			 AC_LANG_POP([C++])
			])
			if test "x$ax_cv_boost_serialization" = "xyes"; then
				AC_DEFINE(HAVE_BOOST_SERIALIZATION,,[define if the Boost::Serialization library is available])
				BN=boost_serialization
				for ax_lib in $BN $BN-$CC $BN-$CC-mt $BN-$CC-mt-s $BN-$CC-s \
                              lib$BN lib$BN-$CC lib$BN-$CC-mt lib$BN-$CC-mt-s lib$BN-$CC-s \
                              $BN-mgw $BN-mgw $BN-mgw-mt $BN-mgw-mt-s $BN-mgw-s ; do
				    AC_CHECK_LIB($ax_lib, main,
                                 [BOOST_SERIALIZATION_LIB="-l$ax_lib" AC_SUBST(BOOST_SERIALIZATION_LIB) link_serialization="yes" break],
                                 [link_serialization="no"])
  				done
				if test "x$link_serialization" = "xno"; then
					AC_MSG_NOTICE(Could not link against $ax_lib !)
				fi

				BN=boost_wserialization
				for ax_lib in $BN $BN-$CC $BN-$CC-mt $BN-$CC-mt-s $BN-$CC-s \
                              lib$BN lib$BN-$CC lib$BN-$CC-mt lib$BN-$CC-mt-s lib$BN-$CC-s \
                              $BN-mgw $BN-mgw $BN-mgw-mt $BN-mgw-mt-s $BN-mgw-s ; do
				    AC_CHECK_LIB($ax_lib, main,
                                 [BOOST_WSERIALIZATION_LIB="-l$ax_lib" AC_SUBST(BOOST_WSERIALIZATION_LIB) link_wserialization="yes" break],
                                 [link_wserialization="no"])
  				done
				if test "x$link_wserialization" = "xno"; then
					AC_MSG_NOTICE(Could not link against $ax_lib !)
				fi
			fi

			AC_CACHE_CHECK(whether the Boost::Signals library is available,
						   ax_cv_boost_signals,
						[AC_LANG_PUSH([C++])
			 AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[@%:@include <boost/signal.hpp>
												]],
                                   [[boost::signal<void ()> sig;
                                     return 0;
                                   ]]),
                   ax_cv_boost_signals=yes, ax_cv_boost_signals=no)
			 AC_LANG_POP([C++])
			])
			if test "x$ax_cv_boost_signals" = "xyes"; then
				AC_DEFINE(HAVE_BOOST_SIGNALS,,[define if the Boost::Signals library is available])
				BN=boost_signals
				for ax_lib in $BN $BN-$CC $BN-$CC-mt $BN-$CC-mt-s $BN-$CC-s \
                              lib$BN lib$BN-$CC lib$BN-$CC-mt lib$BN-$CC-mt-s lib$BN-$CC-s \
                              $BN-mgw $BN-mgw $BN-mgw-mt $BN-mgw-mt-s $BN-mgw-s ; do
				    AC_CHECK_LIB($ax_lib, main, [BOOST_SIGNALS_LIB="-l$ax_lib" AC_SUBST(BOOST_SIGNALS_LIB) link_signals="yes" break],
                                 [link_signals="no"])
  				done
				if test "x$link_signals" = "xno"; then
					AC_MSG_NOTICE(Could not link against $ax_lib !)
				fi
			fi

			AC_CACHE_CHECK(whether the Boost::Date_Time library is available,
						   ax_cv_boost_date_time,
						[AC_LANG_PUSH([C++])
			 AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[@%:@include <boost/date_time/gregorian/gregorian_types.hpp>
												]],
                                   [[using namespace boost::gregorian; date d(2002,Jan,10);
                                     return 0;
                                   ]]),
                   ax_cv_boost_date_time=yes, ax_cv_boost_date_time=no)
			 AC_LANG_POP([C++])
			])
			if test "x$ax_cv_boost_date_time" = "xyes"; then
				AC_DEFINE(HAVE_BOOST_DATE_TIME,,[define if the Boost::Date_Time library is available])
				BN=boost_date_time
				for ax_lib in $BN $BN-$CC $BN-$CC-mt $BN-$CC-mt-s $BN-$CC-s \
                              lib$BN lib$BN-$CC lib$BN-$CC-mt lib$BN-$CC-mt-s lib$BN-$CC-s \
                              $BN-mgw $BN-mgw $BN-mgw-mt $BN-mgw-mt-s $BN-mgw-s ; do
				    AC_CHECK_LIB($ax_lib, main, [BOOST_DATE_TIME_LIB="-l$ax_lib" AC_SUBST(BOOST_DATE_TIME_LIB) link_thread="yes" break],
                                 [link_thread="no"])
  				done
				if test "x$link_thread"="no" = "xno"; then
					AC_MSG_NOTICE(Could not link against $ax_lib !)
				fi
			fi

			AC_CACHE_CHECK(whether the Boost::Regex library is available,
						   ax_cv_boost_regex,
						[AC_LANG_PUSH([C++])
			 AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[@%:@include <boost/regex.hpp>
												]],
                                   [[boost::regex r(); return 0;]]),
                   ax_cv_boost_regex=yes, ax_cv_boost_regex=no)
			 AC_LANG_POP([C++])
			])
			if test "x$ax_cv_boost_regex" = "xyes"; then
				AC_DEFINE(HAVE_BOOST_REGEX,,[define if the Boost::Regex library is available])
				BN=boost_regex
				for ax_lib in $BN $BN-$CC $BN-$CC-mt $BN-$CC-mt-s $BN-$CC-s \
                              lib$BN lib$BN-$CC lib$BN-$CC-mt lib$BN-$CC-mt-s lib$BN-$CC-s \
                              $BN-mgw $BN-mgw $BN-mgw-mt $BN-mgw-mt-s $BN-mgw-s ; do
				    AC_CHECK_LIB($ax_lib, main, [BOOST_REGEX_LIB="-l$ax_lib" AC_SUBST(BOOST_REGEX_LIB) link_regex="yes" break],
                                 [link_regex="no"])
  				done
				if test "x$link_regex" = "xno"; then
					AC_MSG_NOTICE(Could not link against $ax_lib !)
				fi
			fi

			AC_CACHE_CHECK(whether the Boost::UnitTestFramework library is available,
						   ax_cv_boost_unit_test_framework,
						[AC_LANG_PUSH([C++])
			 AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[@%:@include <boost/test/unit_test.hpp>]],
                                    [[using boost::unit_test::test_suite;
					                 test_suite* test= BOOST_TEST_SUITE( "Unit test example 1" ); return 0;]]),
                   ax_cv_boost_unit_test_framework=yes, ax_cv_boost_unit_test_framework=no)
			 AC_LANG_POP([C++])
			])
			if test "x$ax_cv_boost_unit_test_framework" = "xyes"; then
    		AC_DEFINE(HAVE_BOOST_UNIT_TEST_FRAMEWORK,,[define if the Boost::Unit_test_framework library is available])
			BN=boost_unit_test_framework
    		saved_ldflags="${LDFLAGS}"
			for ax_lib in $BN $BN-$CC $BN-$CC-mt $BN-$CC-mt-s $BN-$CC-s \
                          lib$BN lib$BN-$CC lib$BN-$CC-mt lib$BN-$CC-mt-s lib$BN-$CC-s \
                          $BN-mgw $BN-mgw $BN-mgw-mt $BN-mgw-mt-s $BN-mgw-s ; do
                LDFLAGS="${LDFLAGS} -l$ax_lib"
    			AC_CACHE_CHECK(the name of the Boost::UnitTestFramework library,
	      					   ax_cv_boost_unit_test_framework_link,
						[AC_LANG_PUSH([C++])
                   AC_LINK_IFELSE([AC_LANG_PROGRAM([[@%:@include <boost/test/unit_test.hpp>
                                                     using boost::unit_test::test_suite;
                                                     test_suite* init_unit_test_suite( int argc, char * argv[] ) {
                                                     test_suite* test= BOOST_TEST_SUITE( "Unit test example 1" );
                                                     return test;
                                                     }
                                                   ]],
                                 [[ return 0;]])],
                                 link_unit_test_framework="yes",link_unit_test_framework="no")
			 AC_LANG_POP([C++])
               ])
                LDFLAGS="${saved_ldflags}"
			    if test "x$link_unit_test_framework" = "xyes"; then
                    BOOST_UNIT_TEST_FRAMEWORK_LIB="-l$ax_lib"
                    AC_SUBST(BOOST_UNIT_TEST_FRAMEWORK_LIB)
					break
				fi
              done
			    if test "x$link_unit_test_framework" = "xno"; then
				   AC_MSG_NOTICE(Could not link against $ax_lib !)
				fi
			fi
		fi
        CPPFLAGS="$CPPFLAGS_SAVED"
        LDFLAGS="$LDFLAGS_SAVED"
	fi
])
dnl @synopsis AC_PYTHON_DEVEL([version])
dnl
dnl Note: Defines as a precious variable "PYTHON_VERSION". Don't
dnl override it in your configure.ac.
dnl
dnl This macro checks for Python and tries to get the include path to
dnl 'Python.h'. It provides the $(PYTHON_CPPFLAGS) and
dnl $(PYTHON_LDFLAGS) output variables. It also exports
dnl $(PYTHON_EXTRA_LIBS) and $(PYTHON_EXTRA_LDFLAGS) for embedding
dnl Python in your code.
dnl
dnl You can search for some particular version of Python by passing a
dnl parameter to this macro, for example ">= '2.3.1'", or "== '2.4'".
dnl Please note that you *have* to pass also an operator along with the
dnl version to match, and pay special attention to the single quotes
dnl surrounding the version number. Don't use "PYTHON_VERSION" for
dnl this: that environment variable is declared as precious and thus
dnl reserved for the end-user.
dnl
dnl This macro should work for all versions of Python >= 2.1.0. As an
dnl end user, you can disable the check for the python version by
dnl setting the PYTHON_NOVERSIONCHECK environment variable to something
dnl else than the empty string.
dnl
dnl If you need to use this macro for an older Python version, please
dnl contact the authors. We're always open for feedback.
dnl
dnl @category InstalledPackages
dnl @author Sebastian Huber <sebastian-huber@web.de>
dnl @author Alan W. Irwin <irwin@beluga.phys.uvic.ca>
dnl @author Rafael Laboissiere <laboissiere@psy.mpg.de>
dnl @author Andrew Collier <colliera@nu.ac.za>
dnl @author Matteo Settenvini <matteo@member.fsf.org>
dnl @author Horst Knorr <hk_classes@knoda.org>
dnl @version 2006-05-27
dnl @license GPLWithACException

AC_DEFUN([AC_PYTHON_DEVEL],[
	#
	# Allow the use of a (user set) custom python version
	#
	AC_ARG_VAR([PYTHON_VERSION],[The installed Python
		version to use, for example '2.3'. This string
		will be appended to the Python interpreter
		canonical name.])

	AC_PATH_PROG([PYTHON],[python[$PYTHON_VERSION]])
	if test -z "$PYTHON"; then
	   AC_MSG_ERROR([Cannot find python$PYTHON_VERSION in your system path])
	   PYTHON_VERSION=""
	fi

	#
	# Check for a version of Python >= 2.1.0
	#
	AC_MSG_CHECKING([for a version of Python >= '2.1.0'])
	ac_supports_python_ver=`$PYTHON -c "import sys, string; \
		ver = string.split(sys.version)[[0]]; \
		print ver >= '2.1.0'"`
	if test "$ac_supports_python_ver" != "True"; then
		if test -z "$PYTHON_NOVERSIONCHECK"; then
			AC_MSG_RESULT([no])
			AC_MSG_FAILURE([
This version of the AC@&t@_PYTHON_DEVEL macro
doesn't work properly with versions of Python before
2.1.0. You may need to re-run configure, setting the
variables PYTHON_CPPFLAGS, PYTHON_LDFLAGS, PYTHON_SITE_PKG,
PYTHON_EXTRA_LIBS and PYTHON_EXTRA_LDFLAGS by hand.
Moreover, to disable this check, set PYTHON_NOVERSIONCHECK
to something else than an empty string.
])
		else
			AC_MSG_RESULT([skip at user request])
		fi
	else
		AC_MSG_RESULT([yes])
	fi

	#
	# if the macro parameter ``version'' is set, honour it
	#
	if test -n "$1"; then
		AC_MSG_CHECKING([for a version of Python $1])
		ac_supports_python_ver=`$PYTHON -c "import sys, string; \
			ver = string.split(sys.version)[[0]]; \
			print ver $1"`
		if test "$ac_supports_python_ver" = "True"; then
	   	   AC_MSG_RESULT([yes])
		else
			AC_MSG_RESULT([no])
			AC_MSG_ERROR([this package requires Python $1.
If you have it installed, but it isn't the default Python
interpreter in your system path, please pass the PYTHON_VERSION
variable to configure. See ``configure --help'' for reference.
])
			PYTHON_VERSION=""
		fi
	fi

	#
	# Check if you have distutils, else fail
	#
	AC_MSG_CHECKING([for the distutils Python package])
	ac_distutils_result=`$PYTHON -c "import distutils" 2>&1`
	if test -z "$ac_distutils_result"; then
		AC_MSG_RESULT([yes])
	else
		AC_MSG_RESULT([no])
		AC_MSG_ERROR([cannot import Python module "distutils".
Please check your Python installation. The error was:
$ac_distutils_result])
		PYTHON_VERSION=""
	fi

	#
	# Check for Python include path
	#
	AC_MSG_CHECKING([for Python include path])
	if test -z "$PYTHON_CPPFLAGS"; then
		python_path=`$PYTHON -c "import distutils.sysconfig; \
           		print distutils.sysconfig.get_python_inc();"`
		if test -n "${python_path}"; then
		   	python_path="-I$python_path"
		fi
		PYTHON_CPPFLAGS=$python_path
	fi
	AC_MSG_RESULT([$PYTHON_CPPFLAGS])
	AC_SUBST([PYTHON_CPPFLAGS])

	#
	# Check for Python library path
	#
	AC_MSG_CHECKING([for Python library path])
	if test -z "$PYTHON_LDFLAGS"; then
		# (makes two attempts to ensure we've got a version number
		# from the interpreter)
		py_version=`$PYTHON -c "from distutils.sysconfig import *; \
			from string import join; \
			print join(get_config_vars('VERSION'))"`
		if test "$py_version" == "[None]"; then
			if test -n "$PYTHON_VERSION"; then
				py_version=$PYTHON_VERSION
			else
				py_version=`$PYTHON -c "import sys; \
					print sys.version[[:3]]"`
			fi
		fi

		PYTHON_LDFLAGS=`$PYTHON -c "from distutils.sysconfig import *; \
			from string import join; \
			print '-L' + get_python_lib(0,1), \
		      	'-lpython';"`$py_version
	fi
	AC_MSG_RESULT([$PYTHON_LDFLAGS])
	AC_SUBST([PYTHON_LDFLAGS])

	#
	# Check for site packages
	#
	AC_MSG_CHECKING([for Python site-packages path])
	if test -z "$PYTHON_SITE_PKG"; then
		PYTHON_SITE_PKG=`$PYTHON -c "import distutils.sysconfig; \
		        print distutils.sysconfig.get_python_lib(0,0);"`
	fi
	AC_MSG_RESULT([$PYTHON_SITE_PKG])
	AC_SUBST([PYTHON_SITE_PKG])

	#
	# libraries which must be linked in when embedding
	#
	AC_MSG_CHECKING(python extra libraries)
	if test -z "$PYTHON_EXTRA_LIBS"; then
	   PYTHON_EXTRA_LIBS=`$PYTHON -c "import distutils.sysconfig; \
                conf = distutils.sysconfig.get_config_var; \
                print conf('LOCALMODLIBS'), conf('LIBS')"`
	fi
	AC_MSG_RESULT([$PYTHON_EXTRA_LIBS])
	AC_SUBST(PYTHON_EXTRA_LIBS)

	#
	# linking flags needed when embedding
	#
	AC_MSG_CHECKING(python extra linking flags)
	if test -z "$PYTHON_EXTRA_LDFLAGS"; then
		PYTHON_EXTRA_LDFLAGS=`$PYTHON -c "import distutils.sysconfig; \
			conf = distutils.sysconfig.get_config_var; \
			print conf('LINKFORSHARED')"`
	fi
	AC_MSG_RESULT([$PYTHON_EXTRA_LDFLAGS])
	AC_SUBST(PYTHON_EXTRA_LDFLAGS)

	#
	# final check to see if everything compiles alright
	#
	AC_MSG_CHECKING([consistency of all components of python development environment])
	AC_LANG_PUSH([C])
	# save current global flags
	LIBS="$ac_save_LIBS $PYTHON_LDFLAGS"
	CPPFLAGS="$ac_save_CPPFLAGS $PYTHON_CPPFLAGS"
	AC_TRY_LINK([
		#include <Python.h>
	],[
		Py_Initialize();
	],[pythonexists=yes],[pythonexists=no])

	AC_MSG_RESULT([$pythonexists])

        if test ! "$pythonexists" = "yes"; then
	   AC_MSG_ERROR([
  Could not link test program to Python. Maybe the main Python library has been
  installed in some non-standard library path. If so, pass it to configure,
  via the LDFLAGS environment variable.
  Example: ./configure LDFLAGS="-L/usr/non-standard-path/python/lib"
  ============================================================================
   ERROR!
   You probably have to install the development version of the Python package
   for your distribution.  The exact name of this package varies among them.
  ============================================================================
	   ])
	  PYTHON_VERSION=""
	fi
	AC_LANG_POP
	# turn back to default flags
	CPPFLAGS="$ac_save_CPPFLAGS"
	LIBS="$ac_save_LIBS"

	#
	# all done!
	#
])
dnl @synopsis AX_BOOST_PYTHON
dnl
dnl This macro checks to see if the Boost.Python library is installed.
dnl It also attempts to guess the currect library name using several
dnl attempts. It tries to build the library name using a user supplied
dnl name or suffix and then just the raw library.
dnl
dnl If the library is found, HAVE_BOOST_PYTHON is defined and
dnl BOOST_PYTHON_LIB is set to the name of the library.
dnl
dnl This macro calls AC_SUBST(BOOST_PYTHON_LIB).
dnl
dnl In order to ensure that the Python headers are specified on the
dnl include path, this macro requires AX_PYTHON to be called.
dnl
dnl @category InstalledPackages
dnl @category Cxx
dnl @author Michael Tindal <mtindal@paradoxpoint.com>
dnl @version 2005-05-20
dnl @license GPLWithACException

AC_DEFUN([AX_BOOST_PYTHON],
[AC_REQUIRE([AX_PYTHON])dnl
AC_CACHE_CHECK(whether the Boost::Python library is available,
ac_cv_boost_python,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 CPPFLAGS_SAVE=$CPPFLAGS
 if test x$PYTHON_INCLUDE_DIR != x; then
   CPPFLAGS=-I$PYTHON_INCLUDE_DIR $CPPFLAGS
 fi
 AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[
 #include <boost/python/module.hpp>
 using namespace boost::python;
 BOOST_PYTHON_MODULE(test) { throw "Boost::Python test."; }]],
 			   [[return 0;]]),
  			   ac_cv_boost_python=yes, ac_cv_boost_python=no)
 AC_LANG_RESTORE
 CPPFLAGS=$CPPFLAGS_SAVE
])
if test "$ac_cv_boost_python" = "yes"; then
  AC_DEFINE(HAVE_BOOST_PYTHON,,[define if the Boost::Python library is available])
  ax_python_lib=boost_python
  AC_ARG_WITH([boost-python],AS_HELP_STRING([--with-boost-python],[specify the boost python library or suffix to use]),
  [if test "x$with_boost_python" != "xno"; then
     ax_python_lib=$with_boost_python
     ax_boost_python_lib=boost_python-$with_boost_python
   fi])
  for ax_lib in $ax_python_lib $ax_boost_python_lib boost_python; do
    AC_CHECK_LIB($ax_lib, main, [BOOST_PYTHON_LIB=$ax_lib break])
  done
  AC_SUBST(BOOST_PYTHON_LIB)
fi
])dnl
