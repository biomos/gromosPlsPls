dnl Local function, as GROMOS++ depends
dnl on the STL. -> Does the C++ compiler
dnl support the STL to the degree necessary?
dnl 
AC_DEFUN([AC_CV_CXX_VERSION_OK],
  [AC_CACHE_CHECK(whether the compiler supports the STL,
   ac_cv_cxx_version_ok,
     [AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
      AC_TRY_COMPILE([#include <vector>],[
	std::vector<int> v; v.push_back(0);return 0;],
         ac_cv_cxx_version_ok=yes, ac_cv_cxx_version_ok=no)
      AC_LANG_RESTORE
  ])
])

dnl @synopsis AC_CXX_LIB_GROMOSXX([optional-string "required"])
dnl
dnl Check whether MD++ is installed.
dnl GromosXX is available from the
dnl group of computational chemistry of
dnl Prof. W. F. van Gunsteren
dnl Swiss Federal Institute of Technology Zurich
dnl
dnl   Set the path for MD++ with the option
dnl      --with-mdpp[=DIR]
dnl   Gromos md++ headers should be under DIR/include
dnl   Gromos md++ library should be under DIR/lib
dnl   Then try to compile and run a simple program with a Blitz Array
dnl   Optional argument `required' triggers an error if Blitz++ not installed
dnl 
dnl @version $Id$
dnl @author Patrick Guio <patrick.guio@matnat.uio.no>
dnl
AC_DEFUN([AC_MSG_ERROR_GROMOSXX],[
AC_MSG_ERROR([
$PACKAGE_STRING requires the Gromos md++ library
When installed give the directory of installation with the option
  --with-mdpp@<:@=DIR@:>@
])])


AC_DEFUN([AC_CXX_LIB_GROMOSXX],[

AC_ARG_WITH(mdpp,
AS_HELP_STRING([--with-mdpp@<:@=DIR@:>@],[Set the path for GROMOS MD++]),
[],[withval='yes'])

if test "$1" = required -a "$withval" = no ; then
	AC_MSG_ERROR_GROMOSXX
fi

if test "$withval" != no ; then

	saveCPPFLAGS=$CPPFLAGS
	saveLDFLAGS=$LDFLAGS
	saveLIBS=$LIBS

	if test "$withval" != 'yes'; then
		CPPFLAGS="-I$withval/include"
		LDFLAGS="-L$withval/lib -Wl,-R$withval/lib"
	fi
	LIBS="-lgroxx"

	AC_CACHE_CHECK([whether Gromos MD++ is installed],ac_cv_cxx_lib_gromosxx,
	[AC_LANG_SAVE
	AC_LANG_CPLUSPLUS
	AC_RUN_IFELSE(
	[AC_LANG_PROGRAM([[
#include <cmath>
#include <cassert>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <gromosXX/math/gmath.h>
]],[[
math::VArray p(100);
math::boundary_enum b = math::rectangular;
	]])],[ac_cxx_lib_gromosxx=yes],[ac_cxx_lib_gromosxx=no])
	AC_LANG_RESTORE
	])

	CPPFLAGS=$saveCPPFLAGS
	LDFLAGS=$saveLDFLAGS
	LIBS=$saveLIBS

	if test "$ac_cxx_lib_gromosxx" = yes ; then
		if test "$withval" != yes ; then
			CPPFLAGS="$CPPFLAGS -I$withval/include"
			GROMOSXX_LDFLAGS="-L$withval/lib -Wl,-R$withval/lib"
			AC_SUBST(GROMOSXX_LDFLAGS)
		fi
		GROMOSXX_LIB="-lgroxx"
		AC_SUBST(GROMOSXX_LIB)
   		AC_DEFINE_UNQUOTED([HAVE_GROMOSXX],[],[GromosXX library])

	else
		if test "$1" = required ; then
			AC_MSG_ERROR_GROMOSXX
		fi
	fi

fi

])

dnl @synopsis AC_CXX_LIB_GSL([optional-string "required"])
dnl
dnl Check whether Gnu Scientific Library (GSL) is installed.
dnl GSL is available from
dnl www.gnu.org
dnl
dnl   Set the path for GSL  with the option
dnl      --with-gsl[=DIR]
dnl   GSL headers should be under DIR/include
dnl   GSL library should be under DIR/lib
dnl   Then try to compile and run a simple program with a gsl random number
dnl   Optional argument `required' triggers an error if GSL not installed
dnl 
dnl @version $Id$
dnl @author Patrick Guio <patrick.guio@matnat.uio.no>
dnl
AC_DEFUN([AC_MSG_ERROR_GSL],[
AC_MSG_ERROR([
$PACKAGE_STRING requires the Gnu Scientific Library (GSL)
When installed give the directory of installation with the option
  --with-gsl@<:@=DIR@:>@
])])


AC_DEFUN([AC_CXX_LIB_GSL],[

AC_ARG_WITH(gsl,
AS_HELP_STRING([--with-gsl@<:@=DIR@:>@],[Set the path for GSL]),
[],[withval='yes'])

if test "$1" = required -a "$withval" = no ; then
	AC_MSG_ERROR_GSL
fi

if test "$withval" != no ; then

	saveCPPFLAGS=$CPPFLAGS
	saveLDFLAGS=$LDFLAGS
	saveLIBS=$LIBS

	if test "$withval" != 'yes'; then
		CPPFLAGS="-I$withval/include"
		LDFLAGS="-L$withval/lib -Wl,-R$withval/lib"
	fi
	LIBS="-lgsl -lgslcblas"

	AC_CACHE_CHECK([whether Gnu Scientific Library is installed],ac_cv_cxx_lib_gsl,
	[AC_LANG_SAVE
	AC_LANG_CPLUSPLUS
	AC_RUN_IFELSE(
	[AC_LANG_PROGRAM([[
#include <gsl/gsl_matrix.h>
]],[[
gsl_matrix * mat = gsl_matrix_alloc(3,3);
gsl_matrix_set_zero(mat);
gsl_matrix_free(mat);
	]])],[ac_cxx_lib_gsl=yes],[ac_cxx_lib_gsl=no])
	AC_LANG_RESTORE
	])

	CPPFLAGS=$saveCPPFLAGS
	LDFLAGS=$saveLDFLAGS
	LIBS=$saveLIBS

	if test "$ac_cxx_lib_gsl" = yes ; then
		if test "$withval" != yes ; then
			CPPFLAGS="$CPPFLAGS -I$withval/include"
			GSL_LDFLAGS="-L$withval/lib -Wl,-R$withval/lib"
			AC_SUBST(GSL_LDFLAGS)
		fi
		GSL_LIB="-lgsl -lgslcblas"
		AC_SUBST(GSL_LIB)

   		AC_DEFINE_UNQUOTED([HAVE_GSL],[],[Gnu Scientific Library])

	else
		if test "$1" = required ; then
			AC_MSG_ERROR_GSL
		fi
	fi

fi

])
AC_DEFUN([AC_PROG_CXX_SUNCC],
[AC_CACHE_CHECK(whether we are using Sun C++, SUN_cv_CXX,
[cat > conftest.c <<EOF
# if defined(__SUNPRO_CC)
  yes;
#endif
EOF
if AC_TRY_COMMAND(${CXX} -E conftest.c) | egrep yes >/dev/null 2>&1; then
  SUN_CXX=yes
  compiler=suncc
else
  SUN_CXX=no
fi])])

AC_DEFUN([AC_PROG_CXX_INTELCC],
[AC_CACHE_CHECK(whether we are using Intel C++, INTEL_cv_CXX,
[cat > conftest.c <<EOF
# if defined(__ICC)
  yes;
#endif
EOF
if AC_TRY_COMMAND(${CXX} -E conftest.c) | egrep yes >/dev/null 2>&1; then
  INTEL_CXX=yes
  compiler=intelcc
else
  INTEL_CXX=no
fi])])

dnl check for lib CCP4/Clipper
AC_DEFUN([AM_PATH_CCP4_CLIPPER],[
  dnl allow for ccp4 lib directory specification
  AC_ARG_WITH(ccp4,
    [  --with-ccp4=DIR     CCP4 library directory to use],
    [
      [CXXFLAGS="$CXXFLAGS -I${withval}/include -L${withval}/lib"]
      [LDFLAGS="$LDFLAGS -L${withval}/lib"]
    ],
    [
      AC_MSG_WARN([Assuming default paths for CCP4])
    ])
  AC_ARG_WITH(clipper,
    [  --with-clipper=DIR  clipper library directory to use],
    [
      [CXXFLAGS="$CXXFLAGS -I${withval}/include -L${withval}/lib"]
      [LDFLAGS="$LDFLAGS -L${withval}/lib"]
      [CLIPPER_LIB="-lclipper-ccp4 -lccp4c -lclipper-contrib -lclipper-core -lrfftw -lfftw -lpthread"]
      AC_DEFINE_UNQUOTED([HAVE_CLIPPER],[],[Have clipper x-ray library])
    ],
    [
      AC_MSG_WARN([clipper path was not specified. Disabling clipper support])
      [CLIPPER_LIB=""]
    ]
  )
  dnl check for lib with these settings. To be implemented
  AC_SUBST(CLIPPER_LIB)
])


