## Template for automake for libgromos

SUBDIRS = utils test err block const complex vector matrix permutation cblas blas linalg eigen sys

SUBLIBS = block/libgslblock.la blas/libgslblas.la complex/libgslcomplex.la\
          eigen/libgsleigen.la err/libgslerr.la vector/libgslvector.la\
          matrix/libgslmatrix.la utils/libutils.la test/libgsltest.la cblas/libgslcblas.la\
          linalg/libgsllinalg.la permutation/libgslpermutation.la sys/libgslsys.la

pkginclude_HEADERS = header/gsl_math.h header/gsl_mode.h header/gsl_math.h header/gsl_nan.h\
                     header/gsl_machine.h header/gsl_precision.h gsl_version.h

lib_LTLIBRARIES = libgsl.la
libgsl_la_SOURCES = version.c
libgsl_la_LIBADD = $(SUBLIBS)
noinst_HEADERS = header/templates_on.h header/templates_off.h

## will not be built during development
## lib_LTLIBRARIES = libgromos.la

## List to be completed...
## libgromos_la_SOURCES = ...

## libgromos_la_LDFLAGS = -version-info $(GROMOS_MAJOR_VERSION):$(GROMOS_MICRO_VERSION):0


