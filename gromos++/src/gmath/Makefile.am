## Template for automake for libgromos

# noinst_LIBRARIES = libgmath.a

lib_LTLIBRARIES = libgmath.la

libgmath_la_LDFLAGS = -version-info $(GROMOS_MAJOR_VERSION):$(GROMOS_MICRO_VERSION):0

libgmath_la_SOURCES = Vec.cc\
                     Matrix.cc\
                     Distribution.cc\
                     Stat.cc

libgmath_la_LIBADD = $(SUBLIBS)

bin_PROGRAMS = Vec\
               Matrix\
               Distribution


Vec_SOURCES = Vec.t.cc
Matrix_SOURCES = Matrix.t.cc
Distribution_SOURCES = Distribution.t.cc

Vec_LDADD = libgmath.la
Distribution_LDADD = libgmath.la
Matrix_LDADD = libgmath.la
