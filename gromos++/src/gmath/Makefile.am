## Template for automake for libgromos

# noinst_LIBRARIES = libgmath.a

lib_LTLIBRARIES = libgmath.la

libgmath_la_LDFLAGS = -version-info $(GROMOS_MAJOR_VERSION):$(GROMOS_MICRO_VERSION):0

include_HEADERS     = Vec.h \
                      Matrix.h \
                      Distribution.h \
                      Stat.h  \
                      Expression.h \
                      Correlation.h \
                      physics.h

libgmath_la_SOURCES = Vec.cc\
                     Matrix.cc\
                     Distribution.cc\
                     Stat.cc\
                     Expression.cc\
                     Correlation.cc

check_PROGRAMS = Vec\
               Matrix\
               Distribution\
               Expression\
               Correlation


Vec_SOURCES = Vec.t.cc
Matrix_SOURCES = Matrix.t.cc
Distribution_SOURCES = Distribution.t.cc
Expression_SOURCES = Expression.t.cc
Correlation_SOURCES = Correlation.t.cc

Vec_LDADD = libgmath.la
Distribution_LDADD = libgmath.la
Matrix_LDADD = libgmath.la
Expression_LDADD = libgmath.la
Correlation_LDADD = libgmath.la
