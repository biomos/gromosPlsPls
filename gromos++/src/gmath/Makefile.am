## Template for automake for libgromos

noinst_LTLIBRARIES = libgmath.la

gincludedir = $(includedir)/gromos++

ginclude_HEADERS     = Vec.h \
                      Matrix.h \
                      Distribution.h \
                      Stat.h \
                      Stat.cc \
                      Expression.h \
                      Correlation.h \
                      Physics.h

libgmath_la_SOURCES = Vec.cc\
                     Matrix.cc\
                     Distribution.cc\
                     Expression.cc\
                     Correlation.cc\
                     Physics.cc

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

LDADD = libgmath.la -lgslcblas -lgsl
