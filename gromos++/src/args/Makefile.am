## Template for automake for libargs

lib_LTLIBRARIES = libargs.la
libargs_la_SOURCES = Arguments.cc\
                     BoundaryParser.cc\
                     GatherParser.cc

libargs_la_LDFLAGS = -version-info $(GROMOS_MAJOR_VERSION):$(GROMOS_MICRO_VERSION):0

noinst_PROGRAMS = Arguments

LDADD = libargs.la\
        ../bound/libbound.la\
        ../gio/libgio.la\
        ../gcore/libgcore.la\
        ../gmath/libgmath.la

Arguments_SOURCES = Arguments.t.cc
