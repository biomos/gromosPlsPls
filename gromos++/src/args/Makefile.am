
noinst_LTLIBRARIES = libargs.la

libargs_la_SOURCES = Arguments.cc\
                     BoundaryParser.cc\
		     ReferenceParser.cc\
                     GatherParser.cc

gincludedir = $(includedir)/gromos++

ginclude_HEADERS    = Arguments.h \
		     BoundaryParser.h \
	             ReferenceParser.h \
		     GatherParser.h

check_PROGRAMS = Arguments

LDADD = libargs.la

Arguments_SOURCES = Arguments.t.cc
