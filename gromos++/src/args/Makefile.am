AM_CXXFLAGS = @MY_CXXFLAGS@
AM_CPPFLAGS = -I$(top_srcdir)/src
noinst_LTLIBRARIES = libargs.la

libargs_la_SOURCES = Arguments.cc\
	BoundaryParser.cc\
	ReferenceParser.cc\
	GatherParser.cc\
	OutformatParser.cc

gincludedir = $(includedir)/gromos++

ginclude_HEADERS = Arguments.h\
	BoundaryParser.h\
	ReferenceParser.h\
	GatherParser.h\
	OutformatParser.h

check_PROGRAMS = Arguments

LDADD = libargs.la

Arguments_SOURCES = Arguments.t.cc
