
noinst_LTLIBRARIES = libfit.la

gincludedir = $(includedir)/gromos++

ginclude_HEADERS = PositionUtils.h \
		  Reference.h \
		  RotationalFit.h \
		  TranslationalFit.h\
                  FastRotationalFit.h

libfit_la_SOURCES = PositionUtils.cc\
                    Reference.cc\
                    RotationalFit.cc\
                    TranslationalFit.cc\
                    FastRotationalFit.cc

check_PROGRAMS = PositionUtils\
               Reference\
               RotationalFit\
               TranslationalFit

LDADD = ../libgromos.la

PositionUtils_SOURCES = PositionUtils.t.cc
Reference_SOURCES = Reference.t.cc
RotationalFit_SOURCES = RotationalFit.t.cc
TranslationalFit_SOURCES = TranslationalFit.t.cc


PositionUtils_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

