## Template for automake for libfit

lib_LTLIBRARIES = libfit.la

include_HEADERS = PositionUtils.h \
				  Reference.h \
				  RotationalFit.h \
				  TranslationalFit.h

libfit_la_SOURCES = PositionUtils.cc\
                    Reference.cc\
                    RotationalFit.cc\
                    TranslationalFit.cc

libfit_la_LDFLAGS = -version-info $(GROMOS_MAJOR_VERSION):$(GROMOS_MICRO_VERSION):0

libfit_la_LIBADD = $(top_builddir)/src/gmath/libgmath.la\
				   $(top_builddir)/src/gio/libgio.la\
				   $(top_builddir)/src/gcore/libgcore.la\
				   $(top_builddir)/src/gmath/libgmath.la
				   

check_PROGRAMS = PositionUtils\
               Reference\
               RotationalFit\
               TranslationalFit

LDADD = libfit.la

PositionUtils_SOURCES = PositionUtils.t.cc
Reference_SOURCES = Reference.t.cc
RotationalFit_SOURCES = RotationalFit.t.cc
TranslationalFit_SOURCES = TranslationalFit.t.cc
