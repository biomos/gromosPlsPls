## Template for automake for libgromos

# noinst_LIBRARIES = libbound.a

lib_LTLIBRARIES = libbound.la

include_HEADERS     = Boundary.h \
					  Vacuum.h \
					  TruncOct.h \
					  RectBox.h

libbound_la_LDFLAGS = -version-info $(GROMOS_MAJOR_VERSION):$(GROMOS_MICRO_VERSION):0

libbound_la_SOURCES = Boundary.cc\
                     Vacuum.cc\
                     TruncOct.cc\
                     RectBox.cc

check_PROGRAMS = TruncOct RectBox

LDADD = ../gio/libgio.la ../gcore/libgcore.la libbound.la

TruncOct_SOURCES = TruncOct.t.cc 
RectBox_SOURCES = RectBox.t.cc
