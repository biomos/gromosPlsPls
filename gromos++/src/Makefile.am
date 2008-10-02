## Template for automake for libgromos
AM_CXXFLAGS = @MY_CXXFLAGS@
SUBDIRS = gromos gmath gcore gio bound utils fit args .

EXTRA_DIST = Doxygen.h

lib_LTLIBRARIES = libgromos.la

libgromos_la_LDFLAGS = -no-undefined -version-info $(GROMOS_MAJOR_VERSION):$(GROMOS_MICRO_VERSION):0 $(GSL_LDFLAGS)

libgromos_la_SOURCES =

libgromos_la_LIBADD =	gmath/libgmath.la \
			gcore/libgcore.la \
			gio/libgio.la \
			bound/libbound.la \
			utils/libutils.la \
			fit/libfit.la \
			args/libargs.la \
			$(GSL_LIB)




