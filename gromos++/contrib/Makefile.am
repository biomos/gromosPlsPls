# This is the programs automake template.
AM_CXXFLAGS = @MY_CXXFLAGS@
INCLUDES = -I$(top_srcdir)/src

EXTRA_DIST =   0index_contrib.doxy

bin_PROGRAMS = aggregates\
               rottemp\
               ener_rplc

aggregates_SOURCES = aggregates.cc              
rottemp_SOURCES = rottemp.cc
ener_rplc_SOURCES = ener_rplc.cc

LDADD = $(top_builddir)/src/libgromos.la
