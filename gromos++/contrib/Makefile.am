# This is the programs automake template.

INCLUDES = -I$(top_srcdir)/src

EXTRA_DIST =   0index_contrib.doxy

bin_PROGRAMS = aggregates\
               rottemp

aggregates_SOURCES = aggregates.cc              
rottemp_SOURCES = rottemp.cc

LDADD = $(top_builddir)/src/libgromos.la
