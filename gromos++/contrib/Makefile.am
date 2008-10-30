# This is the programs automake template.
AM_CXXFLAGS = @MY_CXXFLAGS@
INCLUDES = -I$(top_srcdir)/src

EXTRA_DIST =   0index_contrib.doxy

bin_PROGRAMS = aggregates\
               rottemp\
               ener_rplc\
               pairener\
               prep_bb\
               cg_potential\
               countbb\
               stacking\
               shake\
               model

aggregates_SOURCES = aggregates.cc              
rottemp_SOURCES = rottemp.cc
ener_rplc_SOURCES = ener_rplc.cc
pairener_SOURCES = pairener.cc
prep_bb_SOURCES = prep_bb.cc
cg_potential_SOURCES = cg_potential.cc
countbb_SOURCES = countbb.cc
stacking_SOURCES = stacking.cc
shake_SOURCES = shake.cc
model_SOURCES = model.cc

LDADD = $(top_builddir)/src/libgromos.la

shake_LDADD = $(GROMOSXX_LIB) $(top_builddir)/src/libgromos.la
shake_LDFLAGS = $(GROMOSXX_LDFLAGS)

