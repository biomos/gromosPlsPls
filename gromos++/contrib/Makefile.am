# This is the programs automake template.
AM_CXXFLAGS = @MY_CXXFLAGS@
INCLUDES = -I$(top_srcdir)/src

EXTRA_DIST =   0index_contrib.doxy

bin_PROGRAMS = aggregates\
               rottemp\
               ener_rplc\
               espmap\
               lb_top\
               m_ener\
               pairener\
               prep_bb\
               cg_potential\
               countbb\
               stacking\
               shake\
               model\
	       gca_ener\
	       renumber\
	       distmat\
	       iwdcf\
	       jepot\
	       solv_int_ener\
               alignment_tensor\
	       tesselate

aggregates_SOURCES = aggregates.cc              
rottemp_SOURCES = rottemp.cc
ener_rplc_SOURCES = ener_rplc.cc
espmap_SOURCES = espmap.cc
lb_top_SOURCES = lb_top.cc
m_ener_SOURCES = m_ener.cc
pairener_SOURCES = pairener.cc
prep_bb_SOURCES = prep_bb.cc
cg_potential_SOURCES = cg_potential.cc
countbb_SOURCES = countbb.cc
stacking_SOURCES = stacking.cc
shake_SOURCES = shake.cc
model_SOURCES = model.cc
gca_ener_SOURCES = gca_ener.cc
renumber_SOURCES = renumber.cc
distmat_SOURCES = distmat.cc
iwdcf_SOURCES = iwdcf.cc
jepot_SOURCES = jepot.cc
solv_int_ener_SOURCES = solv_int_ener.cc
tesselate_SOURCES = tesselate.cc
alignment_tensor_SOURCES = tesselate.cc

LDADD = $(top_builddir)/src/libgromos.la

shake_LDADD = $(GROMOSXX_LIB) $(top_builddir)/src/libgromos.la
shake_LDFLAGS = $(GROMOSXX_LDFLAGS)

