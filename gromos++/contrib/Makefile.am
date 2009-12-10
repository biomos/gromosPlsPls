# This is the programs automake template.
AM_CXXFLAGS = @MY_CXXFLAGS@
INCLUDES = -I$(top_srcdir)/src

EXTRA_DIST =   0index_contrib.doxy

bin_PROGRAMS = \
	aggregates\
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
	put_into_box\
	gca_ener\
	renumber\
	distmat\
	iwdcf\
	duplicate\
	prep_xray\
	structure_factor\
	xrayts\
	r_factor\
	xray_map\
	bb2tex\
	crossdih_top\
	oparam\
	saltbridge\
	bilayer_gel\
        svd_fit

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
put_into_box_SOURCES = put_into_box.cc
duplicate_SOURCES = duplicate.cc
prep_xray_SOURCES = prep_xray.cc
structure_factor_SOURCES = structure_factor.cc
xrayts_SOURCES = xrayts.cc
r_factor_SOURCES = r_factor.cc
xray_map_SOURCES = xray_map.cc
bb2tex_SOURCES = bb2tex.cc
crossdih_top_SOURCES = crossdih_top.cc
oparam_SOURCES = oparam.cc
saltbridge_SOURCES = saltbridge.cc
bilayer_gel_SOURCES = bilayer_gel.cc
svd_fit_SOURCES = svd_fit.cc

LDADD = $(top_builddir)/src/libgromos.la

shake_LDADD = $(GROMOSXX_LIB) $(top_builddir)/src/libgromos.la
shake_LDFLAGS = $(GROMOSXX_LDFLAGS)

structure_factor_LDADD = $(CLIPPER_LIB) $(top_builddir)/src/libgromos.la
structure_factor_LDFLAGS = $(CLIPPER_LDFLAGS)
prep_xray_LDADD = $(CLIPPER_LIB) $(top_builddir)/src/libgromos.la
prep_xray_LDFLAGS = $(CLIPPER_LDFLAGS)
r_factor_LDADD = $(CLIPPER_LIB) $(top_builddir)/src/libgromos.la
r_factor_LDFLAGS = $(CLIPPER_LDFLAGS)
xray_map_LDADD = $(CLIPPER_LIB) $(top_builddir)/src/libgromos.la
xray_map_LDFLAGS = $(CLIPPER_LDFLAGS)
