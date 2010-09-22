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
	duplicate\
	prep_xray\
	structure_factor\
	xrayts\
	r_factor\
	xray_map\
	bb2tex\
	crossdih_top\
	saltbridge\
	bilayer_gel\
	make_rdc_spec\
	svd_fit\
	add_atom\
	r_real_factor\
	prep_xray_le\
	r_real_factor\
	top_to_build\
	helix_axis\
	rdc_sh\
	gromacs2gromos\
	diff_top\
	matrix_overlap\
	close_pair

aggregates_SOURCES = aggregates.cc              
rottemp_SOURCES = rottemp.cc
ener_rplc_SOURCES = ener_rplc.cc
espmap_SOURCES = espmap.cc
lb_top_SOURCES = lb_top.cc
m_ener_SOURCES = m_ener.cc
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
duplicate_SOURCES = duplicate.cc
prep_xray_SOURCES = prep_xray.cc
structure_factor_SOURCES = structure_factor.cc
xrayts_SOURCES = xrayts.cc
r_factor_SOURCES = r_factor.cc
xray_map_SOURCES = xray_map.cc
bb2tex_SOURCES = bb2tex.cc
crossdih_top_SOURCES = crossdih_top.cc
saltbridge_SOURCES = saltbridge.cc
bilayer_gel_SOURCES = bilayer_gel.cc
make_rdc_spec_SOURCES = make_rdc_spec.cc
svd_fit_SOURCES = svd_fit.cc
add_atom_SOURCES = add_atom.cc
r_real_factor_SOURCES = r_real_factor.cc
prep_xray_le_SOURCES = prep_xray_le.cc
top_to_build_SOURCES = top_to_build.cc
helix_axis_SOURCES = helix_axis.cc
rdc_sh_SOURCES = rdc_sh.cc
gromacs2gromos_SOURCES = gromacs2gromos.cc
diff_top_SOURCES = diff_top.cc
matrix_overlap_SOURCES = matrix_overlap.cc
close_pair_SOURCES = close_pair.cc

LDADD = $(top_builddir)/src/libgromos.la

shake_LDADD = $(MDPP_LIB) $(top_builddir)/src/libgromos.la
shake_LDFLAGS = $(MDPP_LDFLAGS)

structure_factor_LDADD = $(CLIPPER_LIB) $(top_builddir)/src/libgromos.la
structure_factor_LDFLAGS = $(CLIPPER_LDFLAGS)
prep_xray_LDADD = $(CLIPPER_LIB) $(top_builddir)/src/libgromos.la
prep_xray_LDFLAGS = $(CLIPPER_LDFLAGS)
r_factor_LDADD = $(CLIPPER_LIB) $(top_builddir)/src/libgromos.la
r_factor_LDFLAGS = $(CLIPPER_LDFLAGS)
xray_map_LDADD = $(CLIPPER_LIB) $(top_builddir)/src/libgromos.la
xray_map_LDFLAGS = $(CLIPPER_LDFLAGS)
r_real_factor_LDADD = $(CLIPPER_LIB) $(top_builddir)/src/libgromos.la
r_real_factor_LDFLAGS = $(CLIPPER_LDFLAGS)

gromacs2gromos_LDADD = $(GMX_LIB) $(top_builddir)/src/libgromos.la
gromacs2gromos_LDFLAGS = $(GMX_LDFLAGS) $(top_builddir)/src/libgromos.la
