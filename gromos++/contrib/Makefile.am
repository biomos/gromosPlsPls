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
	cg_potential_restr\
        cg_energy_cg\
        cg_energy_gromacs\
        cg_coulomb\
        cg_startcrd\
        cg_startcrd2\
	cg_internal\
        cg_internal_restr\
        cg_rdf\
        cg_rdf_disres\
        cg_tser_pos\
        cg_dipole\
        cg_dipole_cg\
        cg_dipole_disres\
        cg_dipole_cg_series\
        cg_epsilon\
        cg_eps_restr\
        cg_eps_cg\
        cg_dist_restr\
        cg_cluster\
        cg_exclusive\
        cg_clust2disres\
        cg_rmsd_traj\
        cg_rmsd \
        cg_diffus\
        cg_int_dip\
	countbb\
	stacking\
	shake\
	model\
	gca_ener\
	renumber\
	distmat\
	iwdcf\
	bb2tex\
	crossdih_top\
	saltbridge\
	bilayer_gel\
	make_rdc_spec\
	svd_fit\
        sasa_new\
	add_atom\
	top_to_build\
	helix_axis\
	rdc_sh\
	gromacs2gromos\
	diff_top\
	matrix_overlap\
	grid_dens\
	pdb2seq\
	nscatt\
	structal\
	dimer\
	decomp_dip\
	cgLJpot\
	hvap\
	depsi2native\
	rep_rewrite

aggregates_SOURCES = aggregates.cc              
rottemp_SOURCES = rottemp.cc
ener_rplc_SOURCES = ener_rplc.cc
espmap_SOURCES = espmap.cc
lb_top_SOURCES = lb_top.cc
m_ener_SOURCES = m_ener.cc
prep_bb_SOURCES = prep_bb.cc
cg_potential_SOURCES = cg_potential.cc
cg_potential_restr_SOURCES = cg_potential_restr.cc
cg_energy_cg_SOURCES = cg_energy_cg.cc
cg_energy_gromacs_SOURCES = cg_energy_gromacs.cc
cg_coulomb_SOURCES = cg_coulomb.cc
cg_startcrd_SOURCES = cg_startcrd.cc
cg_startcrd2_SOURCES = cg_startcrd2.cc
cg_internal_SOURCES = cg_internal.cc
cg_internal_restr_SOURCES = cg_internal_restr.cc
cg_rdf_SOURCES = cg_rdf.cc
cg_rdf_disres_SOURCES = cg_rdf_disres.cc
cg_tser_pos_SOURCES = cg_tser_pos.cc
cg_dipole_SOURCES = cg_dipole.cc
cg_dipole_cg_SOURCES = cg_dipole_cg.cc
cg_dipole_disres_SOURCES = cg_dipole_disres.cc
cg_dipole_cg_series_SOURCES = cg_dipole_cg_series.cc
cg_epsilon_SOURCES = cg_epsilon.cc
cg_eps_restr_SOURCES = cg_eps_restr.cc
cg_eps_cg_SOURCES = cg_eps_cg.cc
cg_dist_restr_SOURCES = cg_dist_restr.cc
cg_cluster_SOURCES = cg_cluster.cc
cg_exclusive_SOURCES = cg_exclusive.cc
cg_clust2disres_SOURCES = cg_clust2disres.cc
cg_rmsd_traj_SOURCES = cg_rmsd_traj.cc
cg_rmsd_SOURCES = cg_rmsd.cc
cg_diffus_SOURCES = cg_diffus.cc
cg_int_dip_SOURCES = cg_int_dip.cc
countbb_SOURCES = countbb.cc
stacking_SOURCES = stacking.cc
shake_SOURCES = shake.cc
model_SOURCES = model.cc
gca_ener_SOURCES = gca_ener.cc
renumber_SOURCES = renumber.cc
distmat_SOURCES = distmat.cc
iwdcf_SOURCES = iwdcf.cc
bb2tex_SOURCES = bb2tex.cc
crossdih_top_SOURCES = crossdih_top.cc
saltbridge_SOURCES = saltbridge.cc
bilayer_gel_SOURCES = bilayer_gel.cc
make_rdc_spec_SOURCES = make_rdc_spec.cc
svd_fit_SOURCES = svd_fit.cc
sasa_new_SOURCES = sasa_new.cc
add_atom_SOURCES = add_atom.cc
top_to_build_SOURCES = top_to_build.cc
helix_axis_SOURCES = helix_axis.cc
rdc_sh_SOURCES = rdc_sh.cc
gromacs2gromos_SOURCES = gromacs2gromos.cc
diff_top_SOURCES = diff_top.cc
matrix_overlap_SOURCES = matrix_overlap.cc
grid_dens_SOURCES = grid_dens.cc
pdb2seq_SOURCES = pdb2seq.cc
nscatt_SOURCES = nscatt.cc
structal_SOURCES = structal.cc
dimer_SOURCES = dimer.cc
decomp_dip_SOURCES = decomp_dip.cc
cgLJpot_SOURCES = cgLJpot.cc
hvap_SOURCES = hvap.cc
depsi2native_SOURCES = depsi2native.cc
rep_rewrite_SOURCES = rep_rewrite.cc

LDADD = $(top_builddir)/src/libgromos.la

shake_LDADD = $(MDPP_LIB) $(top_builddir)/src/libgromos.la
shake_LDFLAGS = $(MDPP_LDFLAGS)

gromacs2gromos_LDADD = $(GMX_LIB) $(top_builddir)/src/libgromos.la
gromacs2gromos_LDFLAGS = $(GMX_LDFLAGS) $(top_builddir)/src/libgromos.la
