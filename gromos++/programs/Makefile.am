# This is the programs automake template.
AM_CXXFLAGS = @MY_CXXFLAGS@
INCLUDES = -I$(top_srcdir)/src

include_HEADERS = mk_script.h
		  
EXTRA_DIST =   0index.doxy

bin_PROGRAMS = rmsd\
               pdb2g96\
               ener\
               dg_ener\
               noe\
               prep_noe\
               post_noe\
               rgyr\
               tser\
	             ditrans\
               diffus\
               tstrip\
               gathtraj\
               edyn\
               rmsf\
               dipole\
               rmsdmat\
               cluster\
               postcluster\
               check_box\
               make_top\
               make_sasa_top\
               con_top\
               com_top\
               red_top\
               pt_top\
               copy_box\
               build_box\
               bin_box\
               ran_box\
               ran_solvation\
               pert_top\
               explode\
               hbond\
               frameout\
               nhoparam\
               rdf\
               mdf\
               iondens\
               mk_script\
               check_top\
               ene_ana\
               int_ener\
               visco\
               epsilon\
               gch\
               sim_box\
               ion\
	             gca\
	             dssp\
               atominfo\
               shake_analysis\
               filter\
	             follow\
               jval\
               jepot\
               sasa\
               sasa_hasel\
	             tcf\
               m_widom\
               cog\
               rot_rel\
               inbox\
	             cry\
	             unify_box\
               make_pt_top\
               solute_entropy\
               pairlist\
               prep_eds\
               dfmult\
               reweight\
       	       stat_ana\
               bilayer_dist

nhoparam_SOURCES = nhoparam.cc
frameout_SOURCES = frameout.cc
rmsd_SOURCES = rmsd.cc
pdb2g96_SOURCES = pdb2g96.cc
ener_SOURCES = ener.cc
dg_ener_SOURCES = dg_ener.cc
hbond_SOURCES = hbond.cc
noe_SOURCES = noe.cc
prep_noe_SOURCES = prep_noe.cc
post_noe_SOURCES = post_noe.cc
rgyr_SOURCES = rgyr.cc
tser_SOURCES = tser.cc
ditrans_SOURCES = ditrans.cc
diffus_SOURCES = diffus.cc
tstrip_SOURCES = tstrip.cc
gathtraj_SOURCES = gathtraj.cc
edyn_SOURCES = edyn.cc
rmsf_SOURCES = rmsf.cc
dipole_SOURCES = dipole.cc
cluster_SOURCES = cluster.cc
rmsdmat_SOURCES = rmsdmat.cc
postcluster_SOURCES = postcluster.cc
check_box_SOURCES = check_box.cc
make_top_SOURCES = make_top.cc
make_sasa_top_SOURCES = make_sasa_top.cc
con_top_SOURCES = con_top.cc
com_top_SOURCES = com_top.cc
red_top_SOURCES = red_top.cc
pt_top_SOURCES = pt_top.cc
copy_box_SOURCES = copy_box.cc
build_box_SOURCES = build_box.cc
bin_box_SOURCES = bin_box.cc
ran_box_SOURCES = ran_box.cc
ran_solvation_SOURCES = ran_solvation.cc
pert_top_SOURCES = pert_top.cc
explode_SOURCES = explode.cc
rdf_SOURCES = rdf.cc
docf_SOURCES = docf.cc
solv_shell_SOURCES = solv_shell.cc
mdf_SOURCES = mdf.cc
iondens_SOURCES = iondens.cc
mk_script_SOURCES = mk_script.h mk_script.cc
check_top_SOURCES = check_top.cc
ene_ana_SOURCES = ene_ana.cc
int_ener_SOURCES = int_ener.cc
visco_SOURCES = visco.cc
epsilon_SOURCES = epsilon.cc
gch_SOURCES = gch.cc
sim_box_SOURCES = sim_box.cc
ion_SOURCES = ion.cc
gca_SOURCES = gca.cc
dssp_SOURCES = dssp.cc
atominfo_SOURCES = atominfo.cc
shake_analysis_SOURCES = shake_analysis.cc
filter_SOURCES = filter.cc
follow_SOURCES = follow.cc
jval_SOURCES = jval.cc
jepot_SOURCES = jepot.cc
sasa_SOURCES = sasa.cc
sasa_hasel_SOURCES = sasa_hasel.cc
tcf_SOURCES = tcf.cc
m_widom_SOURCES = m_widom.cc
cog_SOURCES = cog.cc
rot_rel_SOURCES = rot_rel.cc
inbox_SOURCES = inbox.cc
cry_SOURCES = cry.cc
unify_box_SOURCES = unify_box.cc
make_pt_top_SOURCES = make_pt_top.cc
solute_entropy_SOURCES = solute_entropy.cc
pairlist_SOURCES = pairlist.cc
prep_eds_SOURCES = prep_eds.cc
dfmult_SOURCES = dfmult.cc
reweight_SOURCES = reweight.cc
stat_ana_SOURCES = stat_ana.cc
bilayer_dist_SOURCES = bilayer_dist.cc

LDADD = $(top_builddir)/src/libgromos.la
cry_LDADD = $(CLIPPER_LIB) $(top_builddir)/src/libgromos.la

