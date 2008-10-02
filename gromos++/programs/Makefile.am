# This is the programs automake template.
AM_CXXFLAGS = @MY_CXXFLAGS@
INCLUDES = -I$(top_srcdir)/src

include_HEADERS = mk_script.h
		  
EXTRA_DIST =   0index.doxy

bin_PROGRAMS = rmsd\
               pdb2g96\
               ener\
               pairener\
               dg_ener\
	       m_ener\
               espmap\
               noe\
               prep_noe\
               post_noe\
               rgyr\
               tser\
	       ditrans \
               diffus\
               tstrip\
               gathtraj\
               edyn\
               rmsf\
               dipole\
               rmsdmat\
	       distmat\
               cluster\
               postcluster\
               check_box\
               make_top\
               con_top\
               com_top\
               red_top\
               pt_top\
               m_pt_top\
               LB_top\
               copy_box\
               build_box\
               bin_box\
               ran_box\
               ran_solvation\
               pert_top\
               explode\
               hbond\
               frameout\
               oparam\
               nhoparam\
               rdf\
               mdf\
               iwdcf\
               clusmol\
               iondens\
               mk_script\
               check_top\
               ene_ana\
               visco_ana\
               epsilon\
               gch\
               sim_box\
               ion\
	       gca\
	       gca_ener\
	       dssp\
               atominfo\
               shake_analysis\
               procmt_list\
               filter\
	       follow\
               renumber\
               countbb\
               jval\
               prepbbb\
               sasa\
	       shake\
	       tcf\
               m_widom\
               cog\
               rot_rel\
               inbox\
	       cry\
	       unify_box\
               make_pt_top\
               stacking\
               solute_entropy


shake_LDADD = $(BLITZ_LIB) $(GROMOSXX_LIB) $(top_builddir)/src/libgromos.la
shake_LDFLAGS = $(BLITZ_LDFLAGS) $(GROMOSXX_LDFLAGS)

oparam_SOURCES = oparam.cc              
nhoparam_SOURCES = nhoparam.cc
frameout_SOURCES = frameout.cc
rmsd_SOURCES = rmsd.cc
pdb2g96_SOURCES = pdb2g96.cc
ener_SOURCES = ener.cc
pairener_SOURCES = pairener.cc
dg_ener_SOURCES = dg_ener.cc
m_ener_SOURCES = m_ener.cc
espmap_SOURCES = espmap.cc
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
distmat_SOURCES = distmat.cc
postcluster_SOURCES = postcluster.cc
check_box_SOURCES = check_box.cc
make_top_SOURCES = make_top.cc
con_top_SOURCES = con_top.cc
com_top_SOURCES = com_top.cc
red_top_SOURCES = red_top.cc
pt_top_SOURCES = pt_top.cc
m_pt_top_SOURCES = m_pt_top.cc
LB_top_SOURCES = LB_top.cc
copy_box_SOURCES = copy_box.cc
build_box_SOURCES = build_box.cc
bin_box_SOURCES = bin_box.cc
ran_box_SOURCES = ran_box.cc
ran_solvation_SOURCES = ran_solvation.cc
pert_top_SOURCES = pert_top.cc
explode_SOURCES = explode.cc
rdf_SOURCES = rdf.cc
mdf_SOURCES = mdf.cc
iwdcf_SOURCES = iwdcf.cc
clusmol_SOURCES = clusmol.cc
iondens_SOURCES = iondens.cc
mk_script_SOURCES = mk_script.h mk_script.cc
check_top_SOURCES = check_top.cc
ene_ana_SOURCES = ene_ana.cc
visco_ana_SOURCES = visco_ana.cc
epsilon_SOURCES = epsilon.cc
gch_SOURCES = gch.cc
sim_box_SOURCES = sim_box.cc
ion_SOURCES = ion.cc
gca_SOURCES = gca.cc
gca_ener_SOURCES = gca_ener.cc
dssp_SOURCES = dssp.cc
atominfo_SOURCES = atominfo.cc
shake_analysis_SOURCES = shake_analysis.cc
procmt_list_SOURCES = procmt_list.cc
filter_SOURCES = filter.cc
follow_SOURCES = follow.cc
renumber_SOURCES = renumber.cc
countbb_SOURCES = countbb.cc
jval_SOURCES = jval.cc
prepbbb_SOURCES = prepbbb.cc
sasa_SOURCES = sasa.cc
shake_SOURCES = shake.cc
tcf_SOURCES = tcf.cc
m_widom_SOURCES = m_widom.cc
cog_SOURCES = cog.cc
rot_rel_SOURCES = rot_rel.cc
inbox_SOURCES = inbox.cc
cry_SOURCES = cry.cc
unify_box_SOURCES = unify_box.cc
make_pt_top_SOURCES = make_pt_top.cc
stacking_SOURCES = stacking.cc
solute_entropy_SOURCES = solute_entropy.cc

LDADD = $(top_builddir)/src/libgromos.la
