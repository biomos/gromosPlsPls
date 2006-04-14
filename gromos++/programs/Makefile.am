# This is the programs automake template.

INCLUDES = -I$(top_srcdir)/src

include_HEADERS = mkscript.h
EXTRA_DIST =   0index.doxy

bin_PROGRAMS = rmsd\
               pdb2g96\
               ener\
               nb_ener\
               pairener\
               dg_ener\
               espmap\
               convert\
               noe\
               prepnoe\
               noenew\
               prepnoenew\
               postnoe\
               rgyr\
               vmdam\
               tser\
	       rottemp \
	       ditrans \
               cov_ener\
               diffus\
               tstrip\
               gathtraj\
               edyn\
               rmsf\
               dipole\
               rmsdmat\
               cluster\
               rmsdmat2\
               cluster2\
               postcluster\
               checkbox\
               maketop\
               contop\
               comtop\
               redtop\
               pt_top\
               m_pt_top\
               LB_top\
               copy_box\
               buildbox\
               bin_box\
               ran_box\
               ran_solvation\
               perttop\
               pert_top\
               explode\
               hbond\
               frameout\
               oparam\
               nhoparam\
               rdf\
               mdf\
               iwdcf\
               neighbours\
               clusmol\
               methanes\
               mcti\
               iondens\
               mkscript\
               checktopo\
               check_top\
               ene_ana\
               visco_ana\
               epsilon\
               progch\
               probox\
               proion\
               progca\
	       progca_series \
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
               proarse\
	       shake\
	       protcf\
               m_widom\
               cog\
               rot_rel\
               inbox\
		procry\
		unify_box\
		monoconv


shake_LDADD = $(BLITZ_LIB) $(GROMOSXX_LIB) $(top_builddir)/src/libgromos.la
shake_LDFLAGS = $(BLITZ_LDFLAGS) $(GROMOSXX_LDFLAGS)

rmsdmat2_CXXFLAGS = $(LFS)
rmsdmat2_LDFLAGS = $(LFS)

cluster2_CXXFLAGS = $(LFS)
cluster2_LDFLAGS = $(LFS)

oparam_SOURCES = oparam.cc              
nhoparam_SOURCES = nhoparam.cc
frameout_SOURCES = frameout.cc
rmsd_SOURCES = rmsd.cc
pdb2g96_SOURCES = pdb2g96.cc
ener_SOURCES = ener.cc
nb_ener_SOURCES = nb_ener.cc
pairener_SOURCES = pairener.cc
dg_ener_SOURCES = dg_ener.cc
espmap_SOURCES = espmap.cc
hbond_SOURCES = hbond.cc
convert_SOURCES = convert.cc
noe_SOURCES = noe.cc
prepnoe_SOURCES = prepnoe.cc
noenew_SOURCES = noenew.cc
prepnoenew_SOURCES = prepnoenew.cc
postnoe_SOURCES = postnoe.cc
rgyr_SOURCES = rgyr.cc
vmdam_SOURCES = vmdam.cc
tser_SOURCES = tser.cc
rottemp_SOURCES = rottemp.cc
ditrans_SOURCES = ditrans.cc
cov_ener_SOURCES = cov_ener.cc
diffus_SOURCES = diffus.cc
tstrip_SOURCES = tstrip.cc
gathtraj_SOURCES = gathtraj.cc
edyn_SOURCES = edyn.cc
rmsf_SOURCES = rmsf.cc
dipole_SOURCES = dipole.cc
cluster_SOURCES = cluster.cc
rmsdmat_SOURCES = rmsdmat.cc
rmsdmat2_SOURCES = rmsdmat2.cc
cluster2_SOURCES = cluster2.cc
postcluster_SOURCES = postcluster.cc
checkbox_SOURCES = checkbox.cc
maketop_SOURCES = maketop.cc
contop_SOURCES = contop.cc
comtop_SOURCES = comtop.cc
redtop_SOURCES = redtop.cc
pt_top_SOURCES = pt_top.cc
m_pt_top_SOURCES = m_pt_top.cc
LB_top_SOURCES = LB_top.cc
copy_box_SOURCES = copy_box.cc
buildbox_SOURCES = buildbox.cc
bin_box_SOURCES = bin_box.cc
ran_box_SOURCES = ran_box.cc
ran_solvation_SOURCES = ran_solvation.cc
perttop_SOURCES = perttop.cc
pert_top_SOURCES = pert_top.cc
explode_SOURCES = explode.cc
rdf_SOURCES = rdf.cc
mdf_SOURCES = mdf.cc
iwdcf_SOURCES = iwdcf.cc
neighbours_SOURCES = neighbours.cc
clusmol_SOURCES = clusmol.cc
methanes_SOURCES = methanes.cc
mcti_SOURCES = mcti.cc
iondens_SOURCES = iondens.cc
mkscript_SOURCES = mkscript.h mkscript.cc
checktopo_SOURCES = checktopo.cc
check_top_SOURCES = check_top.cc
ene_ana_SOURCES = ene_ana.cc
visco_ana_SOURCES = visco_ana.cc
epsilon_SOURCES = epsilon.cc
progch_SOURCES = progch.cc
probox_SOURCES = probox.cc
proion_SOURCES = proion.cc
progca_SOURCES = progca.cc
progca_series_SOURCES = progca_series.cc
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
proarse_SOURCES = proarse.cc
shake_SOURCES = shake.cc
protcf_SOURCES = protcf.cc
m_widom_SOURCES = m_widom.cc
cog_SOURCES = cog.cc
rot_rel_SOURCES = rot_rel.cc
inbox_SOURCES = inbox.cc
procry_SOURCES = procry.cc
unify_box_SOURCES = unify_box.cc
monoconv_SOURCES = monoconv.cc

LDADD = $(top_builddir)/src/libgromos.la
