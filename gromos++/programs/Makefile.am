# This is the programs automake template.

INCLUDES = -I@abs_top_srcdir@/src

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
               checkbox\
               maketop\
               contop\
               comtop\
               redtop\
               pt_top\
               copy_box\
               buildbox\
               bin_box\
               perttop\
               explode\
               hbond\
               dist\
               frameout\
               oparam\
               nhoparam\
               rdf\
               iwdcf\
               neighbours\
               clusmol\
	       propertyrmsd\
               mcti\
               iondens\
               mkscript\
               checktopo\
               ene_ana\
               epsilon\
               progch\
               probox\
               proion\
               progca\
	       dssp\
               atominfo\
               shake_analysis\
               procmt_list\
               filter\
	       follow\
               renumber\
               countbb\
               jval\
               prepbbb

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
checkbox_SOURCES = checkbox.cc
maketop_SOURCES = maketop.cc
contop_SOURCES = contop.cc
comtop_SOURCES = comtop.cc
redtop_SOURCES = redtop.cc
pt_top_SOURCES = pt_top.cc
copy_box_SOURCES = copy_box.cc
buildbox_SOURCES = buildbox.cc
bin_box_SOURCES = bin_box.cc
perttop_SOURCES = perttop.cc
explode_SOURCES = explode.cc
dist_SOURCES = dist.cc
rdf_SOURCES = rdf.cc
iwdcf_SOURCES = iwdcf.cc
neighbours_SOURCES = neighbours.cc
clusmol_SOURCES = clusmol.cc
propertyrmsd_SOURCES = propertyrmsd.cc
mcti_SOURCES = mcti.cc
iondens_SOURCES = iondens.cc
mkscript_SOURCES = mkscript.cc
checktopo_SOURCES = checktopo.cc
ene_ana_SOURCES = ene_ana.cc
epsilon_SOURCES = epsilon.cc
progch_SOURCES = progch.cc
probox_SOURCES = probox.cc
proion_SOURCES = proion.cc
progca_SOURCES = progca.cc
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

rmsd_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/fit/libfit.la \
	$(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

convert_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/fit/libfit.la \
	$(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

noe_LDADD = $(top_builddir)/src/args/libargs.la\
            $(top_builddir)/src/utils/libutils.la \
            $(top_builddir)/src/bound/libbound.la \
            $(top_builddir)/src/fit/libfit.la \
	    $(top_builddir)/src/gio/libgio.la\
            $(top_builddir)/src/gcore/libgcore.la\
            $(top_builddir)/src/gmath/libgmath.la

prepnoe_LDADD = $(top_builddir)/src/args/libargs.la\
            $(top_builddir)/src/utils/libutils.la \
            $(top_builddir)/src/bound/libbound.la \
            $(top_builddir)/src/fit/libfit.la \
            $(top_builddir)/src/gio/libgio.la\
            $(top_builddir)/src/gcore/libgcore.la\
            $(top_builddir)/src/gmath/libgmath.la

noenew_LDADD = $(top_builddir)/src/args/libargs.la\
            $(top_builddir)/src/utils/libutils.la \
            $(top_builddir)/src/bound/libbound.la \
            $(top_builddir)/src/fit/libfit.la \
            $(top_builddir)/src/gio/libgio.la\
            $(top_builddir)/src/gcore/libgcore.la\
            $(top_builddir)/src/gmath/libgmath.la

prepnoenew_LDADD = $(top_builddir)/src/args/libargs.la\
            $(top_builddir)/src/utils/libutils.la \
            $(top_builddir)/src/bound/libbound.la \
            $(top_builddir)/src/fit/libfit.la \
            $(top_builddir)/src/gio/libgio.la\
            $(top_builddir)/src/gcore/libgcore.la\
            $(top_builddir)/src/gmath/libgmath.la

postnoe_LDADD = $(top_builddir)/src/args/libargs.la\
            $(top_builddir)/src/utils/libutils.la \
            $(top_builddir)/src/bound/libbound.la \
            $(top_builddir)/src/fit/libfit.la \
            $(top_builddir)/src/gio/libgio.la\
            $(top_builddir)/src/gcore/libgcore.la\
            $(top_builddir)/src/gmath/libgmath.la

rgyr_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

vmdam_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

tser_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la

ditrans_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la

diffus_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la

propertyrmsd_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la

cov_ener_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la

tstrip_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/fit/libfit.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gmath/libgmath.la

gathtraj_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/fit/libfit.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

edyn_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

espmap_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

rmsf_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

dipole_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

rmsdmat_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/bound/libbound.la


cluster_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/utils/libutils.la\
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/bound/libbound.la

checkbox_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la\
        $(top_builddir)/src/fit/libfit.la\
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/utils/libutils.la

dist_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la

contop_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/bound/libbound.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la

comtop_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/bound/libbound.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la


redtop_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/bound/libbound.la\
        $(top_builddir)/src/fit/libfit.la\
        $(top_builddir)/src/gmath/libgmath.la

pt_top_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/bound/libbound.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la

maketop_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/fit/libfit.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/bound/libbound.la

copy_box_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

bin_box_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/bound/libbound.la

buildbox_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/bound/libbound.la

perttop_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/bound/libbound.la

explode_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

hbond_LDADD = $(top_builddir)/src/utils/libutils.la\
		$(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la

frameout_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la

oparam_LDADD = $(top_builddir)/src/args/libargs.la\
      $(top_builddir)/src/bound/libbound.la \
      $(top_builddir)/src/gio/libgio.la\
      $(top_builddir)/src/fit/libfit.la\
      $(top_builddir)/src/gcore/libgcore.la\
      $(top_builddir)/src/gmath/libgmath.la

nhoparam_LDADD = $(top_builddir)/src/args/libargs.la\
      $(top_builddir)/src/bound/libbound.la \
      $(top_builddir)/src/gio/libgio.la\
      $(top_builddir)/src/gcore/libgcore.la\
      $(top_builddir)/src/utils/libutils.la \
      $(top_builddir)/src/fit/libfit.la\
      $(top_builddir)/src/gmath/libgmath.la

dg_ener_LDADD = $(top_builddir)/src/args/libargs.la\
                $(top_builddir)/src/bound/libbound.la \
                $(top_builddir)/src/utils/libutils.la \
                $(top_builddir)/src/gio/libgio.la\
                $(top_builddir)/src/gcore/libgcore.la\
                $(top_builddir)/src/gmath/libgmath.la\
                $(top_builddir)/src/fit/libfit.la

pairener_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la

ener_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la

nb_ener_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la

pdb2g96_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la

rdf_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la

iwdcf_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la

neighbours_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

clusmol_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/utils/libutils.la\
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/bound/libbound.la\
        $(top_builddir)/src/fit/libfit.la

mcti_LDADD = $(top_builddir)/src/args/libargs.la\
             $(top_builddir)/src/utils/libutils.la\
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/bound/libbound.la\
        $(top_builddir)/src/fit/libfit.la

iondens_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

mkscript_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

checktopo_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

ene_ana_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

epsilon_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

progch_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

probox_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

proion_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

progca_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

dssp_LDADD = $(top_builddir)/src/utils/libutils.la\
		$(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la

atominfo_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/utils/libutils.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la

shake_analysis_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/utils/libutils.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la

procmt_list_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/utils/libutils.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la

filter_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/utils/libutils.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la

follow_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/utils/libutils.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la

renumber_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/utils/libutils.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la

countbb_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/utils/libutils.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la

jval_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/utils/libutils.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la

prepbbb_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/utils/libutils.la\
        $(top_builddir)/src/gmath/libgmath.la\
        $(top_builddir)/src/fit/libfit.la
