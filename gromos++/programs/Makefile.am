# This is the programs automake template.

INCLUDES = -I../src

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
               rgyr\
               vmdam\
               tser\
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
               perttop\
               explode\
               Hbond\
               dist\
               frameout\
               oparam\
               nhoparam\
               rdf\
               neighbours\
               clusmol\
	       propertyrmsd\
               mcti\
               iondens\
               mkscript\
               checktopo\
               ene_ana\
               ene_ana2\
               epsilon

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
Hbond_SOURCES = Hbond.cc
convert_SOURCES = convert.cc
noe_SOURCES = noe.cc
prepnoe_SOURCES = prepnoe.cc
rgyr_SOURCES = rgyr.cc
vmdam_SOURCES = vmdam.cc
tser_SOURCES = tser.cc
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
perttop_SOURCES = perttop.cc
explode_SOURCES = explode.cc
dist_SOURCES = dist.cc
rdf_SOURCES = rdf.cc
neighbours_SOURCES = neighbours.cc
clusmol_SOURCES = clusmol.cc
propertyrmsd_SOURCES = propertyrmsd.cc
mcti_SOURCES = mcti.cc
iondens_SOURCES = iondens.cc
mkscript_SOURCES = mkscript.cc
checktopo_SOURCES = checktopo.cc
ene_ana_SOURCES = ene_ana.cc
ene_ana2_SOURCES = ene_ana2.cc
epsilon_SOURCES = epsilon.cc

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

Hbond_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/utils/libutils.la\
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

ene_ana2_LDADD = $(top_builddir)/src/args/libargs.la\
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
