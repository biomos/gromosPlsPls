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
               rdf\
               neighbours\
               clusmol\
	       propertyrmsd\
               mcti

oparam_SOURCES = oparam.cc              
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

rgyr_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
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
        $(top_builddir)/src/gmath/libgmath.la

diffus_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

propertyrmsd_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

cov_ener_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

tstrip_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

gathtraj_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
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
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

rmsdmat_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la


cluster_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/utils/libutils.la

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
        $(top_builddir)/src/gmath/libgmath.la

contop_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la

comtop_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la

redtop_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/gcore/libgcore.la

pt_top_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/gcore/libgcore.la

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
        $(top_builddir)/src/gmath/libgmath.la

perttop_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la

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
        $(top_builddir)/src/gmath/libgmath.la

frameout_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

oparam_LDADD = $(top_builddir)/src/args/libargs.la\
      $(top_builddir)/src/bound/libbound.la \
      $(top_builddir)/src/gio/libgio.la\
      $(top_builddir)/src/gcore/libgcore.la\
      $(top_builddir)/src/gmath/libgmath.la

dg_ener_LDADD = $(top_builddir)/src/args/libargs.la

pairener_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

ener_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

nb_ener_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

pdb2g96_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

rdf_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

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
        $(top_builddir)/src/gcore/libgcore.la 

mcti_LDADD = $(top_builddir)/src/args/libargs.la
