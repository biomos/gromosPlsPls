# This is the programs automake template.

INCLUDES = -I../src

bin_PROGRAMS = rmsd\
               pdb2g96\
               ener\
               pairener\
               dg_ener\
               convert\
               noe\
               rgyr\
               vmdam\
               tser\
               diffus\
	       tser2\
               tstrip\
               edyn\
               rmsf\
               dipole\
               cluster\
               checkbox\
               contop\
               comtop\
               redtop\
               copy_box\
               buildbox\
               perttop\
               explode\
               Hbondintra\
               Hbondinter\
               dist\
	       dist2\
               frameout\
               oparam\
               rdf\
               neighbours\
               clusmol\
	       propertyrmsd\
               cluster_analyze

oparam_SOURCES = oparam.cc              
frameout_SOURCES = frameout.cc
rmsd_SOURCES = rmsd.cc
pdb2g96_SOURCES = pdb2g96.cc
ener_SOURCES = ener.cc
pairener_SOURCES = pairener.cc
dg_ener_SOURCES = dg_ener.cc
Hbondintra_SOURCES = Hbondintra.cc
Hbondinter_SOURCES = Hbondinter.cc
convert_SOURCES = convert.cc
noe_SOURCES = noe.cc
rgyr_SOURCES = rgyr.cc
vmdam_SOURCES = vmdam.cc
tser_SOURCES = tser.cc
diffus_SOURCES = diffus.cc
tser2_SOURCES = tser2.cc
tstrip_SOURCES = tstrip.cc
edyn_SOURCES = edyn.cc
rmsf_SOURCES = rmsf.cc
dipole_SOURCES = dipole.cc
cluster_SOURCES = cluster.cc
checkbox_SOURCES = checkbox.cc
contop_SOURCES = contop.cc
comtop_SOURCES = comtop.cc
redtop_SOURCES = redtop.cc
copy_box_SOURCES = copy_box.cc
buildbox_SOURCES = buildbox.cc
perttop_SOURCES = perttop.cc
explode_SOURCES = explode.cc
cluster_analyze_SOURCES = cluster_analyze.cc
dist_SOURCES = dist.cc
dist2_SOURCES = dist2.cc
rdf_SOURCES = rdf.cc
neighbours_SOURCES = neighbours.cc
clusmol_SOURCES = clusmol.cc
propertyrmsd_SOURCES = propertyrmsd.cc

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

tser2_LDADD = $(top_builddir)/src/args/libargs.la\
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

tstrip_LDADD = $(top_builddir)/src/args/libargs.la\
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

cluster_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

cluster_analyze_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

checkbox_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la\
        $(top_builddir)/src/fit/libfit.la\
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

dist_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/fit/libfit.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/gmath/libgmath.la

dist2_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/fit/libfit.la \
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

Hbondintra_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/gcore/libgcore.la\
        $(top_builddir)/src/utils/libutils.la\
        $(top_builddir)/src/gmath/libgmath.la

Hbondinter_LDADD = $(top_builddir)/src/args/libargs.la\
        $(top_builddir)/src/bound/libbound.la \
        $(top_builddir)/src/gio/libgio.la\
        $(top_builddir)/src/utils/libutils.la \
        $(top_builddir)/src/gcore/libgcore.la\
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
