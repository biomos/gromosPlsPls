# Makefile.in generated automatically by automake 1.4-p4 from Makefile.am

# Copyright (C) 1994, 1995-8, 1999 Free Software Foundation, Inc.
# This Makefile.in is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.


SHELL = /bin/sh

srcdir = .
top_srcdir = ../..

prefix = /home/mika/igc4/g++/CVS/gromos++
exec_prefix = ${prefix}

bindir = ${exec_prefix}/bin
sbindir = ${exec_prefix}/sbin
libexecdir = ${exec_prefix}/libexec
datadir = ${prefix}/share
sysconfdir = ${prefix}/etc
sharedstatedir = ${prefix}/com
localstatedir = ${prefix}/var
libdir = ${exec_prefix}/lib
infodir = ${prefix}/info
mandir = ${prefix}/man
includedir = ${prefix}/include
oldincludedir = /usr/include

DESTDIR =

pkgdatadir = $(datadir)/gromos
pkglibdir = $(libdir)/gromos
pkgincludedir = $(includedir)/gromos

top_builddir = ../..

ACLOCAL = aclocal
AUTOCONF = autoconf
AUTOMAKE = automake
AUTOHEADER = autoheader

INSTALL = /usr/local/bin/install -c
INSTALL_PROGRAM = ${INSTALL} $(AM_INSTALL_PROGRAM_FLAGS)
INSTALL_DATA = ${INSTALL} -m 644
INSTALL_SCRIPT = ${INSTALL}
transform = s,x,x,

NORMAL_INSTALL = :
PRE_INSTALL = :
POST_INSTALL = :
NORMAL_UNINSTALL = :
PRE_UNINSTALL = :
POST_UNINSTALL = :
host_alias = @host_alias@
host_triplet = sparc-sun-solaris2.8
AS = @AS@
CC = gcc
CXX = g++
DLLTOOL = @DLLTOOL@
ECHO = /usr/ucb/echo
EXEEXT = 
GROMOS_MAJOR_VERSION = 0
GROMOS_MICRO_VERSION = 1
GROMOS_MINOR_VERSION = 1
GROMOS_VERSION = 0.1.1
LIBTOOL = $(SHELL) $(top_builddir)/libtool
LN_S = ln -s
MAINT = #
MAKEINFO = makeinfo
OBJDUMP = @OBJDUMP@
OBJEXT = o
PACKAGE = gromos
RANLIB = ranlib
STRIP = strip
VERSION = 0.1.1

lib_LTLIBRARIES = libgcore.la

libgcore_la_LDFLAGS = -version-info $(GROMOS_MAJOR_VERSION):$(GROMOS_MICRO_VERSION):0

# noinst_LIBRARIES = libgcore.la

libgcore_la_SOURCES = Angle.cc                      AngleType.cc                      AtomPair.cc                      AtomTopology.cc                      Bond.cc                      BondType.cc                      Box.cc                      Constraint.cc                      Dihedral.cc                      DihedralType.cc                      Exclusion.cc                      GromosForceField.cc                      Improper.cc                      ImproperType.cc                      LJType.cc                      MoleculeTopology.cc                      Molecule.cc                      SolventTopology.cc                      Solvent.cc                      System.cc


bin_PROGRAMS = AtomTopology                Exclusion                MoleculeTopology                Molecule                SolventTopology                Solvent


LDADD = libgcore.la

AtomTopology_SOURCES = AtomTopology.t.cc
Exclusion_SOURCES = Exclusion.t.cc
MoleculeTopology_SOURCES = MoleculeTopology.t.cc
Molecule_SOURCES = Molecule.t.cc
SolventTopology_SOURCES = SolventTopology.t.cc
Solvent_SOURCES = Solvent.t.cc
mkinstalldirs = $(SHELL) $(top_srcdir)/mkinstalldirs
CONFIG_CLEAN_FILES = 
LTLIBRARIES =  $(lib_LTLIBRARIES)


DEFS = -DHAVE_DLFCN_H=1 -DSTDC_HEADERS=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_UNISTD_H=1 -DHAVE_LIBM=1  -I. -I$(srcdir) 
CPPFLAGS = 
LDFLAGS = 
LIBS = -lm 
libgcore_la_LIBADD = 
libgcore_la_OBJECTS =  Angle.lo AngleType.lo AtomPair.lo AtomTopology.lo \
Bond.lo BondType.lo Box.lo Constraint.lo Dihedral.lo DihedralType.lo \
Exclusion.lo GromosForceField.lo Improper.lo ImproperType.lo LJType.lo \
MoleculeTopology.lo Molecule.lo SolventTopology.lo Solvent.lo System.lo
bin_PROGRAMS =  AtomTopology$(EXEEXT) Exclusion$(EXEEXT) \
MoleculeTopology$(EXEEXT) Molecule$(EXEEXT) SolventTopology$(EXEEXT) \
Solvent$(EXEEXT)
PROGRAMS =  $(bin_PROGRAMS)

AtomTopology_OBJECTS =  AtomTopology.t.$(OBJEXT)
AtomTopology_LDADD = $(LDADD)
AtomTopology_DEPENDENCIES =  libgcore.la
AtomTopology_LDFLAGS = 
Exclusion_OBJECTS =  Exclusion.t.$(OBJEXT)
Exclusion_LDADD = $(LDADD)
Exclusion_DEPENDENCIES =  libgcore.la
Exclusion_LDFLAGS = 
MoleculeTopology_OBJECTS =  MoleculeTopology.t.$(OBJEXT)
MoleculeTopology_LDADD = $(LDADD)
MoleculeTopology_DEPENDENCIES =  libgcore.la
MoleculeTopology_LDFLAGS = 
Molecule_OBJECTS =  Molecule.t.$(OBJEXT)
Molecule_LDADD = $(LDADD)
Molecule_DEPENDENCIES =  libgcore.la
Molecule_LDFLAGS = 
SolventTopology_OBJECTS =  SolventTopology.t.$(OBJEXT)
SolventTopology_LDADD = $(LDADD)
SolventTopology_DEPENDENCIES =  libgcore.la
SolventTopology_LDFLAGS = 
Solvent_OBJECTS =  Solvent.t.$(OBJEXT)
Solvent_LDADD = $(LDADD)
Solvent_DEPENDENCIES =  libgcore.la
Solvent_LDFLAGS = 
CXXFLAGS = -g -O2 -Wall
CXXCOMPILE = $(CXX) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS)
LTCXXCOMPILE = $(LIBTOOL) --mode=compile $(CXX) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS)
CXXLD = $(CXX)
CXXLINK = $(LIBTOOL) --mode=link $(CXXLD) $(AM_CXXFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@
DIST_COMMON =  Makefile.am Makefile.in


DISTFILES = $(DIST_COMMON) $(SOURCES) $(HEADERS) $(TEXINFOS) $(EXTRA_DIST)

TAR = tar
GZIP_ENV = --best
DEP_FILES =  .deps/Angle.P .deps/AngleType.P .deps/AtomPair.P \
.deps/AtomTopology.P .deps/AtomTopology.t.P .deps/Bond.P \
.deps/BondType.P .deps/Box.P .deps/Constraint.P .deps/Dihedral.P \
.deps/DihedralType.P .deps/Exclusion.P .deps/Exclusion.t.P \
.deps/GromosForceField.P .deps/Improper.P .deps/ImproperType.P \
.deps/LJType.P .deps/Molecule.P .deps/Molecule.t.P \
.deps/MoleculeTopology.P .deps/MoleculeTopology.t.P .deps/Solvent.P \
.deps/Solvent.t.P .deps/SolventTopology.P .deps/SolventTopology.t.P \
.deps/System.P
SOURCES = $(libgcore_la_SOURCES) $(AtomTopology_SOURCES) $(Exclusion_SOURCES) $(MoleculeTopology_SOURCES) $(Molecule_SOURCES) $(SolventTopology_SOURCES) $(Solvent_SOURCES)
OBJECTS = $(libgcore_la_OBJECTS) $(AtomTopology_OBJECTS) $(Exclusion_OBJECTS) $(MoleculeTopology_OBJECTS) $(Molecule_OBJECTS) $(SolventTopology_OBJECTS) $(Solvent_OBJECTS)

all: all-redirect
.SUFFIXES:
.SUFFIXES: .S .c .cc .lo .o .obj .s
$(srcdir)/Makefile.in: # Makefile.am $(top_srcdir)/configure.in $(ACLOCAL_M4) 
	cd $(top_srcdir) && $(AUTOMAKE) --gnu src/gcore/Makefile

Makefile: $(srcdir)/Makefile.in  $(top_builddir)/config.status $(BUILT_SOURCES)
	cd $(top_builddir) \
	  && CONFIG_FILES=$(subdir)/$@ CONFIG_HEADERS= $(SHELL) ./config.status


mostlyclean-libLTLIBRARIES:

clean-libLTLIBRARIES:
	-test -z "$(lib_LTLIBRARIES)" || rm -f $(lib_LTLIBRARIES)

distclean-libLTLIBRARIES:

maintainer-clean-libLTLIBRARIES:

install-libLTLIBRARIES: $(lib_LTLIBRARIES)
	@$(NORMAL_INSTALL)
	$(mkinstalldirs) $(DESTDIR)$(libdir)
	@list='$(lib_LTLIBRARIES)'; for p in $$list; do \
	  if test -f $$p; then \
	    echo "$(LIBTOOL)  --mode=install $(INSTALL) $$p $(DESTDIR)$(libdir)/$$p"; \
	    $(LIBTOOL)  --mode=install $(INSTALL) $$p $(DESTDIR)$(libdir)/$$p; \
	  else :; fi; \
	done

uninstall-libLTLIBRARIES:
	@$(NORMAL_UNINSTALL)
	list='$(lib_LTLIBRARIES)'; for p in $$list; do \
	  $(LIBTOOL)  --mode=uninstall rm -f $(DESTDIR)$(libdir)/$$p; \
	done

# FIXME: We should only use cygpath when building on Windows,
# and only if it is available.
.c.obj:
	$(COMPILE) -c `cygpath -w $<`

.s.o:
	$(COMPILE) -c $<

.S.o:
	$(COMPILE) -c $<

mostlyclean-compile:
	-rm -f *.o core *.core
	-rm -f *.$(OBJEXT)

clean-compile:

distclean-compile:
	-rm -f *.tab.c

maintainer-clean-compile:

.s.lo:
	$(LIBTOOL) --mode=compile $(COMPILE) -c $<

.S.lo:
	$(LIBTOOL) --mode=compile $(COMPILE) -c $<

mostlyclean-libtool:
	-rm -f *.lo

clean-libtool:
	-rm -rf .libs _libs

distclean-libtool:

maintainer-clean-libtool:

libgcore.la: $(libgcore_la_OBJECTS) $(libgcore_la_DEPENDENCIES)
	$(CXXLINK) -rpath $(libdir) $(libgcore_la_LDFLAGS) $(libgcore_la_OBJECTS) $(libgcore_la_LIBADD) $(LIBS)

mostlyclean-binPROGRAMS:

clean-binPROGRAMS:
	-test -z "$(bin_PROGRAMS)" || rm -f $(bin_PROGRAMS)

distclean-binPROGRAMS:

maintainer-clean-binPROGRAMS:

install-binPROGRAMS: $(bin_PROGRAMS)
	@$(NORMAL_INSTALL)
	$(mkinstalldirs) $(DESTDIR)$(bindir)
	@list='$(bin_PROGRAMS)'; for p in $$list; do \
	  if test -f $$p; then \
	    echo " $(LIBTOOL)  --mode=install $(INSTALL_PROGRAM) $$p $(DESTDIR)$(bindir)/`echo $$p|sed 's/$(EXEEXT)$$//'|sed '$(transform)'|sed 's/$$/$(EXEEXT)/'`"; \
	    $(LIBTOOL)  --mode=install $(INSTALL_PROGRAM) $$p $(DESTDIR)$(bindir)/`echo $$p|sed 's/$(EXEEXT)$$//'|sed '$(transform)'|sed 's/$$/$(EXEEXT)/'`; \
	  else :; fi; \
	done

uninstall-binPROGRAMS:
	@$(NORMAL_UNINSTALL)
	list='$(bin_PROGRAMS)'; for p in $$list; do \
	  rm -f $(DESTDIR)$(bindir)/`echo $$p|sed 's/$(EXEEXT)$$//'|sed '$(transform)'|sed 's/$$/$(EXEEXT)/'`; \
	done

AtomTopology$(EXEEXT): $(AtomTopology_OBJECTS) $(AtomTopology_DEPENDENCIES)
	@rm -f AtomTopology$(EXEEXT)
	$(CXXLINK) $(AtomTopology_LDFLAGS) $(AtomTopology_OBJECTS) $(AtomTopology_LDADD) $(LIBS)

Exclusion$(EXEEXT): $(Exclusion_OBJECTS) $(Exclusion_DEPENDENCIES)
	@rm -f Exclusion$(EXEEXT)
	$(CXXLINK) $(Exclusion_LDFLAGS) $(Exclusion_OBJECTS) $(Exclusion_LDADD) $(LIBS)

MoleculeTopology$(EXEEXT): $(MoleculeTopology_OBJECTS) $(MoleculeTopology_DEPENDENCIES)
	@rm -f MoleculeTopology$(EXEEXT)
	$(CXXLINK) $(MoleculeTopology_LDFLAGS) $(MoleculeTopology_OBJECTS) $(MoleculeTopology_LDADD) $(LIBS)

Molecule$(EXEEXT): $(Molecule_OBJECTS) $(Molecule_DEPENDENCIES)
	@rm -f Molecule$(EXEEXT)
	$(CXXLINK) $(Molecule_LDFLAGS) $(Molecule_OBJECTS) $(Molecule_LDADD) $(LIBS)

SolventTopology$(EXEEXT): $(SolventTopology_OBJECTS) $(SolventTopology_DEPENDENCIES)
	@rm -f SolventTopology$(EXEEXT)
	$(CXXLINK) $(SolventTopology_LDFLAGS) $(SolventTopology_OBJECTS) $(SolventTopology_LDADD) $(LIBS)

Solvent$(EXEEXT): $(Solvent_OBJECTS) $(Solvent_DEPENDENCIES)
	@rm -f Solvent$(EXEEXT)
	$(CXXLINK) $(Solvent_LDFLAGS) $(Solvent_OBJECTS) $(Solvent_LDADD) $(LIBS)
.cc.o:
	$(CXXCOMPILE) -c $<
.cc.obj:
	$(CXXCOMPILE) -c `cygpath -w $<`
.cc.lo:
	$(LTCXXCOMPILE) -c $<

tags: TAGS

ID: $(HEADERS) $(SOURCES) $(LISP)
	list='$(SOURCES) $(HEADERS)'; \
	unique=`for i in $$list; do echo $$i; done | \
	  awk '    { files[$$0] = 1; } \
	       END { for (i in files) print i; }'`; \
	here=`pwd` && cd $(srcdir) \
	  && mkid -f$$here/ID $$unique $(LISP)

TAGS:  $(HEADERS) $(SOURCES)  $(TAGS_DEPENDENCIES) $(LISP)
	tags=; \
	here=`pwd`; \
	list='$(SOURCES) $(HEADERS)'; \
	unique=`for i in $$list; do echo $$i; done | \
	  awk '    { files[$$0] = 1; } \
	       END { for (i in files) print i; }'`; \
	test -z "$(ETAGS_ARGS)$$unique$(LISP)$$tags" \
	  || (cd $(srcdir) && etags $(ETAGS_ARGS) $$tags  $$unique $(LISP) -o $$here/TAGS)

mostlyclean-tags:

clean-tags:

distclean-tags:
	-rm -f TAGS ID

maintainer-clean-tags:

distdir = $(top_builddir)/$(PACKAGE)-$(VERSION)/$(subdir)

subdir = src/gcore

distdir: $(DISTFILES)
	here=`cd $(top_builddir) && pwd`; \
	top_distdir=`cd $(top_distdir) && pwd`; \
	distdir=`cd $(distdir) && pwd`; \
	cd $(top_srcdir) \
	  && $(AUTOMAKE) --include-deps --build-dir=$$here --srcdir-name=$(top_srcdir) --output-dir=$$top_distdir --gnu src/gcore/Makefile
	@for file in $(DISTFILES); do \
	  d=$(srcdir); \
	  if test -d $$d/$$file; then \
	    cp -pr $$d/$$file $(distdir)/$$file; \
	  else \
	    test -f $(distdir)/$$file \
	    || ln $$d/$$file $(distdir)/$$file 2> /dev/null \
	    || cp -p $$d/$$file $(distdir)/$$file || :; \
	  fi; \
	done

DEPS_MAGIC := $(shell mkdir .deps > /dev/null 2>&1 || :)

-include $(DEP_FILES)

mostlyclean-depend:

clean-depend:

distclean-depend:
	-rm -rf .deps

maintainer-clean-depend:

%.o: %.c
	@echo '$(COMPILE) -c $<'; \
	$(COMPILE) -Wp,-MD,.deps/$(*F).pp -c $<
	@-cp .deps/$(*F).pp .deps/$(*F).P; \
	tr ' ' '\012' < .deps/$(*F).pp \
	  | sed -e 's/^\\$$//' -e '/^$$/ d' -e '/:$$/ d' -e 's/$$/ :/' \
	    >> .deps/$(*F).P; \
	rm .deps/$(*F).pp

%.lo: %.c
	@echo '$(LTCOMPILE) -c $<'; \
	$(LTCOMPILE) -Wp,-MD,.deps/$(*F).pp -c $<
	@-sed -e 's/^\([^:]*\)\.o[ 	]*:/\1.lo \1.o :/' \
	  < .deps/$(*F).pp > .deps/$(*F).P; \
	tr ' ' '\012' < .deps/$(*F).pp \
	  | sed -e 's/^\\$$//' -e '/^$$/ d' -e '/:$$/ d' -e 's/$$/ :/' \
	    >> .deps/$(*F).P; \
	rm -f .deps/$(*F).pp

%.o: %.cc
	@echo '$(CXXCOMPILE) -c $<'; \
	$(CXXCOMPILE) -Wp,-MD,.deps/$(*F).pp -c $<
	@-cp .deps/$(*F).pp .deps/$(*F).P; \
	tr ' ' '\012' < .deps/$(*F).pp \
	  | sed -e 's/^\\$$//' -e '/^$$/ d' -e '/:$$/ d' -e 's/$$/ :/' \
	    >> .deps/$(*F).P; \
	rm .deps/$(*F).pp

%.lo: %.cc
	@echo '$(LTCXXCOMPILE) -c $<'; \
	$(LTCXXCOMPILE) -Wp,-MD,.deps/$(*F).pp -c $<
	@-sed -e 's/^\([^:]*\)\.o[ 	]*:/\1.lo \1.o :/' \
	  < .deps/$(*F).pp > .deps/$(*F).P; \
	tr ' ' '\012' < .deps/$(*F).pp \
	  | sed -e 's/^\\$$//' -e '/^$$/ d' -e '/:$$/ d' -e 's/$$/ :/' \
	    >> .deps/$(*F).P; \
	rm -f .deps/$(*F).pp
info-am:
info: info-am
dvi-am:
dvi: dvi-am
check-am: all-am
check: check-am
installcheck-am:
installcheck: installcheck-am
install-exec-am: install-libLTLIBRARIES install-binPROGRAMS
install-exec: install-exec-am

install-data-am:
install-data: install-data-am

install-am: all-am
	@$(MAKE) $(AM_MAKEFLAGS) install-exec-am install-data-am
install: install-am
uninstall-am: uninstall-libLTLIBRARIES uninstall-binPROGRAMS
uninstall: uninstall-am
all-am: Makefile $(LTLIBRARIES) $(PROGRAMS)
all-redirect: all-am
install-strip:
	$(MAKE) $(AM_MAKEFLAGS) AM_INSTALL_PROGRAM_FLAGS=-s install
installdirs:
	$(mkinstalldirs)  $(DESTDIR)$(libdir) $(DESTDIR)$(bindir)


mostlyclean-generic:

clean-generic:

distclean-generic:
	-rm -f Makefile $(CONFIG_CLEAN_FILES)
	-rm -f config.cache config.log stamp-h stamp-h[0-9]*

maintainer-clean-generic:
mostlyclean-am:  mostlyclean-libLTLIBRARIES mostlyclean-compile \
		mostlyclean-libtool mostlyclean-binPROGRAMS \
		mostlyclean-tags mostlyclean-depend mostlyclean-generic

mostlyclean: mostlyclean-am

clean-am:  clean-libLTLIBRARIES clean-compile clean-libtool \
		clean-binPROGRAMS clean-tags clean-depend clean-generic \
		mostlyclean-am

clean: clean-am

distclean-am:  distclean-libLTLIBRARIES distclean-compile \
		distclean-libtool distclean-binPROGRAMS distclean-tags \
		distclean-depend distclean-generic clean-am
	-rm -f libtool

distclean: distclean-am

maintainer-clean-am:  maintainer-clean-libLTLIBRARIES \
		maintainer-clean-compile maintainer-clean-libtool \
		maintainer-clean-binPROGRAMS maintainer-clean-tags \
		maintainer-clean-depend maintainer-clean-generic \
		distclean-am
	@echo "This command is intended for maintainers to use;"
	@echo "it deletes files that may require special tools to rebuild."

maintainer-clean: maintainer-clean-am

.PHONY: mostlyclean-libLTLIBRARIES distclean-libLTLIBRARIES \
clean-libLTLIBRARIES maintainer-clean-libLTLIBRARIES \
uninstall-libLTLIBRARIES install-libLTLIBRARIES mostlyclean-compile \
distclean-compile clean-compile maintainer-clean-compile \
mostlyclean-libtool distclean-libtool clean-libtool \
maintainer-clean-libtool mostlyclean-binPROGRAMS distclean-binPROGRAMS \
clean-binPROGRAMS maintainer-clean-binPROGRAMS uninstall-binPROGRAMS \
install-binPROGRAMS tags mostlyclean-tags distclean-tags clean-tags \
maintainer-clean-tags distdir mostlyclean-depend distclean-depend \
clean-depend maintainer-clean-depend info-am info dvi-am dvi check \
check-am installcheck-am installcheck install-exec-am install-exec \
install-data-am install-data install-am install uninstall-am uninstall \
all-redirect all-am all installdirs mostlyclean-generic \
distclean-generic clean-generic maintainer-clean-generic clean \
mostlyclean distclean maintainer-clean


# Tell versions [3.59,3.63) of GNU make to not export all variables.
# Otherwise a system limit (for SysV at least) may be exceeded.
.NOEXPORT:
