## Template for automake for libutils

lib_LTLIBRARIES = libutils.la
libutils_la_SOURCES = Rmsd.cc\
                      Neighbours.cc\
                      VirtualAtom.cc\
                      Noe.cc\
                      AtomSpecifier.cc\
                      PropertyContainer.cc\
		      Property.cc\
                      Energy.cc

libutils_la_LDFLAGS = -version-info $(GROMOS_MAJOR_VERSION):$(GROMOS_MICRO_VERSION):0

noinst_PROGRAMS = Rmsd\
                  Neighbours\
                  VirtualAtom\
                  Noe\
                  AtomSpecifier\
		  PropertyContainer\
                  Energy


LDADD = libutils.la\
        ../fit/libfit.la\
        ../gio/libgio.la\
        ../gcore/libgcore.la\
        ../gmath/libgmath.la\
        ../bound/libbound.la

Rmsd_SOURCES = Rmsd.t.cc
Neighbours_SOURCES = Neighbours.t.cc
VirtualAtom_SOURCES = VirtualAtom.t.cc
Noe_SOURCES = Noe.t.cc
AtomSpecifier_SOURCES = AtomSpecifier.t.cc
PropertyContainer_SOURCES = PropertyContainer.t.cc
Energy_SOURCES = Energy.t.cc
