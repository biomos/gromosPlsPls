## Template for automake for libutils

lib_LTLIBRARIES = libutils.la

include_HEADERS = Rmsd.h\
                      Neighbours.h\
                      VirtualAtom.h\
                      Noe.h\
                      AtomSpecifier.h\
                      PropertyContainer.h\
		      Property.h\
                      Energy.h\
                      TrajArray.h\
                      RmsdMat.h\
                      CheckTopo.h\
                      Hbond.h\
                      Hbondcalc.h\
		      Dssp.h\
                      EnergyTraj.h\
                      SimplePairlist.h
libutils_la_SOURCES = Rmsd.cc\
                      Neighbours.cc\
                      VirtualAtom.cc\
                      Noe.cc\
                      AtomSpecifier.cc\
                      PropertyContainer.cc\
	   	      Property.cc\
                      Energy.cc\
                      TrajArray.cc\
                      RmsdMat.cc\
                      CheckTopo.cc\
                      Hbond.cc\
                      Hbondcalc.cc\
		      Dssp.cc\
                      EnergyTraj.cc\
                      SimplePairlist.cc

libutils_la_LDFLAGS = -version-info $(GROMOS_MAJOR_VERSION):$(GROMOS_MICRO_VERSION):0

check_PROGRAMS = Rmsd\
                  Neighbours\
                  VirtualAtom\
                  Noe\
                  AtomSpecifier\
		  PropertyContainer\
                  RmsdMat\
                  Energy\
                  CheckTopo\
                  SimplePairlist


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
RmsdMat_SOURCES = RmsdMat.t.cc
Energy_SOURCES = Energy.t.cc
CheckTopo_SOURCES = CheckTopo.t.cc
SimplePairlist_SOURCES = SimplePairlist.t.cc

VirtualAtom_LDADD = libutils.la\
        ../bound/libbound.la\
        ../fit/libfit.la\
        ../gio/libgio.la\
        ../gcore/libgcore.la\
        ../args/libargs.la\
        ../gmath/libgmath.la

