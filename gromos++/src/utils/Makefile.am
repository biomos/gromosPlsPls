
noinst_LTLIBRARIES = libutils.la

gincludedir = $(includedir)/gromos++

ginclude_HEADERS = Rmsd.h\
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
                      Hbond3c.h\
                      Hbondcalc.h\
		      Dssp.h\
                      EnergyTraj.h\
                      SimplePairlist.h\
                      FfExpert.h\
					Cluster.h\
					maketop.h

libutils_la_SOURCES = Rmsd.cc \
                      Neighbours.cc \
                      VirtualAtom.cc \
                      Noe.cc \
                      AtomSpecifier.cc \
                      PropertyContainer.cc \
	   	      Property.cc \
                      Energy.cc \
                      TrajArray.cc \
                      RmsdMat.cc \
                      CheckTopo.cc \
                      Hbond.cc \
                      Hbond3c.cc \
                      Hbondcalc.cc \
		      Dssp.cc \
                      EnergyTraj.cc \
                      SimplePairlist.cc \
                      FfExpert.cc

check_PROGRAMS = Rmsd \
                  Neighbours \
                  VirtualAtom \
                  Noe \
                  AtomSpecifier \
		  PropertyContainer \
                  RmsdMat \
                  Energy \
                  CheckTopo \
                  SimplePairlist \
                  FfExpert


LDADD = ../libgromos.la

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
FfExpert_SOURCES = FfExpert.t.cc

AM_LDFLAGS = $(GSL_LDFLAGS)

VirtualAtom_LDADD = libutils.la \
        ../bound/libbound.la \
        ../fit/libfit.la \
        ../gio/libgio.la \
        ../gcore/libgcore.la \
        ../args/libargs.la \
        ../gmath/libgmath.la \
	$(GSL_LIB)

FfExpert_LDADD = libutils.la \
        ../bound/libbound.la \
        ../fit/libfit.la \
        ../gio/libgio.la \
        ../gcore/libgcore.la \
        ../args/libargs.la \
        ../gmath/libgmath.la \
	$(GSL_LIB)

