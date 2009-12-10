AM_CXXFLAGS = @MY_CXXFLAGS@
noinst_LTLIBRARIES = libutils.la

gincludedir = $(includedir)/gromos++

ginclude_HEADERS =    parse.h \
		      Rmsd.h\
                      Neighbours.h\
                      VirtualAtom.h\
                      Noe.h\
                      AtomSpecifier.h\
                      VectorSpecifier.h\
                      PropertyContainer.h\
		      Property.h\
                      Energy.h\
                      TrajArray.h\
                      CheckTopo.h\
                      Hbond.h\
                      Hbond3c.h\
                      Hbondcalc.h\
		      Dssp.h\
                      EnergyTraj.h\
                      RestrTraj.h\
                      SimplePairlist.h\
                      FfExpert.h\
                      FfExpertGraph.h\
		      maketop.h \
                      make_top.h \
		      Value.h \
		      ExpressionParser.h \
		      ExpressionParser.cc \
		      Time.h \
		      RdcFuncs.h \
		      JvalFuncs.h

libutils_la_SOURCES = parse.cc \
		      Rmsd.cc \
                      Neighbours.cc \
                      VirtualAtom.cc \
                      Noe.cc \
                      AtomSpecifier.cc \
                      AtomSpecifier_parse.cc \
                      VectorSpecifier.cc \
                      PropertyContainer.cc \
	   	      Property.cc \
                      Energy.cc \
                      TrajArray.cc \
                      CheckTopo.cc \
                      Hbond.cc \
                      Hbond3c.cc \
                      Hbondcalc.cc \
		      Dssp.cc \
                      EnergyTraj.cc \
                      RestrTraj.cc \
                      SimplePairlist.cc \
                      FfExpert.cc\
                      FfExpertGraph.cc\
                      Time.cc\
		      RdcFuncs.cc\
		      JvalFuncs.cc

check_PROGRAMS =  ExpressionParser \
		  Rmsd \
                  Neighbours \
                  VirtualAtom \
                  Noe \
                  AtomSpecifier \
		  PropertyContainer \
                  Energy \
                  CheckTopo \
                  SimplePairlist \
                  FfExpert


LDADD = ../libgromos.la

ExpressionParser_SOURCES = ExpressionParser.t.cc
Rmsd_SOURCES = Rmsd.t.cc
Neighbours_SOURCES = Neighbours.t.cc
VirtualAtom_SOURCES = VirtualAtom.t.cc
Noe_SOURCES = Noe.t.cc
AtomSpecifier_SOURCES = AtomSpecifier.t.cc
PropertyContainer_SOURCES = PropertyContainer.t.cc
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

ExpressionParser_LDADD =  libutils.la \
        ../gio/libgio.la \
        ../gcore/libgcore.la \
        ../args/libargs.la \
        ../gmath/libgmath.la \
	../bound/libbound.la \
        ../fit/libfit.la \
	$(GSL_LIB)
