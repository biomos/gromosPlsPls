## Template for automake for libgromos

lib_LTLIBRARIES = libgio.la

libgio_la_LDFLAGS = -version-info $(GROMOS_MAJOR_VERSION):$(GROMOS_MICRO_VERSION):0

include_HEADERS = Ginstream.h\
                   InTopology.h\
                   OutTopology.h\
                   InG96.h\
                   OutCoordinates.h\
                   OutG96.h\
                   OutG96S.h\
                   OutPdb.h\
                   Outvmdam.h\
                   InParameter.h\
                   InBuildingBlock.h

libgio_la_SOURCES = Ginstream.cc\
                   InTopology.cc\
                   OutTopology.cc\
                   InG96.cc\
                   OutCoordinates.cc\
                   OutG96.cc\
                   OutG96S.cc\
                   OutPdb.cc\
                   Outvmdam.cc\
                   InParameter.cc\
                   InBuildingBlock.cc

check_PROGRAMS = InTopology\
               InG96\
               OutG96S\
               OutPdb\
               Outvmdam\
               InParameter\
               InBuildingBlock

LDADD = libgio.la ../gcore/libgcore.la

InTopology_SOURCES = InTopology.t.cc
InG96_SOURCES = InG96.t.cc
OutG96S_SOURCES = OutG96S.t.cc
OutPdb_SOURCES = OutPdb.t.cc
Outvmdam_SOURCES = Outvmdam.t.cc
InParameter_SOURCES = InParameter.t.cc
InBuildingBlock_SOURCES = InBuildingBlock.t.cc
