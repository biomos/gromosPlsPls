## Template for automake for libgromos

lib_LTLIBRARIES = libgio.la

libgio_la_LDFLAGS = -version-info $(GROMOS_MAJOR_VERSION):$(GROMOS_MICRO_VERSION):0

include_HEADERS = Ginstream.h\
                   BlockInput.h\
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
                    InG96.cc\
                   InTopology.cc\
                    InParameter.cc\
                    InBuildingBlock.cc\
                   OutTopology.cc\
                   OutCoordinates.cc\
                   OutG96.cc\
                   OutG96S.cc\
                   OutPdb.cc\
                   Outvmdam.cc

check_PROGRAMS = Ginstream\
                 InG96\
	       InTopology\
               InParameter\
               InBuildingBlock\
               OutG96S\
               OutPdb\
               Outvmdam

LDADD = libgio.la ../gcore/libgcore.la

InG96_SOURCES = InG96.t.cc
InTopology_SOURCES = InTopology.t.cc
InParameter_SOURCES = InParameter.t.cc
InBuildingBlock_SOURCES = InBuildingBlock.t.cc
Ginstream_SOURCES = Ginstream.t.cc
OutG96S_SOURCES = OutG96S.t.cc
OutPdb_SOURCES = OutPdb.t.cc
Outvmdam_SOURCES = Outvmdam.t.cc
