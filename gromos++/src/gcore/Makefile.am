## Template for automake for libgromos

## will not be built during development
lib_LTLIBRARIES = libgcore.la

libgcore_la_LDFLAGS = -version-info $(GROMOS_MAJOR_VERSION):$(GROMOS_MICRO_VERSION):0


include_HEADERS = Angle.h\
                     AngleType.h\
                     AtomPair.h\
                     AtomTopology.h\
                     MassType.h\
                     Bond.h\
                     BondType.h\
                     Box.h\
                     Constraint.h\
                     Dihedral.h\
                     DihedralType.h\
                     Exclusion.h\
                     GromosForceField.h\
                     Improper.h\
                     ImproperType.h\
                     LJType.h\
                     MoleculeTopology.h\
                     Molecule.h\
                     SolventTopology.h\
                     Solvent.h\
                     BuildingBlock.h\
                     BbSolute.h\
                     System.h

libgcore_la_SOURCES = Angle.cc\
                     AngleType.cc\
                     AtomPair.cc\
                     AtomTopology.cc\
                     MassType.cc\
                     Bond.cc\
                     BondType.cc\
                     Box.cc\
                     Constraint.cc\
                     Dihedral.cc\
                     DihedralType.cc\
                     Exclusion.cc\
                     GromosForceField.cc\
                     Improper.cc\
                     ImproperType.cc\
                     LJType.cc\
                     MoleculeTopology.cc\
                     Molecule.cc\
                     SolventTopology.cc\
                     Solvent.cc\
                     BuildingBlock.cc\
                     BbSolute.cc\
                     System.cc

check_PROGRAMS = AtomTopology\
               Exclusion\
               MoleculeTopology\
               Molecule\
               SolventTopology\
               Solvent\
               BbSolute

LDADD = libgcore.la

AtomTopology_SOURCES = AtomTopology.t.cc
Exclusion_SOURCES = Exclusion.t.cc
MoleculeTopology_SOURCES = MoleculeTopology.t.cc
Molecule_SOURCES = Molecule.t.cc
SolventTopology_SOURCES = SolventTopology.t.cc
Solvent_SOURCES = Solvent.t.cc
BbSolute_SOURCES = BbSolute.t.cc
