## Template for automake for libgromos

## will not be built during development
lib_LTLIBRARIES = libgcore.la

libgcore_la_LDFLAGS = -version-info $(GROMOS_MAJOR_VERSION):$(GROMOS_MICRO_VERSION):0

## List to be completed...

# noinst_LIBRARIES = libgcore.la

libgcore_la_SOURCES = Angle.cc\
                     AngleType.cc\
                     AtomPair.cc\
                     AtomTopology.cc\
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
                     System.cc

bin_PROGRAMS = AtomTopology\
               Exclusion\
               MoleculeTopology\
               Molecule\
               SolventTopology\
               Solvent

LDADD = libgcore.la

AtomTopology_SOURCES = AtomTopology.t.cc
Exclusion_SOURCES = Exclusion.t.cc
MoleculeTopology_SOURCES = MoleculeTopology.t.cc
Molecule_SOURCES = Molecule.t.cc
SolventTopology_SOURCES = SolventTopology.t.cc
Solvent_SOURCES = Solvent.t.cc
