// gcore_Molecule.cc

#include "Molecule.h"
#include "MoleculeTopology.h"
#include "../gmath/Vec.h"

using namespace std;
using gcore::Molecule;
using gcore::MoleculeTopology;
using gmath::Vec;

Molecule::Molecule(const MoleculeTopology &mt):
  d_mt(new MoleculeTopology(mt)),
  d_pos(mt.numAtoms())
{
  for(int i=0; i<mt.numAtoms();++i)
    d_pos[i]=new Vec();
}

Molecule::Molecule(const Molecule &mol):
  d_mt(new MoleculeTopology(*mol.d_mt)),
  d_pos(mol.d_pos.size())
{
  for(int i=0; i<mol.numAtoms();++i)
    d_pos[i]=new Vec(mol.pos(i));
}

Molecule::~Molecule(){
  delete d_mt;
  for(int i=0; i<int(d_pos.size());++i)
    delete d_pos[i];
}
MoleculeTopology &Molecule::topology()
{
  return *d_mt;
}

const MoleculeTopology &Molecule::topology()const{
  return *d_mt;
}

int Molecule::numAtoms()const{
  return d_pos.size();
}

