// gcore_Molecule.cc
/**
 * Class Molecule
 * Addition: velocity configuration added to Molecule definition;
 * Author  : gee          
 */

#include <cassert>
#include <vector>
#include <set>
#include "AtomPair.h"
#include "LJException.h"
#include "MoleculeTopology.h"
#include "../gmath/Vec.h"
#include "Molecule.h"

using namespace std;
using gcore::Molecule;
using gcore::MoleculeTopology;
using gmath::Vec;

Molecule::Molecule(const MoleculeTopology &mt):
  d_mt(new MoleculeTopology(mt)),
  d_pos(0),
  d_vel(0),
  d_bfac(0),
  d_cosDisplacement(0){}

Molecule::Molecule(const Molecule &mol):
  d_mt(new MoleculeTopology(*mol.d_mt)),
  d_pos(mol.d_pos.size()),
  d_vel(mol.d_vel.size()),
  d_bfac(mol.d_bfac.size()),
  d_cosDisplacement(mol.d_cosDisplacement.size())
{
  for(int i=0; i<mol.numPos();++i)
  {
    d_pos[i]=new Vec(mol.pos(i));
  }
  for(int i=0; i<mol.numVel();++i){
    d_vel[i]=new Vec(mol.vel(i));
  }
  for(int i=0; i<mol.numBfac();++i){
    d_bfac[i]=0.0;
  }
  for(int i=0; i<mol.numCosDisplacements();++i){
    d_cosDisplacement[i] = new Vec(mol.cosDisplacement(i));
  }
}

Molecule::~Molecule(){
  delete d_mt;
  for(unsigned int i=0; i<d_pos.size();++i)
    delete d_pos[i];

  for(unsigned int i=0; i<d_vel.size();++i)
    delete d_vel[i];

  for(unsigned int i=0; i<d_cosDisplacement.size();++i)
    delete d_cosDisplacement[i];
}

void Molecule::initPos(){
  d_pos.resize(numAtoms());
  for(int i=0; i < numAtoms(); ++i){
    d_pos[i]=new Vec();
  }
}
void Molecule::initVel(){
  d_vel.resize(numAtoms());
  for(int i=0; i < numAtoms(); ++i){
    d_vel[i]=new Vec();
  }
}
void Molecule::initBfac(){
  d_bfac.resize(numAtoms());
  for(int i=0; i < numAtoms(); ++i){
    d_bfac[i]=0.0;
  }
}
void Molecule::setBfac(int i, double b){
    d_bfac[i]=b;
}
void Molecule::initCosDisplacements(){
  d_cosDisplacement.resize(numAtoms());
  for(int i=0; i < numAtoms(); ++i){
    d_cosDisplacement[i]=new Vec();
  }
}
MoleculeTopology &Molecule::topology()
{
  return *d_mt;
}

const MoleculeTopology &Molecule::topology()const{
  return *d_mt;
}

int Molecule::numAtoms()const{
  return d_mt->numAtoms();
}

