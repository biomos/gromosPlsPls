// gcore_System.cc

#include <cassert>
#include "System.h"
#include "Molecule.h"
#include "MoleculeTopology.h"
#include "Solvent.h"
#include "Box.h"
#include <new>

using gcore::System;

System::System():
  d_mol(),
  d_sol()
{
  d_box=new Box();
  hasPos = false;
  hasBox = false;
  hasVel = false;
}
 

System::System(const System &sys):
  d_mol(sys.d_mol.size()),
  d_sol(sys.d_sol.size())
{
  for (unsigned int i=0; i<d_mol.size();++i){
    d_mol[i]= new Molecule(sys.mol(i));
  }
  for (unsigned int i=0; i<d_sol.size();++i){
    d_sol[i]=new Solvent(sys.sol(i));
  }
  d_box = new Box(sys.box());
  hasBox = sys.hasBox;
  hasPos = sys.hasPos;
  hasVel = sys.hasVel;
}

System::~System(){
  for (unsigned int i=0; i<d_mol.size();++i){
    delete d_mol[i];
  }
  for (unsigned int i=0; i<d_sol.size();++i){
    delete d_sol[i];
  }
  delete d_box;
}

System &System::operator=(const System &sys){
  if(this != &sys){
    // delete this;
    this->~System();
    new(this) System(sys);
  }
  return *this;
}

void System::addMolecule(const Molecule &mol){
  d_mol.push_back(new Molecule(mol));
}

void System::addSolvent(const Solvent &sol){
  d_sol.push_back(new Solvent(sol));
}




