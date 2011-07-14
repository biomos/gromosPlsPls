// gcore_System.cc

#include <cassert>
#include <set>
#include <new>
#include "Molecule.h"
#include "LJException.h"
#include "MoleculeTopology.h"
#include "Solvent.h"
#include "Box.h"
#include "Remd.h"
#include "../gmath/Vec.h"
#include "System.h"
#include "Weight.h"

using gcore::System;
using gmath::Vec;

System::System():
  d_mol(),
  d_sol(),
  d_temperatureGroup(),
  d_pressureGroup(),
  primlist()
{
  d_box=new Box();
  d_remd=new Remd();
  hasPos = false;
  hasBox = false;
  hasVel = false;
  hasCosDisplacements = false;
  hasRemd = false;
  d_weight = new Weight();
}
 

System::System(const System &sys):
  d_mol(sys.d_mol.size()),
  d_sol(sys.d_sol.size()),
  d_temperatureGroup(sys.d_temperatureGroup.size()),
  d_pressureGroup(sys.d_pressureGroup.size()),
  primlist(sys.d_mol.size())
{
  for (unsigned int i=0; i<d_mol.size();++i){
    d_mol[i]= new Molecule(sys.mol(i));
  }
  for (unsigned int i=0; i<d_sol.size();++i){
    d_sol[i]=new Solvent(sys.sol(i));
  }
  for (unsigned int i=0; i<d_mol.size();++i){
      primlist[i][0] = 0;
      primlist[i][1] = i-1;
      primlist[i][2] = 0;
  }
  for (unsigned int i=0; i<d_temperatureGroup.size();++i){
    d_temperatureGroup[i]= new int(sys.temperatureGroup(i));
  }
  for (unsigned int i=0; i<d_pressureGroup.size();++i){
    d_pressureGroup[i]=new int(sys.pressureGroup(i));
  }
  d_box = new Box(sys.box());
  d_remd = new Remd(sys.remd());
  hasBox = sys.hasBox;
  hasPos = sys.hasPos;
  hasVel = sys.hasVel;
  hasCosDisplacements = sys.hasCosDisplacements;
  hasRemd = sys.hasRemd;
  d_weight = new Weight(sys.weight());
}

System::~System(){
  for (unsigned int i=0; i<d_mol.size();++i){
    delete d_mol[i];
  }
  for (unsigned int i=0; i<d_sol.size();++i){
    delete d_sol[i];
  }
  for (unsigned int i=0; i<d_temperatureGroup.size();++i){
    delete d_temperatureGroup[i];
  }
  for (unsigned int i=0; i<d_pressureGroup.size();++i){
    delete d_pressureGroup[i];
  }
  delete d_box;
  delete d_remd;
  delete d_weight;
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

void System::addTemperatureGroup(const int &tg){
  d_temperatureGroup.push_back(new int(tg));
}

void System::addPressureGroup(const int &pg){
  d_pressureGroup.push_back(new int(pg));
}


