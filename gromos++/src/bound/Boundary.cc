// bound_Boundary.cc

#include "Boundary.h"
#include "../gmath/Vec.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"

#include <vector>
#include <cassert>
#include <iostream>

using gmath::Vec;
using gcore::System;
using bound::Boundary;
using bound::Boundary_i;
using namespace std;

class Boundary_i{
  friend class Boundary;
  gcore::System *d_sys;
  vector<const gmath::Vec *> d_ref;
  char d_type;
  Boundary_i(): d_ref(){}
  ~Boundary_i(){}
};

Boundary::Boundary(System *sys): d_this(new Boundary_i()){
  d_this->d_sys=sys;
  for(int i=0; i<d_this->d_sys->numMolecules();++i)
    d_this->d_ref.push_back(&d_this->d_sys->mol(i).pos(0));
  // cout << d_this->d_sys->numSolvents() ;
  //  for(int i=0; i<d_this->d_sys->numSolvents();++i){
  // d_this->d_ref.push_back(&d_this->d_sys->sol(0).pos(0));//}
}

void Boundary::setReference(int i, const Vec &vec){
  assert (i<int(d_this->d_ref.size()));
  d_this->d_ref[i]=&vec;
}


bound::Boundary::~Boundary(){
  delete d_this;
}

const gmath::Vec &Boundary::reference(int i)const{
  assert(i<int(d_this->d_ref.size()));
  return *d_this->d_ref[i];
}

char Boundary::type()
{
  return d_this->d_type;
}

void Boundary::setType(char t)
{
  d_this->d_type=t;
}

System &Boundary::sys(){
  return *d_this->d_sys;
}
