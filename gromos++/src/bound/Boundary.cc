// bound_Boundary.cc

#include "Boundary.h"
#include "../gmath/Vec.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"

#ifndef INCLUDED_GCORE_BOX
#include "../gcore/Box.h"
#define INCLUDED_GCORE_BOX
#endif

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
  for(int i=0; i<d_this->d_sys->numMolecules();++i){
    //the system might not have coordinates yet!
    if(d_this->d_sys->mol(i).numPos()==0)
      d_this->d_sys->mol(i).initPos();
    d_this->d_ref.push_back(&d_this->d_sys->mol(i).pos(0));
  }
}

void Boundary::setReference(int i, const Vec &vec){
  assert (i<int(d_this->d_ref.size()));
  d_this->d_ref[i]=&vec;
}

void Boundary::setReference(System const & sys)
{
  assert(int(d_this->d_ref.size()) == sys.numMolecules());
  
  for(int i=0; i < sys.numMolecules(); ++i){

    d_this->d_ref[i] = &sys.mol(i).pos(0);
  }
}

bound::Boundary::~Boundary(){
  delete d_this;
}

const gmath::Vec &Boundary::reference(int i)const{
  assert(d_this != NULL);
  
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

//bool Boundary::isInBox(const gmath::Vec &r, const gcore::Box &box) const {

//  gmath::Vec boxh(box[0], box[1], box[2]);
//  return (r == nearestImage(boxh, r, box));
//}
