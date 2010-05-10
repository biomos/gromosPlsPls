// bound_Boundary.cc

#include "Boundary.h"
#include "../gmath/Vec.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gio/InG96.h"

#ifndef INCLUDED_GCORE_BOX
#include "../gcore/Box.h"
#define INCLUDED_GCORE_BOX
#endif

#include <vector>
#include <cassert>
#include <iostream>
#include <sstream>

using gmath::Vec;
using gcore::System;
using bound::Boundary;
using bound::Boundary_i;
using namespace std;

class Boundary_i{
  friend class bound::Boundary;
  gcore::System *d_sys;
  gcore::System *d_refSys;
  vector<const gmath::Vec *> d_ref;
  char d_type;
  Boundary_i(): d_ref(){}
  ~Boundary_i(){}
};

Boundary::Boundary(System *sys): d_this(new Boundary_i()) {
  d_this->d_sys=sys;
  d_this->d_refSys = NULL;
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

void Boundary::setReferenceFrame(std::string file) {
  gio::InG96 in(file);
  in.select("ALL");
  d_this->d_refSys = new System(sys());
  in >> refSys();
}

bound::Boundary::~Boundary(){
  if (d_this->d_refSys != NULL)
    delete d_this->d_refSys;
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

void Boundary::setType(char t) {
  switch (t) {
    case 'v':
      d_this->d_sys->box().setNtb(gcore::Box::vacuum);
      break;
    case 'r':
      d_this->d_sys->box().setNtb(gcore::Box::rectangular);
      break;
    case 't':
      d_this->d_sys->box().setNtb(gcore::Box::truncoct);
      break;
    case 'c':
      d_this->d_sys->box().setNtb(gcore::Box::triclinic);
      break;
    default:
      stringstream msg;
      msg << "periodic boundary condition '" << t << "' unknow. Known "
              "boundaries are r (rectangular), t (truncated octahedron), c (triclinic) and v (vacuum)";
      throw gromos::Exception("Boundary", msg.str());
      break;
  }
  d_this->d_type = t;

}

System &Boundary::sys(){
  return *d_this->d_sys;
}

System &Boundary::refSys(){
  return *d_this->d_refSys;
}

//bool Boundary::isInBox(const gmath::Vec &r, const gcore::Box &box) const {

//  gmath::Vec boxh(box[0], box[1], box[2]);
//  return (r == nearestImage(boxh, r, box));
//}
