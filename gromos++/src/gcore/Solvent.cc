// gcore_Solvent.cc
/**
 * Class Solvent
 * Addition: velocity configuration added to Molecule definition;
 * Author  : gee          
 */

#include <cassert>
#include "Solvent.h"
#include "SolventTopology.h"
#include "../gmath/Vec.h"
#include <vector>

using namespace std;
using gcore::Solvent;
using gcore::SolventTopology;
using gmath::Vec;

Solvent::Solvent(const SolventTopology &mt) :
d_mt(new SolventTopology(mt)),
d_pos(),
d_numPos(0),
d_vel(),
d_numVel(0),
d_cosDisplacement(),
d_numCosDisplacements(0) {
}

Solvent::Solvent(const Solvent &solv):
  d_mt(new SolventTopology(*solv.d_mt)),
  d_pos(solv.d_pos.size()),
  d_vel(solv.d_vel.size()),
  d_cosDisplacement(solv.d_cosDisplacement.size())
{
  for(int i=0; i<solv.numPos();++i){
    d_pos[i]=new Vec(solv.pos(i));
  }
  d_numPos=solv.d_numPos;
  for(int i=0; i<solv.numVel();++i){
    d_vel[i]=new Vec(solv.vel(i));
  }
  d_numVel=solv.d_numVel;
  for(int i=0; i<solv.numCosDisplacements();++i){
    d_cosDisplacement[i]=new Vec(solv.cosDisplacement(i));
  }
  d_numCosDisplacements=solv.d_numCosDisplacements;
}

Solvent::~Solvent(){
  delete d_mt;
  for(unsigned int i=0; i<d_pos.size();++i)
    delete d_pos[i];
  
  for(unsigned int i=0; i<d_vel.size();++i)
    delete d_vel[i];
  
  for(unsigned int i=0; i<d_cosDisplacement.size();++i)
    delete d_cosDisplacement[i];
}

void Solvent::addPos(Vec v){
  d_pos.push_back(new Vec(v));
  d_numPos++;
  return;
}

void Solvent::addVel(Vec v){
  d_vel.push_back(new Vec(v));
  d_numVel++;
  return;
}

void Solvent::addCosDisplacement(Vec v){
  d_cosDisplacement.push_back(new Vec(v));
  d_numCosDisplacements++;
  return;
}

void Solvent::setNumPos(int i){
  d_numPos = i;
  for(int j=d_pos.size();j > i;--j)
    delete d_pos[j-1];

  d_pos.resize(i);
}

void Solvent::setNumVel(int i){
  d_numVel = i;
  for(int j=d_vel.size();j > i;--j)
    delete d_vel[j-1];

  d_vel.resize(i);
}

void Solvent::setNumCosDisplacements(int i){
  d_numCosDisplacements = i;
  for(int j=d_cosDisplacement.size();j > i;--j)
    delete d_cosDisplacement[j-1];

  d_cosDisplacement.resize(i);
}


Solvent &Solvent::operator=(const Solvent &s){
  if(this != &s){
    this->~Solvent();
    new(this) Solvent(s);
  }
  return *this;
}

const SolventTopology &Solvent::topology()const{
  return *d_mt;
}
SolventTopology &Solvent::topology(){
  return *d_mt;
}


int Solvent::numPos()const{
  return d_numPos;
}

int Solvent::numVel()const{
  return d_numVel;
}

int Solvent::numCosDisplacements()const{
  return d_numCosDisplacements;
}

int Solvent::numAtoms()const
{
  return std::max(std::max(d_numPos, d_numVel), d_numCosDisplacements);
}


