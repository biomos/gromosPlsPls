// gcore_Solvent.cc

#include "Solvent.h"
#include "SolventTopology.h"
#include "../gmath/Vec.h"
#include <vector>

using namespace std;
using gcore::Solvent;
using gcore::SolventTopology;
using gmath::Vec;

Solvent::Solvent(const SolventTopology &mt):
  d_mt(new SolventTopology(mt)),
  d_pos(),
  d_numCoords()
{
  d_numCoords = 0;
}

Solvent::Solvent(const Solvent &solv):
  d_mt(new SolventTopology(*solv.d_mt)),
  d_pos(solv.d_pos.size()),
  d_numCoords()
{
  for(int i=0; i<solv.numCoords();++i)
    d_pos[i]=new Vec(solv.pos(i));
  d_numCoords=solv.d_numCoords;
}

Solvent::~Solvent(){
  delete d_mt;
  for(int i=0; i<int(d_pos.size());++i)
    delete d_pos[i];
}

void Solvent::addCoord(Vec v){

    d_pos.push_back(new Vec(v));
  d_numCoords++;
   return;
  }

void Solvent::setnumCoords(int i){
  d_numCoords = i;
  for(int j=d_pos.size();j > i;--j){
    delete d_pos[j-1];
  }
  d_pos.resize(i);
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


int Solvent::numCoords()const{
  return d_numCoords;
}

