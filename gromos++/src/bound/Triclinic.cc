// bound_Triclinic.cc

#include "Triclinic.h"
#include "../gmath/Vec.h"
#include "../gcore/System.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/Molecule.h"
#include "../gcore/Box.h"
#include <cmath>
//#include <iostream>;

using bound::Triclinic;
using gmath::Vec;
using namespace gcore;

static Vec nim(const Vec &r1,const  Vec &r2, const Box &box){

  Vec P = r2 - r1;
  int k,l,m;
  k = int(rint(box.cross_K_L_M()[0].dot(P)));
  l = int(rint(box.cross_K_L_M()[1].dot(P)));
  m = int(rint(box.cross_K_L_M()[2].dot(P)));
  
  P += box.K() * k + box.L() * l + box.M() * m;
  
  return r1 + P;
}

Vec Triclinic::nearestImage(const Vec &v1, const Vec &v2, const Box &box)const{
  return nim(v1, v2, box);
}

void Triclinic::gathergr(){
    for(int i=0; i<sys().numMolecules();++i){
    Molecule &mol=sys().mol(i);
    mol.pos(0)=nim(reference(i),mol.pos(0),sys().box());
    for(int j=1;j<mol.numAtoms();++j)
      mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());
  }
}


void Triclinic::gather(){
  for(int i=0; i<sys().numMolecules();++i){
    Molecule &mol=sys().mol(i);
    mol.pos(0)=nim(reference(i),mol.pos(0),sys().box());
    for(int j=1;j<mol.numAtoms();++j)
      mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());
  }
 // do the solvent 
   Solvent &sol=sys().sol(0);
   for(int i=0;i<sol.numPos();i+= sol.topology().numAtoms()){
   sol.pos(i)=nim(reference(0),sol.pos(i),sys().box());
   for (int j=i+1;j < (i+sol.topology().numAtoms());++j){
    sol.pos(j)=nim(sol.pos(j-1),sol.pos(j),sys().box());
   }
  }
}

void Triclinic::coggather(){
   Molecule &mol=sys().mol(0);
   Solvent &sol=sys().sol(0);
  
   Vec ref(0.0,0.0,0.0);
   Vec cog;
   int atoms=0;

   // do mol(0) with respect to ref (0,0,0)
   mol.pos(0)=nim(ref,mol.pos(0),sys().box());
   for(int j=1;j<mol.numAtoms();++j){
     mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());}

   // calculate COG of mol(0)
    for (int i=0;i < mol.numAtoms(); i++) {
      cog = cog + mol.pos(i);
      ++atoms;
    }
   cog = (1.0/double(atoms))*cog;

   // do the rest of the molecules
 for(int i=1;i<sys().numMolecules();++i){
      Molecule &mol=sys().mol(i);
      mol.pos(0)=nim(cog,mol.pos(0),sys().box());      
   for(int j=1;j<mol.numAtoms();++j){
    mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());
   }
 }

 // do the solvent 
   for(int i=0;i<sol.numPos();i+= sol.topology().numAtoms()){
   sol.pos(i)=nim(cog,sol.pos(i),sys().box());
   for (int j=i+1;j < (i+sol.topology().numAtoms());++j){
    sol.pos(j)=nim(sol.pos(j-1),sol.pos(j),sys().box());
   }
  }
} 

