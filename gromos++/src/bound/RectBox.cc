// bound_RectBox.cc

#include "RectBox.h"
#include "../gmath/Vec.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/Box.h"
#include <cmath>
#include <iostream>

using namespace std;
using bound::RectBox;
using gmath::Vec;
using namespace gcore;

static Vec nim(const Vec &r1,const  Vec &r2, const Box &box){

  Vec diff=r2-r1;
  Vec a;
  
  a[0] = diff[0] - box[0] * rint(diff[0]/box[0]);
  a[1] = diff[1] - box[1] * rint(diff[1]/box[1]);
  a[2] = diff[2] - box[2] * rint(diff[2]/box[2]);


  Vec rec = r1 + a;
  
  return rec;
  
}


Vec RectBox::nearestImage(const Vec &v1, const Vec &v2, const Box &box)const{
  return nim(v1, v2, box);
}

void RectBox::gather(){

  if (!sys().hasBox) throw gromos::Exception("Gather problem",  
                              "System does not contain Box block! Abort!");

  if (sys().box()[0] == 0 || sys().box()[1] == 0 || sys().box()[2] == 0) 
    throw gromos::Exception("Gather problem",  
			    "Box block contains element(s) of value 0.0! Abort!");  

  for(int i=0; i<sys().numMolecules();++i){
    Molecule &mol=sys().mol(i);
    mol.pos(0)=nim(reference(0),mol.pos(0),sys().box());
    for(int j=1;j<mol.numPos();++j)
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

void RectBox::gathergr(){

  if (!sys().hasBox) 
    throw gromos::Exception("Gather problem",  
			    "System does not contain Box block! Abort!");

   if (sys().box()[0] == 0 || sys().box()[1] == 0 || sys().box()[2] == 0)
     throw gromos::Exception("Gather problem",  
			     "Box block contains element(s) of value 0.0! Abort!");
   
   for(int i=0; i<sys().numMolecules();++i){
     Molecule &mol=sys().mol(i);
     mol.pos(0)=nim(reference(i),mol.pos(0),sys().box());
     for(int j=1;j<mol.numAtoms();++j)
       mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());
   }
}

void RectBox::coggather(){

  if (!sys().hasBox)
    throw gromos::Exception("Gather problem",  
			    "System does not contain Box block! Abort!");
  
  if (sys().box()[0] == 0 || sys().box()[1] == 0 || sys().box()[2] == 0)
    throw gromos::Exception("Gather problem",  
			    "Box block contains element(s) of value 0.0! Abort!");
  
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
    for(int j=1;j<mol.numPos();++j){
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

