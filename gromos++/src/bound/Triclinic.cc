// bound_Triclinic.cc

#include "Triclinic.h"
#include "../gmath/Vec.h"
#include "../gcore/System.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Bond.h"
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

void Triclinic::nogather(){

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

void Triclinic::crsgather(){

  if (!sys().hasBox)
    throw gromos::Exception("Gather problem",
               "System does not contain Box block! Abort!");

  if (sys().box()[0] == 0 || sys().box()[1] == 0 || sys().box()[2] == 0)
    throw gromos::Exception("Gather problem",
                "Box block contains element(s) of value 0.0! Abort!");

  // Reconstruct the connectivity of the submolecules
  for(int i=0; i<sys().numMolecules();++i){
    Molecule &mol=sys().mol(i);
    mol.pos(0)=nim(reference(i),mol.pos(0),sys().box());
    for(int j=1;j<mol.numAtoms();++j)
      mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());
  }

  // Determine the positions of the centres of geometry of the gathered molecules 
  // and store them in vcog
  std::vector<Vec> vcog;
  for(int i=0; i<sys().numMolecules();++i){
    Vec cog(0.0,0.0,0.0);
    int numat=0;
    for (int j=0; j<sys().mol(i).numAtoms(); ++j) {
      cog = cog + sys().mol(i).pos(j);
      ++numat;
    }
    cog = (1.0/double(numat))*cog;
    vcog.push_back(cog);
  }

  // Gather nearest image of cog of molecule 1 w.r.t. origin
  // vcog[0]=nim((0.0,0.0,0.0),vcog[0],sys().box());

  // Now gather cog's w.r.t. cog of previous molecule
  // ocog: buffer to store the overall cog of already gathered cog's
  Vec ocog=vcog[0];
  for(int i=0; i<sys().numMolecules()-1;++i){
    Vec nimj=nim(vcog[i],vcog[i+1],sys().box());
    vcog[i+1]=nimj;
    ocog+=vcog[i+1];
  }

  // Now gather the atoms of the solute molecules with
  // the newly determined cog's of the respective molecule
  for(int i=0;i<sys().numMolecules();++i){
    Molecule &mol=sys().mol(i);
    for(int j=0;j<mol.numAtoms();++j){
      mol.pos(j)=nim(vcog[i],mol.pos(j),sys().box());
    }
  }

  // Gather the solvent molecules with ocog as a reference
  ocog=ocog/double(sys().numMolecules());
  Solvent &sol=sys().sol(0);
  for(int i=0;i<sol.numPos();i+=sol.topology().numAtoms()){
    sol.pos(i)=nim(ocog,sol.pos(i),sys().box());
    for(int j=i+1;j<(i+sol.topology().numAtoms());++j){
      sol.pos(j)=nim(sol.pos(j-1),sol.pos(j),sys().box());
    }
  }
}

void Triclinic::seqgather(){

  if (!sys().hasBox)
    throw gromos::Exception("Gather problem",
               "System does not contain Box block! Abort!");

  if (sys().box()[0] == 0 || sys().box()[1] == 0 || sys().box()[2] == 0)
    throw gromos::Exception("Gather problem",
                "Box block contains element(s) of value 0.0! Abort!");

  // Reconstruct the connectivity of the submolecules
  for(int i=0; i<sys().numMolecules();++i){
    Molecule &mol=sys().mol(i);
    mol.pos(0)=nim(reference(i),mol.pos(0),sys().box());
    for(int j=1;j<mol.numAtoms();++j)
      mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());
  }

  // Determine the positions of the centres of geometry of the gathered molecules 
  // and store them in vcog
  std::vector<Vec> vcog;
  for(int i=0; i<sys().numMolecules();++i){
    Vec cog(0.0,0.0,0.0);
    int numat=0;
    for (int j=0; j<sys().mol(i).numAtoms(); ++j) {
      cog = cog + sys().mol(i).pos(j);
      ++numat;
    }
    cog = (1.0/double(numat))*cog;
    vcog.push_back(cog);
  }

  // Gather nearest image of cog of molecule 1 w.r.t. origin
  // vcog[0]=nim((0.0,0.0,0.0),vcog[0],sys().box());

  // Now gather cog's w.r.t. cog of previous molecule
  // ocog: buffer to store the overall cog of already gathered cog's
  Vec ocog=vcog[0];
  for(int i=0; i<sys().numMolecules()-1;++i){
    // crs:
    // Vec nimj=nim(vcog[i],vcog[i+1],sys().box());
    // seq:
    Vec nimj=nim(ocog/double(i+1),vcog[i+1],sys().box());
    vcog[i+1]=nimj;
    ocog+=vcog[i+1];
  }

  // Now gather the atoms of the solute molecules with
  // the newly determined cog's of the respective molecule
  for(int i=0;i<sys().numMolecules();++i){
    Molecule &mol=sys().mol(i);
    for(int j=0;j<mol.numAtoms();++j){
      mol.pos(j)=nim(vcog[i],mol.pos(j),sys().box());
    }
  }

  // Gather the solvent molecules with ocog as a reference
  ocog=ocog/double(sys().numMolecules());
  Solvent &sol=sys().sol(0);
  for(int i=0;i<sol.numPos();i+=sol.topology().numAtoms()){
    sol.pos(i)=nim(ocog,sol.pos(i),sys().box());
    for(int j=i+1;j<(i+sol.topology().numAtoms());++j){
      sol.pos(j)=nim(sol.pos(j-1),sol.pos(j),sys().box());
    }
  }
}

void Triclinic::gengather(){

  if (!sys().hasBox)
    throw gromos::Exception("Gather problem",
               "System does not contain Box block! Abort!");

  if (sys().box()[0] == 0 || sys().box()[1] == 0 || sys().box()[2] == 0)
    throw gromos::Exception("Gather problem",
                "Box block contains element(s) of value 0.0! Abort!");

  // Reconstruct the connectivity of the submolecules
  for(int i=0; i<sys().numMolecules();++i){
    Molecule &mol=sys().mol(i);
    mol.pos(0)=nim(reference(i),mol.pos(0),sys().box());
    for(int j=1;j<mol.numAtoms();++j)
      mol.pos(j)=nim(mol.pos(j-1),mol.pos(j),sys().box());
  }

  // Determine the positions of the centres of geometry of the gathered molecules 
  // and store them in vcog
  std::vector<Vec> vcog;
  for(int i=0; i<sys().numMolecules();++i){
    Vec cog(0.0,0.0,0.0);
    int numat=0;
    for (int j=0; j<sys().mol(i).numAtoms(); ++j) {
      cog = cog + sys().mol(i).pos(j);
      ++numat;
    }
    cog = (1.0/double(numat))*cog;
    vcog.push_back(cog);
  }

  // Gather nearest image of cog of molecule 1 w.r.t. origin
  // vcog[0]=nim((0.0,0.0,0.0),vcog[0],sys().box());

  // Use vcogi to make the graph connecting the closest cog's
  std::vector<int> vcogi;
  for(int i=0; i<sys().numMolecules();++i){
    vcogi.push_back(i);
  }

  // Now gather cog's w.r.t. each other
  // ocog: buffer to store the overall cog of already gathered cog's
  Vec ocog=vcog[0];
  for(int i=0; i<sys().numMolecules()-1;++i){
    // Determine closest nim to i among remaining molecules (using vcogi)
    int bufi=vcogi[i];
    int inimcogi=vcogi[i+1];
    Vec nimcogi=nim(vcog[bufi],vcog[inimcogi],sys().box());
    int jclose=i+1;
    for(int j=i+2; j<sys().numMolecules();++j){
      int bufj=vcogi[j];
      if( (nim(vcog[bufi],vcog[bufj],sys().box())-vcog[bufi]).abs()<(nimcogi-vcog[bufi]).abs()){
        nimcogi=nim(vcog[bufi],vcog[bufj],sys().box());
        inimcogi=bufj;
        jclose=j;
      }
    }
    // Now swap inimcogi with i+1 in vcogi
    int bufci=vcogi[i+1];
    vcogi[i+1]=inimcogi;
    vcogi[jclose]=bufci;
    
    // Set vcog[i+1] either to its nim to vcog[i], or to
    // nim to overall cog of molecules[1 ... i], depending
    // on what corresponds with the closest distance
    Vec nic1=nimcogi;
    Vec nic2=nim(ocog/double(i+1),nimcogi,sys().box());
    if((nic1-vcog[bufi]).abs()<(nic2-ocog/double(i+1)).abs()){
      vcog[inimcogi]=nic1;
    }
    else{
      vcog[inimcogi]=nic2;
    }
    ocog+=vcog[inimcogi];
  }

  // Now gather the atoms of the solute molecules with
  // the newly determined cog's of the respective molecule
  // as a reference
  for(int i=0;i<sys().numMolecules();++i){
    Molecule &mol=sys().mol(i);
    for(int j=0;j<mol.numAtoms();++j){
      mol.pos(j)=nim(vcog[i],mol.pos(j),sys().box());
    }
  }

  // Gather the solvent molecules with ocog as a reference
  ocog=ocog/double(sys().numMolecules());
  Solvent &sol=sys().sol(0);
  for(int i=0;i<sol.numPos();i+=sol.topology().numAtoms()){
    sol.pos(i)=nim(ocog,sol.pos(i),sys().box());
    for(int j=i+1;j<(i+sol.topology().numAtoms());++j){
      sol.pos(j)=nim(sol.pos(j-1),sol.pos(j),sys().box());
    }
  }
}

void Triclinic::bondgather(){

  if (!sys().hasBox) throw gromos::Exception("Gather problem",  
                              "System does not contain Box block! Abort!");

  if (sys().box()[0] == 0 || sys().box()[1] == 0 || sys().box()[2] == 0) 
    throw gromos::Exception("Gather problem",  
			    "Box block contains element(s) of value 0.0! Abort!");  

  for(int i=0; i<sys().numMolecules();++i){
    Molecule &mol=sys().mol(i);
    mol.pos(0)=nim(reference(0),mol.pos(0),sys().box());
    for(int j=1;j<mol.numPos();++j){
      
      //find a previous atom to which we are bound
      BondIterator bi(mol.topology());
      int k=0;
      for(;bi;++bi)
	if(bi()[1]==j) { k = bi()[0]; break; }
      mol.pos(j)=nim(mol.pos(k),mol.pos(j),sys().box());
    }
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

