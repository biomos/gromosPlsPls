#include <iostream>
#include <string>
#include <cassert>
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gcore/Exclusion.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/Box.h"
#include "../bound/Boundary.h"
#include "../gmath/Vec.h"
#include "AtomSpecifier.h"
#include "SimplePairlist.h"

namespace utils{

  SimplePairlist::SimplePairlist(gcore::System &sys, 
				 bound::Boundary &pbc, 
				 double c){
    setSystem(sys);
    d_pbc = &pbc;
    d_cut2 = c*c;
    d_atom = -1;
    d_chargeGroupBased=false;
  }
  
  void SimplePairlist::setPbc(bound::Boundary &pbc)
  {
    d_pbc = &pbc;
  }
  
  void SimplePairlist::setCutOff(double c)
  {
    d_cut2 = c*c;
  }
  
  void SimplePairlist::setType(std::string s)
  {
    if(s=="ATOMIC")
      d_chargeGroupBased = false;
    else if(s=="CHARGEGROUP")
      d_chargeGroupBased = true;
    else
      throw gromos::Exception("SimplePairlist", 
			      "Pairlist type " + s + " unknown");
  }

  void SimplePairlist::setAtom(int m, int a)
  {
    if(m<0){
      d_mol=-1;
      d_atom=a;
      return;
    }
    if(m<sys()->numMolecules()){
      d_mol=m;
      if(a<sys()->mol(m).numAtoms())
	d_atom=a;
      else
	throw gromos::Exception("SimplePairlist",
				"Not enough atoms in molecule");
    }
    else
      throw gromos::Exception("SimplePairlist",
			      "Not enough molecules in system");
  }

  void SimplePairlist::calc()
  {
    if(d_chargeGroupBased)
      calcCgb();
    else
      calcAtomic();
  }

  void SimplePairlist::calcCgb()
  {
    if(d_atom < 0)
      throw gromos::Exception("SimplePairlist",
			      "No atom set");
    gmath::Vec atom_i;
    if(d_mol<0 && d_atom >= sys()->sol(0).numPos())
      throw gromos::Exception("SimplePairlist",
				"Not enough solvent atoms in system");
    
    atom_i = chargeGroupPosition(d_mol, d_atom);
    
    // gather the system to get the charge groups connected
    d_pbc->gather();

    // now loop over all charge gropus and add those atoms that belong to
    // a charge group that is within d_cut2
    double d2;
    gmath::Vec v;
    
    //first solute
    for(int m=0; m<sys()->numMolecules(); m++){
      int a1=-1;
      int a2=0;
      
      
      while(a2<sys()->mol(m).numAtoms()){
	// search for the first chargeGroup==1
	for(; sys()->mol(m).topology().atom(a2).chargeGroup()!=1; a2++);
	v=chargeGroupPosition(m,a2);
	v=d_pbc->nearestImage(atom_i, v, sys()->box());
	d2=(atom_i - v).abs2();
	if(d2 <=d_cut2)
	  // now add everything from a1+1 to a2;
	  for(int i=a1; i<a2; i++) addAtom(m,i+1);
	a1 = a2;
	a2++;
      }
    }

    //and solvent
    int nsa = sys()->sol(0).topology().numAtoms();
    for(int i=0; i< sys()->sol(0).numPos(); i+=nsa){
      v = d_pbc->nearestImage(atom_i, sys()->sol(0).pos(i), sys()->box());
      d2 = (atom_i - v).abs2();
      if(d2 <= d_cut2)
	for(int j=0; j<nsa; j++) addAtom(-1,i+j);
    }
    
    // and remove the atom itself
    removeAtom(d_mol, d_atom);
  }
  
  void SimplePairlist::calcAtomic()
  {
    if(d_atom < 0)
      throw gromos::Exception("SimplePairlist",
			      "No atom set");
    gmath::Vec atom_i;
    if(d_mol<0)
      if(d_atom < sys()->sol(0).numPos())
	atom_i = sys()->sol(0).pos(d_atom);
      else
	throw gromos::Exception("SimplePairlist",
				"Not enough solvent atoms in system");
    else
      atom_i = sys()->mol(d_mol).pos(d_atom);

    // now loop over all atoms and add those atoms that are within d_cut2
    double d2;
    gmath::Vec v;
    
    // first solute
    for(int m=0; m<sys()->numMolecules(); m++){
      for(int a=0; a<sys()->mol(m).numPos(); a++){
	v=d_pbc->nearestImage(atom_i, sys()->mol(m).pos(a), sys()->box());
	d2=(atom_i - v).abs2();
	if(d2<=d_cut2)
	  addAtom(m,a);
      }
    }
    // and solvent
    for(int i=0; i<sys()->sol(0).numPos(); i++){
      v=d_pbc->nearestImage(atom_i, sys()->sol(0).pos(i), sys()->box());
      d2=(atom_i - v).abs2();
      if(d2<=d_cut2)
	addAtom(-1,i);
    }
    // now remove the atom itself
    removeAtom(d_mol, d_atom);
  }
  
  gmath::Vec SimplePairlist::chargeGroupPosition(int m, int a)
  {
    gmath::Vec v(0.0,0.0,0.0);
    if(m<0){
      int i=a/sys()->sol(0).topology().numAtoms();
      i*=sys()->sol(0).topology().numAtoms();
      return sys()->sol(0).pos(i);
    }
    
    int begin=a-1, end=a;
    if(a>0)
      for(begin=a-1; 
	  begin>=0 
	    && sys()->mol(m).topology().atom(begin).chargeGroup()!=1; 
	  begin--);
    for(end=a; 
	sys()->mol(m).topology().atom(end).chargeGroup()!=1; 
	end++);
    
    // charge group goes from begin+1 to end
    for(int k=begin+1; k<=end; k++)
      v += sys()->mol(m).pos(k);
    return v/(end-begin);
  }

  void SimplePairlist::removeExclusions()
  {
    // of course it is not effective to first add and later remove
    // the excluded atoms, but if you only want a list of atoms within a
    // cutoff then you just do not call this function

    // check whether we are looking at a solvent
    if(d_mol<0){
      int nsa=sys()->sol(0).topology().numAtoms();
      int first=d_atom/nsa;
      first *= nsa;
      for(int i=0; i<nsa; i++) removeAtom(-1, first+i);
    }
    else{
      
    // loop over all solute atoms before d_atom
      for(int a=0; a<d_atom; a++)
	for(int i=0; 
	    i< sys()->mol(d_mol).topology().atom(a).exclusion().size();
	    i++)
	  if(d_atom == 
	     sys()->mol(d_mol).topology().atom(a).exclusion().atom(i))
	    removeAtom(d_mol,a);
      // and remove all excluded atoms of d_a
      for(int i=0; 
	  i < sys()->mol(d_mol).topology().atom(d_atom).exclusion().size();
	  i++)
	removeAtom(d_mol, 
	   sys()->mol(d_mol).topology().atom(d_atom).exclusion().atom(i));
    }
  }

  void SimplePairlist::remove14Exclusions()
  {
    // of course it is not effective to first add and later remove
    // the excluded atoms, but if you only want a list of atoms within a
    // cutoff then you just do not call this function

    // check whether we are looking at a solvent
    if(d_mol<0){
      return;
    }
    else{
      
    // loop over all solute atoms before d_atom
      for(int a=0; a<d_atom; a++)
	for(int i=0; 
	    i< sys()->mol(d_mol).topology().atom(a).exclusion14().size();
	    i++)
	  if(d_atom == 
	     sys()->mol(d_mol).topology().atom(a).exclusion14().atom(i))
	    removeAtom(d_mol,a);
      // and remove all excluded atoms of d_a
      for(int i=0; 
	  i < sys()->mol(d_mol).topology().atom(d_atom).exclusion14().size();
	  i++)
	removeAtom(d_mol, 
	   sys()->mol(d_mol).topology().atom(d_atom).exclusion14().atom(i));
    }
  }

}



     
