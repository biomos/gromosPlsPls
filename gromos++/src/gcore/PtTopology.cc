//gcore_PtTopology.cc
#include <vector>
#include <string>
#include <sstream>
#include <cassert>
#include "../gromos/Exception.h"
#include "AtomTopology.h"
#include "MoleculeTopology.h"
#include "Molecule.h"
#include "System.h"
#include "PtTopology.h"


namespace gcore
{
  
  void PtTopology::setSize(int a, int p)
  {
    d_iac.resize(p);
    d_charge.resize(p);
    d_pertnames.resize(p);
    
    for(int i=0;i<p;i++){
      d_iac[i].resize(a);
      d_charge[i].resize(a);
    }
    
    d_atomnames.resize(a);
    d_atomnum.resize(a);
  }
  
  void PtTopology::apply(System &sys, int iipt)
  {
    int counter=0;
    int mol, atom;
    
    //determine molecule and atom for this counter (the first)
    findAtom(sys, mol, atom, counter);
    
    // loop over the molecules 
    for(int m=0;m<sys.numMolecules();m++){
      
      if(!(m<mol || counter==numAtoms())){
	
	// loop over the atoms in this molecule
	for(int aa=0;aa<sys.mol(m).topology().numAtoms();aa++){
	  if(aa==atom){
	    if(atomName(counter)!=sys.mol(m).topology().atom(aa).name()){
	      std::ostringstream os;
	      os << "Atom names in (perturbation) topologies do not match\n"
		 << "Topology: " << sys.mol(m).topology().atom(aa).name() 
		 << " (" << m+1 << ":" << aa+1 << ")"
		 << "\tPerturbation topology: " << atomName(counter)
		 << " (" << counter << ")";
	      
	      
	      throw gromos::Exception("PtTopology", os.str());
	    }
	    
	    sys.mol(m).topology().atom(aa).setIac(iac(counter, iipt));
	    sys.mol(m).topology().atom(aa).setCharge(charge(counter, iipt));
	    
	    // loop to next atom in the perturbation list
	    counter++;
	    findAtom(sys, mol, atom, counter);
	  }
	}
      } 
    }
  }
  
  void PtTopology::findAtom(System &sys, int &mol, int &atom, int counter)
  {
    if(counter>=numAtoms()){
      mol=-1;
      atom=-1;
    }
    else{
      int at_cnt=0;
      for(mol=0;at_cnt<=atomNum(counter)&&mol<sys.numMolecules();mol++)
	at_cnt+=sys.mol(mol).topology().numAtoms();
      mol--;
      atom=atomNum(counter) - at_cnt + sys.mol(mol).topology().numAtoms();
    }
  }
  
  void PtTopology::setIac(int a, int p, int iac)
  {
    d_iac[p][a]=iac;
  }

  void PtTopology::setCharge(int a, int p, double q)
  {
    d_charge[p][a]=q;
  }
  
  void PtTopology::setAtomName(int a, std::string name)
  {
    d_atomnames[a]=name;
  }
  void PtTopology::setPertName(int p, std::string name)
  {
    d_pertnames[p]=name;
  }

  void PtTopology::setAtomNum(int a, int num)
  {
    d_atomnum[a]=num;
  }

}



