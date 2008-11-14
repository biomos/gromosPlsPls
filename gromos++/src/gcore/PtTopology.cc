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
    d_mass.resize(p);
    d_charge.resize(p);
    d_pertnames.resize(p, "STATE");
    d_polarisability.resize(p);
    d_dampingLevel.resize(p);
    
    for(int i=0;i<p;i++){
      d_iac[i].resize(a, 0);
      d_mass[i].resize(a, 0.0);
      d_charge[i].resize(a, 0.0);
      d_polarisability[i].resize(a, 0.0);
      d_dampingLevel[i].resize(a, 0.0);
    }
    
    d_atomnames.resize(a, "AT");
    d_atomnum.resize(a, 0);
    d_alphaLJ.resize(a, 0.0);
    d_alphaCRF.resize(a, 0.0);
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
            sys.mol(m).topology().atom(aa).setMass(mass(counter, iipt));
            if (sys.mol(m).topology().atom(aa).isPolarisable() && hasPolarisationParameters()) {
              sys.mol(m).topology().atom(aa).setPolarisability(polarisability(counter, iipt));
              sys.mol(m).topology().atom(aa).setDampingLevel(dampingLevel(counter, iipt));
            }
	    
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
  
  void PtTopology::setMass(int a, int p, double q)
  {
    d_mass[p][a]=q;
  }
  
  void PtTopology::setPolarisability(int a, int p, double q)
  {
    d_polarisability[p][a]=q;
  }
  
  void PtTopology::setDampingLevel(int a, int p, double q)
  {
    d_dampingLevel[p][a]=q;
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
  
  void PtTopology::setAlphaLJ(int a, double alpha)
  {
    d_alphaLJ[a]=alpha;
  }
  
  void PtTopology::setAlphaCRF(int a, double alpha)
  {
    d_alphaCRF[a]=alpha;
  }
  
  void PtTopology::setHasPolarisationParameters(bool b)
  {
    d_hasPolaristaionParams=b;
  }

}



