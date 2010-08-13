// utils_Rmsd.cc

#include <cassert>
#include <sstream>
#include <vector>
#include <set>
#include "../fit/Reference.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/LJException.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gmath/Vec.h"
#include "../gmath/Stat.h"
#include "../utils/AtomSpecifier.h"
#include "../utils/Value.h"
#include "../utils/Property.h"
#include "../utils/PropertyContainer.h"
#include "Rmsd.h"

using gcore::System;
using utils::Rmsd;
using fit::Reference;

using namespace std;

Rmsd::Rmsd(const Reference *ref){
  d_ref=ref;
}



double Rmsd::rmsd(const System &sys){
  double rmsd2=0;
  for(int m=0;m<sys.numMolecules();++m)
    for(int n=0;n<sys.mol(m).numAtoms();++n)
      if(d_ref->weight(m,n))
	rmsd2+=d_ref->weight(m,n)*
	  (d_ref->sys().mol(m).pos(n) 
	   - sys.mol(m).pos(n)).abs2();

  return sqrt(rmsd2);
}


/**
 * @TODO should include periodicity. But this depends on the property?
 */
double Rmsd::rmsdproperty(const System &sys){

  double rmsd2=0;
  PropertyContainer prop_ref = *d_prop_ref;
  PropertyContainer prop_sys = *d_prop_sys;

  for(unsigned int i=0; i < d_prop_ref->size(); i++){
    
    Value res = abs2(prop_ref[i]->calc() - prop_sys[i]->calc());
    rmsd2 += res.scalar();
  }
  
  return sqrt(rmsd2/d_prop_ref->size());
}

void Rmsd::addproperty(const PropertyContainer *prop_ref, const PropertyContainer *prop_sys) {  
  d_prop_ref = prop_ref;
  d_prop_sys = prop_sys;
}
