// utils_Rmsd.cc

#include "Rmsd.h"
#include "../fit/Reference.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gmath/Vec.h"
#include "../utils/Property.h"
#include "../utils/PropertyContainer.h"
#include <vector>

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


double Rmsd::rmsdproperty(const System &sys){

  double rmsd2=0;
  PropertyContainer prop_ref = *d_prop_ref;
  PropertyContainer prop_sys = *d_prop_sys;

  for(unsigned int i=0; i< d_prop_ref->size(); i++) rmsd2 += pow((prop_ref[i]->calc() - prop_sys[i]->calc()),2);

  return sqrt(rmsd2);
}

void Rmsd::addproperty(const PropertyContainer *prop_ref, const PropertyContainer *prop_sys) {  
  d_prop_ref = prop_ref;
  d_prop_sys = prop_sys;
}
