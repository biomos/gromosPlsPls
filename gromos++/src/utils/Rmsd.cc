// utils_Rmsd.cc

#include "Rmsd.h"
#include "../fit/Reference.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gmath/Vec.h"
#include <vector>

using gcore::System;
using utils::Rmsd;
using fit::Reference;

Rmsd::Rmsd(const Reference *ref){
  d_ref=ref;
}

double Rmsd::rmsd(const System &sys)const{
  double rmsd2=0;
  for(int m=0;m<sys.numMolecules();++m)
    for(int n=0;n<sys.mol(m).numAtoms();++n)
      if(d_ref->weight(m,n))
	rmsd2+=d_ref->weight(m,n)*
	  (d_ref->sys().mol(m).pos(n) 
	   - sys.mol(m).pos(n)).abs2();

  return sqrt(rmsd2);
}
