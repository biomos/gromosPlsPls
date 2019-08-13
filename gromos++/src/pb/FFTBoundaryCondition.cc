// pb_FFTBoundaryCondition.cc


#include <new>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <set>
#include "../fit/PositionUtils.h"
#include "../utils/AtomSpecifier.h"
#include "../gmath/Physics.h"
#include "../gcore/System.h"


#include "FFTBoundaryCondition.h"
#include "FFTInteractionTypeCodes.h"


using pb::FFTBoundaryCondition;

FFTBoundaryCondition::FFTBoundaryCondition(int type, std::string stype,
					   double alpha1, double alpha2, int nalias1, int nalias2, double cutoff, double epsRF, ofstream &os):ppp(os){
  // from public top
  this->tinynum = ppp.getTiny_real();
    
  this->type=type;
  this->stype=stype;
  this->alpha1=alpha1;
  this->alpha2=alpha2;
  this->nalias1=nalias1;
  this->nalias2=nalias2;
  this->cutoff=cutoff;

  this->eps=epsRF;

}



void FFTBoundaryCondition::dumpparameters(ofstream &os) {
  os << "# BOUNDARY PARAMETERS" << endl;
  os << "# -------------------" << endl;
  os << "# TYPE " << type << endl;
  os << "# STYPE " << stype << endl;
  os << "# HAT CHARGE SHAPING FUNCTION" << endl;
  os << "# ALPHA1: " <<  alpha1 << endl;
  os << "# ALPHA2: " << alpha2 << endl;
  os << "# NALIAS1: " << nalias1 << endl;
  os << "# NALIAS2: " << nalias2 << endl;
}
