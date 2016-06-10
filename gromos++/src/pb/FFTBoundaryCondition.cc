// pb_FFTBoundaryCondition.cc


#include <new>
#include <iostream>
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
        double alpha1, double alpha2, int nalias1, int nalias2, double cutoff, double epsRF  ){

    this->tinynum = ppp.getTiny_real();
    
    this->type=type;
    this->stype=stype;
    this->alpha1=alpha1;
    this->alpha2=alpha2;
    this->nalias1=nalias1;
    this->nalias2=nalias2;
    this->cutoff=cutoff;

    this->eps=epsRF;




          /*      try{
                    FFTInteractionTypeCodes iii;
		if (  (
                       fabs(eps-1)   < tinynum  )  &&  (type == iii.rfType)  ) {
			  throw gromos::Exception("FFTBoundaryCondition","You want RF but have a permittivity of 1. Choose SC for eps==1.");
                    }// endif
}// end of try

 catch (const gromos::Exception &e){
                                         cerr << e.what() << endl;
                                         exit(1);
                                }

*/

}


	void FFTBoundaryCondition::dumpparameters() {
		std::cout << "# BOUNDARY PARAMETERS" << endl;
		std::cout << "# -------------------" << endl;
                std::cout << "# TYPE " << type << endl;
                std::cout << "# STYPE " << stype << endl;
		std::cout << "# HAT CHARGE SHAPING FUNCTION" << endl;
		std::cout << "# ALPHA1: " <<  alpha1 << endl;
		std::cout << "# ALPHA2: " << alpha2 << endl;
		std::cout << "# NALIAS1: " << nalias1 << endl;
		std::cout << "# NALIAS2: " << nalias2 << endl;
        }