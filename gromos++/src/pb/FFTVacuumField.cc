// pb_FFTVacuumField.cc


#include <new>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <set>
#include "../fit/PositionUtils.h"
#include "../utils/AtomSpecifier.h"
#include "../gmath/Physics.h"
#include "../gcore/System.h"


#include "FFTVacuumField.h"
#include "FFTBoundaryCondition.h"
#include "FFTGridType.h"
#include "PB_Parameters.h"


using pb::FFTVacuumField;

FFTVacuumField::FFTVacuumField(utils::AtomSpecifier atoms, FFTGridType gt, FFTBoundaryCondition bc){

		this->atoms = atoms;
		this->gt = gt;
		this->bc = bc;
                this->tinynum=ppp.getTiny_real();
                this->csfunc=ppp.get_default_chshape();
                this->pi2=2*ppp.getPI();
                this->eps0=1.0/(ppp.getFPEPSI()*4*ppp.getPI());
                this->fpepsi=ppp.getFPEPSI();
	}
	
	/* FFTVacuumField FFTVacuumField::getVF(
			int type,
			AtomSpecifier atoms,
			FFTGridType gt, FFTBoundaryCondition bc){

               try{
		if (FFTInteractionTypeCodes.lsType == type)
			return new FFTVacuumField_LS(atoms,  gt, bc);

                else if (FFTInteractionTypeCodes.rfType == type)
			return new FFTVacuumField_RF(atoms,  gt, bc);

		else{
		throw gromos::Exception("FFTVacuumField","Invalid interaction type code. Exiting.");
	} // endif
    } //end of try

                                catch (const gromos::Exception &e){
                                         cerr << e.what() << endl;
                                         exit(1);
                                }

}*/
