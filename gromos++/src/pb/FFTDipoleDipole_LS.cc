// pb_FFTDipoleDipole_LS.cc


#include <new>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <set>
#include "../fit/PositionUtils.h"
#include "../utils/AtomSpecifier.h"
#include "../gmath/Physics.h"
#include "../gcore/System.h"


#include "FFTDipoleDipole_LS.h"
#include "FFTDipoleDipole.h"



using pb::FFTDipoleDipole_LS;
using pb::FFTDipoleDipole;

FFTDipoleDipole_LS::FFTDipoleDipole_LS(double epsBoundary,double epsS): FFTDipoleDipole(epsS,epsBoundary){
//FFTDipoleDipole_LS::FFTDipoleDipole_LS(double epsBoundary){

		this->es_1_es = (epsS - 1) / epsS;

                this->epsB=epsBoundary;

  try{
		if (epsB < 1.0 &&  (fabs(epsB) >  tinynum)  ) {
			  throw gromos::Exception("FFTDipoleDipole",
                                  "Invalid permittivity for boundary (should be 0, 1 or <1) ...");

		                }// endif
}// end of try
 catch (const gromos::Exception &e){
                                         cerr << e.what() << endl;
                                         exit(1);
                                }
  

		if (fabs(epsB) < tinynum) {
			// conducting boundary
			k0EpsilonFactor = 1.0;
		} else {
			k0EpsilonFactor = (2 * epsB + 1) / (2 * epsB + epsS);
		}




}



	void FFTDipoleDipole_LS::updateTensorK0(double ( & tensor) [3][3]) {
		
		for (int ii = 0; ii < 3; ii++) 
			for (int jj = 0; jj < 3; jj++)
				tensor[ii][jj] = 0.0;
		
		for (int ii = 0; ii < 3; ii++){
			tensor[ii][ii] = k0EpsilonFactor;
                
                }
		
	}
	
	double FFTDipoleDipole_LS::computeAFactor(double k2) {
		return 1.0;
	}
	
	double FFTDipoleDipole_LS::computeBFactor(double k2) {
		return - es_1_es;
	}
