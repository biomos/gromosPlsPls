// pb_FFTChargeDipole_RF.cc


#include <new>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <set>
#include "../fit/PositionUtils.h"
#include "../utils/AtomSpecifier.h"
#include "../gmath/Physics.h"
#include "../gcore/System.h"


#include "FFTChargeDipole_RF.h"
#include "FFTDipoleDipole_RF.h"
#include "FFTChargeDipole.h"



using pb::FFTChargeDipole_RF;
using pb::FFTChargeDipole;
using pb::FFTDipoleDipole_RF;

FFTChargeDipole_RF::FFTChargeDipole_RF(
			double epsilonBoundary,
			double cutoff,double epssolvent) : FFTChargeDipole(epssolvent) {


                this->epsB=epsilonBoundary;
		this->cut = cutoff;


                 std::cout << "# FFTChargeDipole_RF: epssolvent " << epssolvent << endl;
                 std::cout << "# FFTChargeDipole_RF: epsB " << epsB << endl;;


                try{
		if (epsB < 1.0 &&   (fabs(epsB) >  tinynum) ) {
			  throw gromos::Exception("FFTChargeDipole_RF","Invalid permittivity for boundary (should be 0, 1 or >1) ...");

		                }// endif
}// end of try

                         
 catch (const gromos::Exception &e){
                                         cerr << e.what() << endl;
                                         exit(1);
                                }       
                
                
                if (fabs(epsB)<tinynum) {
			this->first = 1.0;
			this->second = 0.0;
		}
		else {
			double twoEpsPrimePlusOne = 2 * epsB + 1;
			double twoEpsPrimeMinusOne = 2 * (epsB - 1);
			
			this->first = twoEpsPrimeMinusOne / twoEpsPrimePlusOne;
			this->second = 3 / twoEpsPrimePlusOne;
		}			

}


		
		
	
	
	/* kr the product of the norm of the k-vector and the cutoff distance */

	double FFTChargeDipole_RF::fkr(double kr) {
            double res;
            //return FFTDipoleDipole_RF.fkr(kr);


		double kri = 1 / kr;
		double kr2i = kri * kri;

		return 3 * kr2i * ( - cos(kr) + sin(kr) * kri);



	}
	
	/*  we assume that the A-factor (and therefore fckr) is up-to-date */
	 
	double FFTChargeDipole_RF::polarization(double k2) {
		
		double kr = sqrt(k2) * cut;		
		
		return 1.0 - first * fkr(kr) - second * sin(kr) / kr;
	}



