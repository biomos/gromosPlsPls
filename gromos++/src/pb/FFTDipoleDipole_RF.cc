// pb_FFTDipoleDipole_RF.cc


#include <new>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <set>
#include "../fit/PositionUtils.h"
#include "../utils/AtomSpecifier.h"
#include "../gmath/Physics.h"
#include "../gcore/System.h"


#include "FFTDipoleDipole_RF.h"
#include "FFTDipoleDipole_LS.h"
#include "FFTDipoleDipole.h"



using pb::FFTDipoleDipole_RF;
using pb::FFTDipoleDipole;

FFTDipoleDipole_RF::FFTDipoleDipole_RF(
			double epsilonBoundary,
			double cutoff, double epsS): FFTDipoleDipole(epsS,epsilonBoundary), ew(epsilonBoundary,epsS){
	
/*    FFTDipoleDipole_RF::FFTDipoleDipole_RF(
			double epsilonBoundary,
			double cutoff){*/

	
                this->epsB=epsilonBoundary;
		this->twoEpsPrimePlusOne = 2 * epsB + 1;
		this->epsSminusOne = epsS - 1;
		this->upper = epsSminusOne * twoEpsPrimePlusOne * twoEpsPrimePlusOne;
		this->twoEs1ePrime = 2 * (epsS - 1) * epsB;
		this->esTwoEprimePlusOne = epsS * twoEpsPrimePlusOne;
		this->cut = cutoff;

            //    if (ppp.get_debugvar()==1){
                    std::cout << "# FFTDipoleDipole_RF : epsB " << epsB << endl;
                    std::cout << "# FFTDipoleDipole_RF : epsS " << epsS << endl;
              //  }


		// Vincent: we keep one of those around, because
		// their expressions for the k0-vector are the same

		//ew = new pb::FFTDipoleDipole_LS(epsB);
	}
	
	
	double FFTDipoleDipole_RF::fkr(double kr) {
		
		double kri = 1 / kr;
		double kr2i = kri * kri;
		
		return 3 * kr2i * ( - cos(kr) + sin(kr) * kri);
	}
	
	double FFTDipoleDipole_RF::computeAFactor(double k2) {
		fckr = fkr(sqrt(k2) * cut);
		if (0.0 == epsB) {
			// conducting boundary
			return 1.0;
		} else {
			return twoEpsPrimePlusOne / (twoEpsPrimePlusOne + epsSminusOne * fckr);
		}
	}
	
	double FFTDipoleDipole_RF::computeBFactor(double k2) {
		// this assumes computeAFactor(k2) has already been called
		// and fckr is up-to-date. ugly, i agree.
		
		if (0.0 == epsB) {
			// conducting boundary
			return - epsSminusOne * (1 - fckr) / (1 + epsSminusOne * (1 - fckr));
		}
		else {
			return - upper * (fckr - 1) / (
					(twoEpsPrimePlusOne + epsSminusOne * fckr) *
					(twoEs1ePrime * fckr - esTwoEprimePlusOne)
			);
		}
	}
	
	void FFTDipoleDipole_RF::updateTensorK0(  double  (& tensor)[3][3])  {
		ew.updateTensorK0(tensor);			
	}






