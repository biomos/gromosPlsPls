// pb_FFTDipoleDipole_RF.h

#ifndef INCLUDED_PB_FFTDipoleDipole_RF
#define INCLUDED_PB_FFTDipoleDipole_RF
#ifndef INCLUDED_PB_FFTDipoleDipole
#include "FFTDipoleDipole.h"
#endif
#ifndef INCLUDED_PB_FFTDipoleDipole_LS
#include "FFTDipoleDipole_LS.h"
#endif

namespace pb{



class FFTDipoleDipole_RF: virtual public FFTDipoleDipole{

	
	double twoEpsPrimePlusOne;
	double epsSminusOne;
	double upper;
	double twoEs1ePrime; 
	double esTwoEprimePlusOne;
	double cut;
	double fckr;
	
	FFTDipoleDipole_LS ew;

public:
    // constructor
	 FFTDipoleDipole_RF(
			double epsilonBoundary, 
			double cutoff, double epsS);
    // deconstructor
         ~FFTDipoleDipole_RF(){}

	
	double fkr(double kr);

	
	double computeAFactor(double k2);
	
	double computeBFactor(double k2);
	
	void updateTensorK0( double  (& tensor)[3][3] );



}; // class
} // namespace


#endif
