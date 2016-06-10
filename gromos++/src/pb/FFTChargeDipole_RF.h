// pb_FFTChargeDipole_RF.h

#ifndef INCLUDED_PB_FFTChargeDipole_RF
#define INCLUDED_PB_FFTChargeDipole_RF
#ifndef INCLUDED_PB_FFTChargeDipole
#include "FFTChargeDipole.h"
#endif


namespace pb{



class FFTChargeDipole_RF: virtual public FFTChargeDipole{


    double cut;
    double epsB;
    double first;
    double second;

	
public:
    // constructor
	 FFTChargeDipole_RF(
			double epsilonBoundary,
			double cutoff,double epssolvent);
    // deconstructor
         ~FFTChargeDipole_RF(){}


	double fkr(double kr);
	double polarization(double k2);





}; // class
} // namespace


#endif
