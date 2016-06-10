// pb_FFTChargeDipole_LS.h

#ifndef INCLUDED_PB_FFTChargeDipole_LS
#define INCLUDED_PB_FFTChargeDipole_LS
#ifndef INCLUDED_PB_FFTChargeDipole
#include "FFTChargeDipole.h"
#endif


namespace pb{



class FFTChargeDipole_LS: virtual public FFTChargeDipole{


public:
        //constructor

	 FFTChargeDipole_LS(double epssolvent);


        //deconstructor

        ~FFTChargeDipole_LS(){}


        // methods

	double polarization(double k2);



}; // class
} // namespace


#endif


