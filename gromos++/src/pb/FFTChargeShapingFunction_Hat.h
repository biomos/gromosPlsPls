// pb_FFTChargeShapingFunction_Hat.h

#ifndef INCLUDED_PB_FFTChargeShapingFunction_Hat
#define INCLUDED_PB_FFTChargeShapingFunction_Hat
#ifndef INCLUDED_PB_FFTChargeShapingFunction
#include "FFTChargeShapingFunction.h"
#endif


namespace pb{



class FFTChargeShapingFunction_Hat: virtual public FFTChargeShapingFunction{


public:
        //constructor

	FFTChargeShapingFunction_Hat();


        //deconstructor

        ~FFTChargeShapingFunction_Hat(){};


        // methods

	double calc(double x, double y, double z,
				double alpha, double eps0);



}; // class
} // namespace


#endif
	