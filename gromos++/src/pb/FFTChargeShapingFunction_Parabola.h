// pb_FFTChargeShapingFunction_Parabola.h

#ifndef INCLUDED_PB_FFTChargeShapingFunction_Parabola
#define INCLUDED_PB_FFTChargeShapingFunction_Parabola
#ifndef INCLUDED_PB_FFTChargeShapingFunction
#include "FFTChargeShapingFunction.h"
#endif


namespace pb{



class FFTChargeShapingFunction_Parabola: virtual public FFTChargeShapingFunction{


public:
        //constructor

	FFTChargeShapingFunction_Parabola();


        //deconstructor

        ~FFTChargeShapingFunction_Parabola(){};


        // methods

	double calc(double x, double y, double z,
				double alpha, double eps0);



}; // class
} // namespace


#endif
