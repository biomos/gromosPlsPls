// pb_FFTChargeShaping_Function_Parabola.cc


#include <new>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <set>
#include "../fit/PositionUtils.h"
#include "../utils/AtomSpecifier.h"
#include "../gmath/Physics.h"
#include "../gcore/System.h"


#include "FFTChargeShapingFunction_Parabola.h"
#include "FFTChargeShapingFunction.h"



using pb::FFTChargeShapingFunction_Parabola;
using pb::FFTChargeShapingFunction;

FFTChargeShapingFunction_Parabola::FFTChargeShapingFunction_Parabola():FFTChargeShapingFunction(){
//FFTChargeShapingFunction_Parabola::FFTChargeShapingFunction_Parabola(){
}

	

		

		
		 /* k-space version of the potential due to
		  a spherical parabola function plus homogeneous background charge */
		 
	double FFTChargeShapingFunction_Parabola::calc(
				double x, double y, double z,
				double alpha, double eps0) {
                        double k2 = x*x + y*y + z*z;
			if (k2 >= tinynum) {
				double alphak = alpha*sqrt(k2);
				return 1.0/(eps0*k2)*15.0*(-3.0*alphak*cos(alphak)
						+(3.0-alphak*alphak)*sin(alphak))
						/(alphak*alphak*alphak*alphak*alphak);
			} else {
				return 0.0;
			}
		}
	
        