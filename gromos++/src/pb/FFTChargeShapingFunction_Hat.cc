// pb_FFTChargeShaping_Function_Hat.cc


#include <new>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <set>
#include "../fit/PositionUtils.h"
#include "../utils/AtomSpecifier.h"
#include "../gmath/Physics.h"
#include "../gcore/System.h"


#include "FFTChargeShapingFunction_Hat.h"
#include "FFTChargeShapingFunction.h"



using pb::FFTChargeShapingFunction_Hat;
using pb::FFTChargeShapingFunction;

FFTChargeShapingFunction_Hat::FFTChargeShapingFunction_Hat():FFTChargeShapingFunction(){
// FFTChargeShapingFunction_Hat::FFTChargeShapingFunction_Hat(){
}

	

		
		 /* k-space version of the potential due to
		 a spherical hat function at the origin plus homogeneous background charge 
		 Phi(k) = 1/(eps0*k^2) * 12(2 - 2 Math.cos ak - ak Math.sin ak)/(ak)^4 */
		 
	double FFTChargeShapingFunction_Hat::calc(
				double x, double y, double z, 
				double alpha, double eps0) {
			
			
		double k2 = x*x + y*y + z*z;
			if (k2 >= tinynum) {
			double alphak = alpha*sqrt(k2);
				return 1.0/(eps0*k2)*12.0*(2.0-2.0*cos(alphak)-alphak*sin(alphak))/(alphak*alphak*alphak*alphak);
			} else {
				return 0.0;
			}
		}
	


