// pb_FFTChargeShapingFunction.h

#ifndef INCLUDED_PB_FFTChargeShapingFunction
#define INCLUDED_PB_FFTChargeShapingFunction
#ifndef INCLUDED_PB_PB_Parameters
#include "PB_Parameters.h"
#endif

namespace pb{


class FFTChargeShapingFunction{
public:
	
	int HatType;
	int ParabolaType;
//	int OtherParabolaType;
     
        double tinynum;
        pb::PB_Parameters ppp;


public:
    //constructor
      FFTChargeShapingFunction();
    //deconstructor
     virtual ~FFTChargeShapingFunction(){}
    //methods
	/* FFTChargeShapingFunction getCSF(int type);*/
	virtual double calc(double xx, double yy, double zz,
		            double alpha, double eps0){}




}; // class
} // namespace


#endif

