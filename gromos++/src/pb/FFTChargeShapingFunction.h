// pb_FFTChargeShapingFunction.h

#ifndef INCLUDED_PB_FFTChargeShapingFunction
#define INCLUDED_PB_FFTChargeShapingFunction
#ifndef INCLUDED_PB_PB_Parameters
#include "PB_Parameters.h"
#endif

namespace pb{


class FFTChargeShapingFunction{
public:
  pb::PB_Parameters ppp;	
	int HatType;
	int ParabolaType;
//	int OtherParabolaType;
     
        double tinynum;



public:
    //constructor
      FFTChargeShapingFunction(ofstream &os);
    //deconstructor
     virtual ~FFTChargeShapingFunction(){}
    //methods
	/* FFTChargeShapingFunction getCSF(int type);*/
	virtual double calc(double xx, double yy, double zz,
		            double alpha, double eps0) {}




}; // class
} // namespace


#endif

