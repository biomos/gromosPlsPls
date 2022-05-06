// pb_FFTChargeDipole.h

#ifndef INCLUDED_PB_FFTChargeDipole
#define INCLUDED_PB_FFTChargeDipole
#ifndef INCLUDED_PB_PB_Parameters
#include "PB_Parameters.h"
#endif


namespace pb{



class FFTChargeDipole{
protected:


double tinynum;
double epssolvent;

 pb::PB_Parameters ppp;

public:

    // constructor
 FFTChargeDipole(double epssolvent, ofstream &os);

    // deconstructor
  virtual ~FFTChargeDipole(){}

    // methods
   
  virtual double polarization(double k2){return 0.0;}
   
	




}; // class
} // namespace


#endif




