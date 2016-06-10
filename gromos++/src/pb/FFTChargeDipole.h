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
  FFTChargeDipole(double epssolvent);

    // deconstructor
  virtual ~FFTChargeDipole(){}

    // methods
   
  virtual double polarization(double k2){}
   
	




}; // class
} // namespace


#endif




