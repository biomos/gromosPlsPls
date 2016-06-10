// pb_FFTDipoleDipole.h

#ifndef INCLUDED_PB_FFTDipoleDipole
#define INCLUDED_PB_FFTDipoleDipole
#ifndef INCLUDED_PB_PB_Parameters
#include "PB_Parameters.h"
#endif


namespace pb{



class FFTDipoleDipole{

protected:
    
double tinynum;
double epsS;
double epsB;
pb::PB_Parameters ppp;

public:

    //constructor
     FFTDipoleDipole(double epsilonSolvent, double epsilonB);

    // deconstructor
   virtual ~FFTDipoleDipole(){}

    //methods
    
    
     /* Determine the inverse k-space dipole-dipole interaction
	  tensor for the k=0-vector. you shouldn't need this directly*/

	virtual void updateTensorK0(double (& tensor)[3][3]){}


	 /* determine the scalar multiplier A of the identity matrix
	  component of the interaction tensor. you shouldn't need
	  this directly*/

	virtual double computeAFactor(double k2){}

	/* determine the scalar multiplier B of the outer product
	  component of the interaction tensor. you shouldn't need
	  this directly*/

	virtual double computeBFactor(double k2){}


	/* Determine the inverse k-space dipole-dipole interaction
	 tensor*/

	void updateTensor(		double k2,     double  (& tensor)[3][3]
	) ;
    
    
    
    
}; // class
} // namespace


#endif
