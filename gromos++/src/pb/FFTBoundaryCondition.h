// pb_FFTBoundaryCondition.h

#ifndef INCLUDED_PB_FFTBoundaryCondition
#define INCLUDED_PB_FFTBoundaryCondition
#ifndef INCLUDED_PB_PB_Parameters
#include "PB_Parameters.h"
#endif


namespace pb{



class FFTBoundaryCondition{

public:


        pb::PB_Parameters ppp;

        double tinynum;
        
	int type;
	
	std::string stype; 
	/* 
	  Ewald parameter alpha1
	  width of the charge-shaping function
	 */
	double alpha1;
	/*
	  Ewald parameter alpha2
	  width of the charge-shaping function
	 */
	double alpha2;


	/*
	  Ewald parameter charge_shape
	  1 : hat function
          2 : parabola function
          3 : other parabola function
	 
	int charge_shape;
         */

	/*
	  Ewald parameter nalias1
	  number of alias vectors for vacuum potential
	 */

	int nalias1;
	/* 
	  Ewald parameter nalias2
	  number of alias vectors for vacuum potential
	 */
	int nalias2;     
	/*
	  Ewald parameter rd_field
	  controls reading k-space potntial from file
	 
	boolean rd_field = false;*/

	/*
	  cutoff radius
	 */
	double cutoff;          


 /* epsilon (for checking) */
        double eps;

 public:
  // constructor

     FFTBoundaryCondition(int type, std::string stype, double alpha1, double alpha2, int nalias1, int nalias2, double cutoff , double eps );
     FFTBoundaryCondition(){}

   // deconstructor
    ~FFTBoundaryCondition(){}


  // methods


	void dumpparameters(); 

}; // class
} // namespace


#endif