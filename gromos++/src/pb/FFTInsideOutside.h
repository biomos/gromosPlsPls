// pb_FFTInsideOutside.h

#ifndef INCLUDED_PB_FFTInsideOutside
#define INCLUDED_PB_FFTInsideOutside
#ifndef INCLUDED_PB_PB_Parameters
#include "PB_Parameters.h"
#endif



namespace pb{



  class FFTInsideOutside{


    double tinynum;
    int ngridpoints;
    //std::vector<double>  in;



  public:
    //constructor
    FFTInsideOutside(int ngridpoints, ofstream &os);



    // deconstructor
    ~FFTInsideOutside(){}

    //methods

    /* determine the grid points in the solute,
       This is identical to inside_sol_orig, but
       the outer loop is the loop over the solute atoms, this way, only grid points
       which may actually be inside the solute need to be investigated, greatly speeding
       up the algorithm */


    void inside_sol(  int gridN[3],  double  gridD[3],
		      int nc,
		      utils::AtomSpecifier & atoms, std::vector<double> & in);


    /* Integrates the inside grid and sets values to 1 or 0,
       computes integral of inside and boundary and returns nr of points inside */
    void integr_inside(std::vector <double>  &  inside, ofstream &os);


  }; // class
} // namespace


#endif
