// pb_FFTDipoleDipole.cc


#include <new>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <set>
#include "../fit/PositionUtils.h"
#include "../utils/AtomSpecifier.h"
#include "../gmath/Physics.h"
#include "../gcore/System.h"


#include "FFTDipoleDipole.h"
#include "PB_Parameters.h"




using pb::FFTDipoleDipole;

FFTDipoleDipole::FFTDipoleDipole(double epsilonSolvent, double epsilonB, ofstream &os):ppp(epsilonSolvent, os){
	
  this->tinynum = ppp.getTiny_real();
  this->epsS = epsilonSolvent;
  this->epsB=epsilonB;
                

  os << "# FFTDipoleDipole: epsS " << epsS << endl;
  os << "# FFTDipoleDipole: epsB " << epsB << endl;;

                 
  try{
    if (  (epsB < 1.0 && fabs(epsB) >  tinynum)   )  {
      throw gromos::Exception("FFTDipoleDipole","Invalid permittivity for boundary (should be 0, 1 or >1) ...");
    }// endif
  }// end of try

  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }

}
	

	
/* Determine the inverse k-space dipole-dipole interaction
   tensor*/
	
void FFTDipoleDipole::updateTensor(
				   double k2,
				   double (& tensor)[3][3],
				   ofstream &os
				   ) {


  // os << "# FFTDipoleDipole::updateTensor ... " << endl;
		
  if (k2 < tinynum) {
    // this is the zero-vector
    updateTensorK0(tensor);
    //    os << "# was zero vec ..." << endl;
    return;
  }
		
  double A = computeAFactor(k2);
		
  double B = computeBFactor(k2) / k2;
  if (ppp.get_debugvar()==1){
    os << "# A and B " << A << " " << B << endl;
  }

  for (int ii = 0; ii < 3; ii++){
    for (int jj = 0; jj < 3; jj++){
      tensor[ii][jj] *= B;
      if (ppp.get_debugvar()==1){
	os << "# Bfac ... tensor [ " << ii << " ][ " << jj << "] = " << tensor[ii][jj] << endl;
      }
    }}
		
  for (int ii = 0; ii < 3; ii++){
    tensor[ii][ii] += A;
    if (ppp.get_debugvar()==1){
      os << "# Afac ... tensor [ " << ii << " ][ " << ii << "] = " << tensor[ii][ii] << endl;
    }
  }
                
                
}

