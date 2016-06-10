// pb_FFTVacuumField.h

#ifndef INCLUDED_PB_FFTVacuumField
#define INCLUDED_PB_FFTVacuumField
#ifndef INCLUDED_PB_PB_Parameters
#include "PB_Parameters.h"
#endif
#ifndef INCLUDED_PB_FFTGridType
#include "FFTGridType.h"
#endif
#ifndef INCLUDED_PB_FFTBoundaryCondition
#include "FFTBoundaryCondition.h"
#endif

namespace pb{


class FFTVacuumField{
protected:
	utils::AtomSpecifier atoms;
	
	FFTGridType gt;
	FFTBoundaryCondition bc;

        PB_Parameters ppp;
        double tinynum;
        double csfunc;
        double pi2;
        double eps0;
        double fpepsi;
        
public:
    //constructor
    FFTVacuumField(utils::AtomSpecifier atoms, FFTGridType gt, FFTBoundaryCondition bc);
    //deconstructor
     virtual  ~FFTVacuumField(){}
    //methods
	
   /* FFTVacuumField getVF(
			int type,
			utils::AtomSpecifier atoms,
			FFTGridType gt, FFTBoundaryCondition bc);
*/
	virtual void calcVacField(
			std::vector <double> &  Vx, std::vector <double> & Vy, std::vector <double> &  Vz){}



}; // class
} // namespace


#endif

