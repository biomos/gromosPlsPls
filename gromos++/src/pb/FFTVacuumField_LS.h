// pb_FFTVacuumField_LS.h

#ifndef INCLUDED_PB_FFTVacuumField_LS
#define INCLUDED_PB_FFTVacuumField_LS
#ifndef INCLUDED_PB_FFTVacuumField
#include "FFTVacuumField.h"
#endif
#ifndef INCLUDED_PB_FFTChargeShapingFunction
#include "FFTChargeShapingFunction.h"
#endif
#ifndef INCLUDED_PB_PB_Parameters
#include "PB_Parameters.h"
#endif


namespace pb{



class FFTVacuumField_LS: virtual public FFTVacuumField{
    
        FFTChargeShapingFunction *csfunc;
        PB_Parameters ppp;
        int cstype;

	std::vector <int> ion_list1;
	std::vector <int> ion_list2;
	
	int ion_count1;
	int ion_count2;
	
	
	 /* arrays for storing configuration independent part of the vacuum field */

	std::vector <double> Vx1_store;
	std::vector <double> Vy1_store;
	std::vector <double> Vz1_store;
	std::vector <double> Vx2_store;
	std::vector <double> Vy2_store;
	std::vector <double> Vz2_store;
	
	int ngrdx;
	int ngrdy;
	int ngrdz;
	double dkx;
	double dky;
	double dkz;
	

public:
        //constructor
        FFTVacuumField_LS(utils::AtomSpecifier atoms, FFTGridType gt, FFTBoundaryCondition bc);
       //deconstructor
        ~FFTVacuumField_LS(){}
       //methods

	
	
	void complexFromDouble(
			std::vector <double> & doubleArray, std::vector <double> & complexArray);
		
		
	 /* make lists of "big" and "small" atoms */
	 
	void makeatomlist();
	
	
	void calcVacField(
			std::vector <double> &fldx_k,
                        std::vector <double> & fldy_k,
                        std::vector <double> & fldz_k);

        void positionIndependentVacField(
			std::vector <double> &  destX, std::vector <double> &  destY, std::vector <double> &  destZ,
			double kax, double kay, double kaz,
			int nAlias, double alpha);

        void recyclefield(
			std::vector <double> &  destX, std::vector <double> & destY, std::vector <double> &  destZ	);

	void updateMultipliers(
			std::vector <double> &  multipliers,
			std::vector <int> & atomIndices, int numAtoms);

}; // class
} // namespace


#endif
		