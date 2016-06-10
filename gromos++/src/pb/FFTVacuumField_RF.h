// pb_FFTVacuumField_RF.h

#ifndef INCLUDED_PB_FFTVacuumField_RF
#define INCLUDED_PB_FFTVacuumField_RF
#ifndef INCLUDED_PB_FFTVacuumField
#include "FFTVacuumField.h"
#endif


namespace pb{



class FFTVacuumField_RF: virtual public FFTVacuumField{
	
	public:
        //constructor
        FFTVacuumField_RF(utils::AtomSpecifier atoms, FFTGridType gt, FFTBoundaryCondition bc);
       //deconstructor
        ~FFTVacuumField_RF(){}
       //methods


        void calcVacField(
			std::vector<double>& fldx,
			std::vector<double>& fldy,
			std::vector<double>& fldz);


        
        
	void calcVacFieldVincent(
			std::vector<double> & fldx,
			std::vector<double> & fldy,
			std::vector<double> & fldz,
			double epsRf);



}; // class
} // namespace


#endif