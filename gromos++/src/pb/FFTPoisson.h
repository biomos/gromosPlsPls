// pb_FFTPoisson.h

#ifndef INCLUDED_PB_FFTPoisson
#define INCLUDED_PB_FFTPoisson
#ifndef INCLUDED_PB_PB_Parameters
#include "PB_Parameters.h"
#endif
#ifndef INCLUDED_PB_FFTPoissonIterator
#include "FFTPoissonIterator.h"
#endif
#ifndef INCLUDED_PB_FFTBoundaryCondition
#include "FFTBoundaryCondition.h"
#endif
#ifndef INCLUDED_PB_FFTGridType
#include "FFTGridType.h"
#endif
#ifndef INCLUDED_PB_FFTInteractionTypeCodes
#include "FFTInteractionTypeCodes.h"
#endif
#ifndef INCLUDED_PB_FFTVacuumField
#include "FFTVacuumField.h"
#endif


namespace pb{



class FFTPoisson{
	
	utils::AtomSpecifier atoms;
        utils::AtomSpecifier atoms_to_charge;
        
        pb::PB_Parameters ppp;
	pb::FFTPoissonIterator pbiterator;
	pb::FFTBoundaryCondition bc;
        pb::FFTGridType gt;
	
	double tinynum;
        double convergence;
        int maxsteps;
        double lambda;
        double epssolvent;
        bool split_potentialbool;
       
	//static j3DFFT j3DFFT;
	
	
	
	//plans for FFTW(V3)
	//fftw_plan my_planV3_f; //forward
	//fftw_plan my_planV3_br; //backward
public:
    // constructor
        FFTPoisson(utils::AtomSpecifier atoms,utils::AtomSpecifier atoms_to_charge, FFTGridType gt, FFTBoundaryCondition bc, int maxsteps, double convergence, double lambda,
                double epssolvent, bool split_potentialbool);



    // deconstructor
    ~FFTPoisson(){}

    // methods
    // void setFFTPlans(
    //		int numFFTwThreads,
    //			int nx, int ny, int nz);

	
     void solve_poisson();
		
     void setupVacuumField(
			std::vector<double> &  inside,
			std::vector<double> & vx,std::vector<double> &  vy, std::vector<double> &  vz);

      void center_atoms_on_grid(utils::AtomSpecifier  & atoms, double gridcenterx, double gridcentery, double gridcenterz);

      void gridcheck();

            }; // class
} // namespace


#endif
