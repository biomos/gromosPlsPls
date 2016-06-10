// pb_FFTPoissonIterator.h

#ifndef INCLUDED_PB_FFTPoissonIterator
#define INCLUDED_PB_FFTPoissonIterator
#ifndef INCLUDED_PB_FFTGridType
#include "FFTGridType.h"
#endif
#ifndef INCLUDED_PB_DipoleDipole
#include "FFTDipoleDipole.h"
#endif
#ifndef INCLUDED_PB_FFTChargeDipole
#include "FFTChargeDipole.h"
#endif
#ifndef INCLUDED_PB_FFTBoundaryCondition
#include "FFTBoundaryCondition.h"
#endif
#ifndef INCLUDED_PB_PB_Parameters
#include "PB_Parameters.h"
#endif

namespace pb{



class FFTPoissonIterator{
	
	
	utils::AtomSpecifier atoms;
        utils::AtomSpecifier atoms_to_charge;

        pb::PB_Parameters ppp;
	pb::FFTBoundaryCondition bc;
	pb::FFTGridType gt;
	//pb::FFTDipoleDipole *ddTensor;
	//pb::FFTChargeDipole *cdTensor;
        

	double tinynum;
        double convergence;
        int maxsteps;
        double lambda;
        double epssolvent;
        double ls_boundary_eps;
        double rf_boundary_eps;
        double sc_boundary_eps;

        bool split_potentialbool;

        fftw_plan my_planV3_br;
        fftw_plan my_planV3_f;

        FFTDipoleDipole *ddTensor;
	FFTChargeDipole *cdTensor;



public:
    //constructor
        FFTPoissonIterator(utils::AtomSpecifier atoms, utils::AtomSpecifier atoms_to_charge, int maxsteps, double convergence, double lambda,
        FFTGridType gt, FFTBoundaryCondition bc, double epssolvent, bool split_potential);



    // deconstructor
    ~FFTPoissonIterator(){}

        //methods
       void setFFTPlans(
			int numFFTwThreads,
			int nx, int ny, int nz);


        int iterate_poisson(
			 std::vector<double> & Vx, std::vector<double> & Vy,  std::vector<double> &Vz,
			 std::vector<double> & Ex, std::vector<double> & Ey,  std::vector<double> & Ez,
			 std::vector<double> & inside,
			 std::vector<double> &  RPot);


        void postIteration(
			std::vector<double> & Vx, std::vector<double> & Vy, std::vector<double> & Vz,
			std::vector<double> & Ex, std::vector<double> & Ey, std::vector<double> &Ez,
			std::vector<double> & RPot,
			int nx, int ny, int nz,
			std::vector<double> & k_vecX,std::vector<double> & k_vecY,std::vector<double> & k_vecZ,
			int steps, bool converged);


        void updateVacuumField(
			std::vector<double> & Vx, std::vector<double> &Vy, std::vector<double> & Vz,
			std::vector<double> & Ex, std::vector<double> & Ey, std::vector<double> & Ez,
			std::vector<double> & inside,
			double deltaSigma);


        void realSpaceElectricField(
			std::vector<double> & Ex,
			std::vector<double> & Ey,
			std::vector<double> & Ez);

        void enforceZeroImaginaryComponent(std::vector<double> & complexVector);

        void reactionFieldHat(
			std::vector<double> & Ex,
			std::vector<double> & Ey,
			std::vector<double> & Ez,
			std::vector<double> & RPot,
			int nx,
			int ny,
			int nz,
			std::vector<double> & k_vecX,
			std::vector<double> & k_vecY,
			std::vector<double> & k_vecZ);
        
       void computeEfieldFromVacuumField(
			std::vector<double> & Ex,
			std::vector<double> & Ey,
			std::vector<double> & Ez,
			int nx,
			int ny,
			int nz,
			std::vector<double> & k_vecX,
			std::vector<double> & k_vecY,
			std::vector<double> & k_vecZ);

       	void fourierTransformedVacuumField(
			std::vector<double> & Vx,
			std::vector<double> & Vy,
			std::vector<double> & Vz,
			std::vector<double> & Ex,
			std::vector<double> & Ey,
			std::vector<double> & Ez);


        double computeResidualField(
			std::vector<double> & Ex, std::vector<double> & Ey, std::vector<double> & Ez,
			std::vector<double> & inside);


       double free_energy(std::vector<double> & pot);
       double free_energy_restricted(std::vector<double> & pot);


        void split_potential(
			std::vector<double> & Vx, std::vector<double> & Vy, std::vector<double> & Vz,
			std::vector<double> & Ex, std::vector<double> & Ey, std::vector<double> & Ez);


        void fft_grid(std::vector<double> & r_space, std::vector<double> & k_space);

        void bft_grid(std::vector<double> & rk_space);

}; // class
} // namespace


#endif
