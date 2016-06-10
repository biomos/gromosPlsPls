// pb_FDPoisson.cc


#include <new>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <set>
#include <fftw3.h>

#include "../fit/PositionUtils.h"
#include "../utils/AtomSpecifier.h"
#include "../gmath/Physics.h"
#include "../gcore/System.h"


#include "FFTPoisson.h"
#include "FFTPoissonIterator.h"
#include "FFTGridType.h"
#include "FFTBoundaryCondition.h"
#include "FFTDipoleDipole.h"
#include "PB_Parameters.h"
#include "FFTInteractionTypeCodes.h"
#include "FFTVacuumField.h"
#include "FFTVacuumField_LS.h"
#include "FFTVacuumField_RF.h"
#include "FFTInsideOutside.h"


using pb::FFTPoisson;
using pb::FFTVacuumField;
using pb::FFTInteractionTypeCodes;
using pb::FFTBoundaryCondition;
using pb::FFTInsideOutside;
using pb::FFTDipoleDipole;

FFTPoisson::FFTPoisson(utils::AtomSpecifier atoms,utils::AtomSpecifier atoms_to_charge, FFTGridType gt, FFTBoundaryCondition bc, int maxsteps, double convergence, double lambda,
        double epssolvent, bool split_potentialbool):
        ppp(epssolvent),pbiterator(atoms, atoms_to_charge,maxsteps, convergence, lambda,  gt, bc, epssolvent,split_potentialbool)
      {

        this->atoms=atoms;
        this->atoms_to_charge=atoms_to_charge;
    	this->tinynum=ppp.getTiny_real();
        this->convergence=convergence;
        this->maxsteps=maxsteps;
        this->lambda=lambda;

        this->gt=gt;
        this->bc=bc;

        this->epssolvent=epssolvent;
        this->split_potentialbool=split_potentialbool;



        std::cout << "# epssolvent " << epssolvent << endl;
        for (int i=0; i<atoms_to_charge.size();i++){
        std::cout << "# atom " << i << " charge " << atoms_to_charge.charge(i) << endl;
        std::cout << "# atom " << i << " radius " << atoms_to_charge.radius(i) << endl;
         }

		
	//pbiterator = new pb::FFTPoissonIterator(atoms, maxsteps, convergence, lambda,
        //gt, bc, epssolvent, split_potential);
		
	       // std::cout << "# FFTPoisson setup ... make fftplans " << endl;
		//setFFTPlans(1, gt.ngrdx, gt.ngrdy, gt.ngrdz);

             

                std::cout << "# FFTPoisson setup ... center atoms " << endl;
                center_atoms_on_grid(atoms, gt.centerx,gt.centery,gt.centerz);



	}
	

/* void FFTPoisson::setFFTPlans(
			int numFFTwThreads,
			int nx, int ny, int nz) {
		
		
				// use fftw version 3.x
				// 
				// we adopt the array layout as required by FFTW, i.e.
				// element (i) is REAL, wheras
				// element (i+1) is the corresponding COMPLEX to (i).
				// this makes the iteration sometimes a bit cumbersome...
     
				//std::vector<double>  test1;
                              //  std::vector<double>  test2;
                                
                               // test1.resize(gt.ngr3 * 2);
                               // test2.resize(gt.ngr3 * 2);

   

                               // Maria, use fftw_complex;
                    fftw_complex *test1, *test2;
                    test1=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*gt.ngr3);
                    test2=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*gt.ngr3);
                         //    fftw_complex test1[gt.ngr3];
                         //   fftw_complex test2[gt.ngr3];


				
			//	my_planV3_f = new FFTW3DPlanV3(
			//			nx, ny, nz,
			//			test1, test2,
        		//			FFTW3DPlanV3.FFTW_FORWARD,
			//			FFTW3DPlanV3.FFTW_MEASURE|FFTW3DPlanV3.FFTW_OUT_OF_PLACE);
				
			//	my_planV3_br = new FFTW3DPlanV3(
			//			nx, ny, nz,
			//			test1, test2,
			//      		FFTW3DPlanV3.FFTW_BACKWARD,
			//			FFTW3DPlanV3.FFTW_MEASURE|FFTW3DPlanV3.FFTW_OUT_OF_PLACE);
                                
                                fftw_plan my_planV3_f;
                                fftw_plan my_planV3_br;


                                my_planV3_f = fftw_plan_dft_3d(nx,ny,nz,test1,test2,FFTW_FORWARD,FFTW_MEASURE);
                                // Maria does it not have to be in place if we reuse the plan??
                               // my_planV3_br = fftw_plan_dft_3d(nx,ny,nz,test1,test2,FFT_BACKWARD,FFTW_MEASURE);
                                my_planV3_br = fftw_plan_dft_3d(nx,ny,nz,test1,test1,FFTW_BACKWARD,FFTW_MEASURE);

                                // free the space
                                fftw_free(test1);
                                fftw_free(test2);


 } */
	
	 void FFTPoisson::solve_poisson() {
		
		//check grid dimension....
		gridcheck();
	
		
	     
		int gridN[3];
                double gridD[3];
                gridN[0]=gt.ngrdx;
                gridN[1]=gt.ngrdy;
                gridN[2]=gt.ngrdz;
                gridD[0]=gt.drx;
                gridD[1]=gt.dry;
                gridD[2]=gt.drz;

                
                FFTInsideOutside inandout(gt.ngr3);
           
		std::vector<double>  inside;
                inside.resize(gt.ngr3);


           //     for (int iii =0; iii < inside.size() ;iii++) {
          //          inside[iii] = 0.0 ;
          //      }

                std::cout << "# FFTPoisson setup ...  insideoutside" << endl;
                inandout.inside_sol(gridN, gridD, gt.ncubes, atoms, inside);
		// force garbage collection to free surface points hash
		//System.gc();
		
		std::cout << "# FFTPoisson setup ... done insideoutside" << endl;

		//we adopt the array layout as required by FFTW, i.e.
		//element (i) is REAL, wheras
		//element (i+1) is the corresponding COMPLEX to (i).
		//this makes the iteration sometimes a bit cumbersome...
		std::cout << "# FFTPoisson setup ... grid" << endl;
                std::vector<double> Vx;
		std::vector<double> Vy;
		std::vector<double> Vz;

                Vx.resize(gt.ngr3 * 2);
                Vy.resize(gt.ngr3 * 2);
                Vz.resize(gt.ngr3 * 2);

                std::cout << "# FFTPoisson setup ... done gridresizing" << endl;
	     
		std::cout << "# FFTPoisson setup ... setup vacfield " << endl;
		setupVacuumField(inside, Vx, Vy, Vz);
		// force garbage collection to free vacuum field recycling stuff if necessary
		//System.gc();
		
		
		
		/* Now eliminate everything from the inside that is not 1.0 */
		inandout.integr_inside(inside);
		
		//init grids...
		//we adopt the array layout as required by FFTW, i.e.
		//element (i) is REAL, wheres 
		//element (i+1) is the corresponding COMPLEX to (i).
		//this makes the iteration sometimes a bit cumbersome...


                std::vector<double> Pot;
                std::vector<double> fldx;
                std::vector<double> fldy;
                std::vector<double> fldz;


                Pot.resize(gt.ngr3 * 2);
                fldx.resize(gt.ngr3 * 2);
                fldy.resize(gt.ngr3 * 2);
                fldz.resize(gt.ngr3 * 2);


		
		pbiterator.iterate_poisson(
				Vx, Vy, Vz,
				fldx, fldy, fldz, inside, 
				Pot);
		
	}

	/* set up initial (unmodified) vacuum field */
	  
	 /* inside non-zero if the corresponding grid cell is inside the solute
	 vx vacuum field - x component
	 vy vacuum field - y component
	 vz vacuum field - z component
	 */


         
	 void FFTPoisson::setupVacuumField(
			std::vector<double> & inside,
			std::vector<double> & vx, std::vector<double> & vy, std::vector<double> & vz) {

               FFTVacuumField *vacfield;//(atoms, gt, bc);
               FFTInteractionTypeCodes interx_codes;

               std::cout << "# setting up vacuum field ... " << endl;
               std::cout << "# bc.type is: " << bc.type << endl;

               try{
		if (bc.type == interx_codes.lsType) {
                    std::cout << "# need LS field " <<  endl;
			vacfield = new FFTVacuumField_LS(atoms,gt,bc);
		} else if (bc.type == interx_codes.rfType) {
                     std::cout << "# need RF field " <<  endl;
			vacfield = new FFTVacuumField_RF(atoms,gt,bc);
		}else if (bc.type == interx_codes.scType) {
                     std::cout << "# need SC field" <<  endl;
			vacfield = new FFTVacuumField_RF(atoms,gt,bc);
		} else{
                throw gromos::Exception("FFTPoisson","Boundary type not known ...");
                } //end of if
                } // end of try
                catch (const gromos::Exception &e){
                                         cerr << e.what() << endl;
                                         exit(1);
                }

		if (   bc.type !=  interx_codes.lsType )
                { 
                        vacfield->calcVacField(vx, vy, vz);
                }

                else {

		 // ewald-type vacuum field
                        vacfield->calcVacField(vx, vy, vz);

                        // we obtain the fourier transformed
                        // vacuum field in the ewald case, but we want
                        // the r-space vacuum field, so transform:
                        pbiterator.bft_grid(vx);
                        pbiterator.bft_grid(vy);
                        pbiterator.bft_grid(vz);

                }
			






		/* initial guess for the modified vacuum field -
                 * set V zero inside of solute */
                int size = gt.ngr3 * 2;
                int shift = 0;
                for (int index=0;index < size;index+=2) {
                        double heaviside = (1-inside[index - shift]);



                        vx[index] *= heaviside;
                        vy[index] *= heaviside;
                        vz[index] *= heaviside;
                        vx[index+1] *= heaviside;
                        vy[index+1] *= heaviside;
                        vz[index+1] *= heaviside;


 //  std::cout << "# @@@ vvv " <<  vx[index] << " " <<  vy[index] << " "  << vz[index] << endl;
  //                      std::cout << "# @@@ vvv " <<  vx[index+1] << " " <<  vy[index+1] << " "  << vz[index+1] << endl;


                        ++shift;
                }

                std::cout<< "# after setupfield: vx[0] " << vx[0] << endl;
                std::cout<< "# after setupfield: vy[0] " << vy[0] << endl;
                std::cout<< "# after setupfield: vz[0] " << vz[0] << endl;
                std::cout<< "# after setupfield: vx[100] " << vx[100] << endl;
                std::cout<< "# after setupfield: vy[100] " << vy[100] << endl;
                std::cout<< "# after setupfield: vz[100] " << vz[100] << endl;



               
        }
	
	void FFTPoisson::gridcheck() {
		
		int size = atoms.size();
		double rad;
		
		for (int i=0; i < size; ++i) {
			rad = atoms.radius(i);
                        try{
			if (   (  (atoms.pos(i))[0] + rad  ) > (gt.xlen)
					|| ((atoms.pos(i))[0] - rad) < (0.0)
					|| ((atoms.pos(i))[1] + rad) > (gt.ylen)
					|| ((atoms.pos(i))[1] - rad) < (0.0)
					|| ((atoms.pos(i))[2]+ rad) > (gt.zlen)
					|| ((atoms.pos(i))[2] - rad) < (0.0) ) {




                               std::cout << "# ATOM EXTENDING GRID!" << endl;
                               std::cout << "# ATOM ID:     " << (i+1)  << endl;
                               std::cout << "# ATOM POS:    " << (atoms.pos(i))[0] << " " << (atoms.pos(i))[1] << " " << (atoms.pos(i))[2]<< endl;
                               std::cout << "# ATOM RADIUS: "  <<  rad << endl;
                               std::cout << "# GRID START   "  <<  0.0 << " " << 0.0 << " " <<  0.0 << endl;
                               std::cout << "# GRID BOX     "  <<  (gt.xlen) << " " << (gt.ylen) << " " << (gt.zlen) << endl;


                       throw gromos::Exception("FDPoissonBolzmann","Atom extending grid. Exiting ...");

                        } // if endif
                        } // end of try

                         catch (const gromos::Exception &e){
                                         cerr << e.what() << endl;
                                         exit(1);
                                }

                } // loop over atoms

          }



        // shift atoms so that their cog coincides with the grid center
        void FFTPoisson::center_atoms_on_grid(utils::AtomSpecifier  & atoms,
                double gridcenterx, double gridcentery, double gridcenterz){

    double cogAtoms[3];
    // get cog
    cogAtoms[0]=fit::PositionUtils::cog(*atoms.sys(), atoms)[0];
    cogAtoms[1]=fit::PositionUtils::cog(*atoms.sys(), atoms)[1];
    cogAtoms[2]=fit::PositionUtils::cog(*atoms.sys(), atoms)[2];


    std::cout << "# initial cog x,y,z = " <<  cogAtoms[0] << " "<< cogAtoms[1] << " " << cogAtoms[2] << endl;

    // get the translation vector
    double tvec[3];
    tvec[0]=gridcenterx-cogAtoms[0];
    tvec[1]=gridcentery-cogAtoms[1];
    tvec[2]=gridcenterz-cogAtoms[2];


    // now shift


    for (int i=0;i<atoms.size(); i++){
        atoms.pos(i)[0] += tvec[0];
        atoms.pos(i)[1] += tvec[1];
        atoms.pos(i)[2] += tvec[2];

    }


    std::cout << "# gridcenter x,y,z = " << gridcenterx << " " << gridcentery << " " << gridcenterz << endl;

     // check the new cog
    cogAtoms[0]=fit::PositionUtils::cog(*atoms.sys(), atoms)[0];
    cogAtoms[1]=fit::PositionUtils::cog(*atoms.sys(), atoms)[1];
    cogAtoms[2]=fit::PositionUtils::cog(*atoms.sys(), atoms)[2];
    std::cout << "# final cog (after centering on grid) x,y,z = " <<  cogAtoms[0] << " "<< cogAtoms[1] << " " << cogAtoms[2] << endl;


        }
