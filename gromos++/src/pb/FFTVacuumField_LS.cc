// pb_FFTVacuumField_LS.cc


#include <new>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <set>
#include "../fit/PositionUtils.h"
#include "../utils/AtomSpecifier.h"
#include "../gmath/Physics.h"
#include "../gcore/System.h"


#include "FFTVacuumField_LS.h"
#include "FFTVacuumField.h"
#include "FFTChargeShapingFunction.h"
#include "FFTChargeShapingFunction_Hat.h"
#include "FFTChargeShapingFunction_Parabola.h"


using pb::FFTVacuumField_LS;
using pb::FFTVacuumField;
using pb::FFTChargeShapingFunction;


 FFTVacuumField_LS::FFTVacuumField_LS(utils::AtomSpecifier atoms, FFTGridType gt, FFTBoundaryCondition bc):FFTVacuumField(atoms, gt, bc){
//   FFTVacuumField_LS::FFTVacuumField_LS(){
	
	this->ion_count1 = 0;
	this->ion_count2 = 0;
	this->ngrdx = gt.ngrdx;
	this->ngrdy = gt.ngrdy;
	this->ngrdz = gt.ngrdz;
	this->dkx = gt.dkx;
	this->dky = gt.dky;
	this->dkz = gt.dkz;

        this->cstype=ppp.get_default_chshape();

        csfunc=new FFTChargeShapingFunction();
      //  int tmp1 = csfunc->HatType;
      //  std::cout << "# TMP1 " << tmp1 << endl;

        std::cout << "# FFTVacuumField : cstype " << cstype << endl;
        std::cout << "# FFTVacuumField : csfunc->HatType " << (csfunc->HatType) << endl;
        std::cout << "# FFTVacuumField : csfunc->ParabolaType " << (csfunc->ParabolaType) << endl;

      try{
       if ( cstype ==  csfunc->HatType ){
          csfunc = new FFTChargeShapingFunction_Hat();
       }
       else if ( cstype == csfunc->ParabolaType){
           csfunc = new FFTChargeShapingFunction_Parabola();
       }
       else{
         throw gromos::Exception("FFTVacuumField_LS","Invalid charge shaping function code. Exiting.");
       }
       }// end of try
        catch (const gromos::Exception &e){
                                         cerr << e.what() << endl;
                                         exit(1);
                                }
       

        

		
		
		ion_list1.resize(atoms.size());
		ion_list2.resize(atoms.size());
		
	

		Vx1_store.resize(gt.ngr3);
		Vx2_store.resize(gt.ngr3);
		Vy1_store.resize(gt.ngr3);
                Vy2_store.resize(gt.ngr3);
		Vz1_store.resize(gt.ngr3);
		Vz2_store.resize(gt.ngr3);

		

		//create list with ions
		makeatomlist();





                


}
	
        void FFTVacuumField_LS::complexFromDouble(
			std::vector <double> & doubleArray, std::vector <double> & complexArray) {

		for (int ii = 0; ii < doubleArray.size(); ii++) {
			complexArray[2 * ii] =      doubleArray[ii];
			complexArray[2 * ii + 1] =  doubleArray[ii] ; // 0;  // Maria, why not 0??
		}
		}

	
	
	void FFTVacuumField_LS::makeatomlist() {
		
		
		/* determine the charged atoms and select the appropriate cs */

            try{
		ion_count1 = 0;ion_count2 = 0;
		for (int ion=0;ion< atoms.size();ion++) {
			if ((atoms.radius(ion) > (bc.alpha1+tinynum))){


					if (fabs(atoms.charge(ion))>tinynum) {
				ion_list1[ion_count1]=ion;
				ion_count1++;
			}
                        }
                        else{
                    
		//	if ((atoms.radius(ion) < (bc.alpha1-tinynum)))
				if	(fabs(atoms.charge(ion))>tinynum) {
				
                                if (atoms.radius(ion) < bc.alpha2){
                                std::cout << "Ion " << ion << " : radius " << atoms.radius(ion) << endl;
                                throw gromos::Exception("FFTVacuumField_LS","Ion too small. Exiting.");
                                }

				ion_list2[ion_count2]=ion;
				ion_count2++;
			}    
		}

                 } // end of for-loop

                 }// end of try
                                         catch (const gromos::Exception &e){
                                         cerr << e.what() << endl;
                                         exit(1);
                                         }

                std::cout << "# IONLIST 1: NUMBER OF ENTRIES: " << ion_count1 <<endl;
                for (int ion=0;ion<ion_count1;ion++){
                    std::cout << "# list 1 : entry " << ion << " has ion_list1[ "<< ion << " ] = " << ion_list1[ion] <<  " and radius " << atoms.radius(ion_list1[ion]) << endl;
                }
                std::cout << "# IONLIST 2: NUMBER OF ENTRIES: " << ion_count2 <<endl;
                for (int ion=0;ion<ion_count2;ion++){
                    std::cout << "# list 2 : entry " << ion << " has ion_list2[ "<< ion << " ] = " << ion_list2[ion] <<  " and radius " << atoms.radius(ion_list2[ion]) << endl;
                }

               

                

	}

	
	
	void FFTVacuumField_LS::calcVacField(
			std::vector <double> & fldx_k,
                        std::vector <double> & fldy_k,
                        std::vector <double> & fldz_k) {
		
	std::cout<< "# FFTVacuumField_LS::calcVacField :  Calculating vac field from scratch ..."  << endl;
		
		double kax = pi2/gt.drx;
		double kay = pi2/gt.dry;
		double kaz = pi2/gt.drz;
		
		if (ion_count1 > 0) {
		
		std::cout<< "# Computing position-independent part of the vacuum field due to big ions, with alpha " << bc.alpha1 << endl;
	
			positionIndependentVacField(
					Vx1_store, Vy1_store, Vz1_store,
					kax, kay, kaz, 
					bc.nalias1, bc.alpha1);
		}


                std::cout<< "# after positionIndependentVacField: Vx1_store[0] " << Vx1_store[0] << endl;
                std::cout<< "# after positionIndependentVacField: Vy1_store[0] " << Vy1_store[0] << endl;
                std::cout<< "# after positionIndependentVacField: Vz1_store[0] " << Vz1_store[0] << endl;
                std::cout<< "# after positionIndependentVacField: Vx1_store[100] " << Vx1_store[100] << endl;
                std::cout<< "# after positionIndependentVacField: Vy1_store[100] " << Vy1_store[100] << endl;
                std::cout<< "# after positionIndependentVacField: Vz1_store[100] " << Vz1_store[100] << endl;



		
		if (ion_count2 > 0) {
			std::cout<< "# Computing position-independent part of the vacuum field due to small ions, with alpha " << bc.alpha2 << endl;

                        positionIndependentVacField(
					Vx2_store, Vy2_store, Vz2_store,
					kax, kay, kaz, 
					bc.nalias2, bc.alpha2);
		}
		
		recyclefield(fldx_k, fldy_k, fldz_k);		
		

                std::cout<< "# after recyclefield: fldx_k[0] " << fldx_k[0] << endl;
                std::cout<< "# after recyclefield: fldy_k[0] " << fldy_k[0] << endl;
                std::cout<< "# after recyclefield: fldz_k[0] " << fldz_k[0] << endl;
                std::cout<< "# after recyclefield: fldx_k[100] " << fldx_k[100] << endl;
                std::cout<< "# after recyclefield: fldy_k[100] " << fldy_k[100] << endl;
                std::cout<< "# after recyclefield: fldz_k[100] " << fldz_k[100] << endl;




	}


	void FFTVacuumField_LS::positionIndependentVacField(
			std::vector <double> &  destX, std::vector <double> &  destY, std::vector <double> &  destZ,
			double kax, double kay, double kaz, 
			int nAlias, double alpha) {
		
		double fld[3];

                

		for (int i=0;i< ngrdx;i++) {

                    std::cout << "# FFTVacuumField_LS::positionIndependentVacField : slice " << i << endl;



			double kx = dkx*(((i+ngrdx/2)%ngrdx)-ngrdx/2);
			for (int j=0;j<ngrdy;j++) {
				double ky = dky*(((j+ngrdy/2)%ngrdy)-ngrdy/2);
				for (int k=0;k<ngrdz;k++) {
					double kz = dkz*(((k+ngrdz/2)%ngrdz)-ngrdz/2);
					
					/* k-space field using aliases */
					fld[0] = 0.0;
                                        fld[1] = 0.0;
                                        fld[2] = 0.0;
                                        
					for (int iax=-nAlias; iax<=nAlias; iax++) {
						for (int iay=-nAlias; iay<=nAlias; iay++) {
							for (int iaz=-nAlias; iaz<=nAlias; iaz++) {									
								double psi_k_t = csfunc->calc(kx+iax*kax,
                                                                        ky+iay*kay,
                                                                        kz+iaz*kaz,
                                                                        alpha,
                                                                        eps0);



 //std::cout << "# FFTVacuumField_LS::positionIndependentVacField :  psi_k_t " << psi_k_t << endl;

								fld[0] += - (kx+iax*kax) * psi_k_t;
								fld[1] += - (ky+iay*kay) * psi_k_t;
								fld[2] += - (kz+iaz*kaz) * psi_k_t;
							}
						}
					}
					
					//real component at (i), imaginary component at (i+1)
					//final int index = 2 * (k + ngrdz * ( j + ngrdy * i ));
					
					// actually, they're the same, so we can save quite a bit
					// here...
					int index = k + ngrdz * ( j + ngrdy * i );
					destX[index] 	= fld[0];
					//destX[index+1] 	= fld[0];
					destY[index] 	= fld[1];
					//destY[index+1] 	= fld[1];
					destZ[index] 	= fld[2];
					//destZ[index+1] 	= fld[2];
				}
			}
		}
	}
	


	
	
	void FFTVacuumField_LS::recyclefield(
			std::vector <double> &  destX, std::vector <double> & destY, std::vector <double> &  destZ	) {
		
            int tmpsizex=destX.size();
            destX.resize(tmpsizex,0.0);
            int tmpsizey=destY.size();
            destY.resize(tmpsizey,0.0);
            int tmpsizez=destZ.size();
            destZ.resize(tmpsizez,0.0);

            std::vector<double> multipliers;
            std::vector<double> complexTmp;
            multipliers.resize(2 * ngrdx * ngrdy * ngrdz);
            complexTmp.resize(gt.ngr3 * 2);

		std::cout << "# FFTVacuumField_LS::recyclefield: Recycling vacuum field" << endl;
		
		updateMultipliers(multipliers,
				ion_list1, ion_count1);		
		
                complexFromDouble(Vx1_store, complexTmp);
		for (int i=0;i<destX.size();i++){
                    destX[i]+= complexTmp[i]*multipliers[i];
                    }
                complexFromDouble(Vy1_store, complexTmp);
                for (int i=0;i<destY.size();i++){
                    destY[i]+=complexTmp[i]*multipliers[i];
                    }
                complexFromDouble(Vz1_store, complexTmp);
                for (int i=0;i<destZ.size();i++){
                    destZ[i]+=complexTmp[i]*multipliers[i];
                    }
		
		updateMultipliers(multipliers,
				ion_list2, ion_count2);		

                complexFromDouble(Vx2_store, complexTmp);
                for (int i=0;i<destX.size();i++){
                   destX[i]+= complexTmp[i]*multipliers[i];
                    }
                complexFromDouble(Vy2_store, complexTmp);
                for (int i=0;i<destY.size();i++){
                   destY[i]+=complexTmp[i]*multipliers[i];
                    }
                complexFromDouble(Vz2_store, complexTmp);
                for (int i=0;i<destZ.size();i++){
                   destZ[i]+=complexTmp[i]*multipliers[i];
                    }

	}
	
	
	void FFTVacuumField_LS::updateMultipliers(
			std::vector <double> &  multipliers,
			std::vector <int> & atomIndices, int numAtoms) {
				
		int index;
		
		double kx,ky,kz;
		
		int shift = 0;
		for (int i=0;i< ngrdx;i++) {
			kx = dkx*(((i+ngrdx/2)%ngrdx)-ngrdx/2);
			for (int j=0;j<ngrdy;j++) {
				ky = dky*(((j+ngrdy/2)%ngrdy)-ngrdy/2);
				for (int k=0;k<ngrdz;k++) {
					kz = dkz*(((k+ngrdz/2)%ngrdz)-ngrdz/2);
					index = k + ngrdz * ( j + ngrdy * i );
					index += shift;
					double sumIm = 0.0, sumRe = 0.0;
					for (int count=0;count< numAtoms;count++) {
						int ion = atomIndices[count];
						double kr = kx * (atoms.pos(ion))[0] + ky * (atoms.pos(ion))[1] + kz * (atoms.pos(ion))[2];
						double charge = atoms.charge(ion);
						sumIm += charge * cos(kr);
						sumRe += charge * sin(kr);					
					}
					multipliers[index] = sumRe;
					multipliers[index+1] = sumIm;
					
					++shift;
				}
			}
		}	
	}	

