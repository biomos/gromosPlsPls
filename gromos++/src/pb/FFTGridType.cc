// pb_FFTGridType.cc


#include <new>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <set>
#include "../fit/PositionUtils.h"
#include "../utils/AtomSpecifier.h"
#include "../gmath/Physics.h"
#include "../gcore/System.h"


#include "FFTGridType.h"
#include "PB_Parameters.h"



using pb::FFTGridType;

FFTGridType::FFTGridType(int ngrdx, int ngrdy, int ngrdz, double xlen, double ylen, double zlen, int ncubes) {

     // grid point numbers
         this->ngrdx = ngrdx;
	 this->ngrdy = ngrdy;
	 this->ngrdz = ngrdz;
	 this->ngrdxy = ngrdx * ngrdy;
	 this->ngr3 = ngrdx * ngrdy * ngrdz;
      // box edges
	 this->xlen = xlen;
	 this->ylen = ylen;
	 this->zlen = zlen;
      // grid center
         this->centerx=0.5*xlen;
         this->centery=0.5*ylen;
         this->centerz=0.5*zlen;


	 double pi2 = 2 * ppp.getPI();

         // volume of box
	 this->vol=xlen*ylen*zlen;

         /* spacing between r-space grid points */
	 this->drx=xlen/ngrdx;
	 this->dry=ylen/ngrdy;
	 this->drz=zlen/ngrdz;
	 /* spacing between k-space grid points */
	 this->dkx=pi2/xlen;
	 this->dky=pi2/ylen;
	 this->dkz=pi2/zlen;

 	 this->ncubes = ncubes;

	}
void FFTGridType::dumpparameters () {
	std::cout << "# GRID PARAMETERS" << endl;
	std::cout << "# ---------------" << endl;
        std::cout << "# DIM (GPX GPY GPZ): " << ngrdx << " " << ngrdy << " " << ngrdz << endl;
        std::cout << "# CELL: " << xlen << " " << ylen << " " << zlen << endl;
        std::cout << "# SPACE(R): " << drx << " " << dry << " " << drz << endl;
        std::cout << "# SPACE(K): " << dkx << " " << dky << " " << dkz << endl;
        std::cout << "# CUBES: " << ncubes << endl;
        std::cout << "# VOLUME: " << vol << " nm^3" << endl;

	}