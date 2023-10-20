/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

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
	FFTVacuumField(utils::AtomSpecifier atoms, FFTGridType gt, FFTBoundaryCondition bc, ofstream &os);
    //deconstructor
     virtual  ~FFTVacuumField(){}
    //methods
	
   /* FFTVacuumField getVF(
			int type,
			utils::AtomSpecifier atoms,
			FFTGridType gt, FFTBoundaryCondition bc);
*/
	virtual void calcVacField(
				  std::vector <double> &  Vx, std::vector <double> & Vy, std::vector <double> &  Vz, ofstream &os){}



}; // class
} // namespace


#endif

