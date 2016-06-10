// pb_FFTChargeDipole.cc


#include <new>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <set>
#include "../fit/PositionUtils.h"
#include "../utils/AtomSpecifier.h"
#include "../gmath/Physics.h"
#include "../gcore/System.h"


#include "FFTChargeDipole.h"
#include "PB_Parameters.h"





using pb::FFTChargeDipole;

FFTChargeDipole::FFTChargeDipole(double epssolvent):ppp(epssolvent){

        this->epssolvent=epssolvent;
	this->tinynum = ppp.getTiny_real();

	
	/* k2 the squared norm of the wave vector */
	// is virtual double FFTChargeDipole::polarization(double k2){}
}





