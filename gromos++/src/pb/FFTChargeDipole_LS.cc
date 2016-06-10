// pb_FFTChargeDipole_LS.cc


#include <new>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <set>
#include "../fit/PositionUtils.h"
#include "../utils/AtomSpecifier.h"
#include "../gmath/Physics.h"
#include "../gcore/System.h"


#include "FFTChargeDipole_LS.h"
#include "FFTChargeDipole.h"



using pb::FFTChargeDipole_LS;
using pb::FFTChargeDipole;

FFTChargeDipole_LS::FFTChargeDipole_LS(double epssolvent) : FFTChargeDipole(epssolvent) {
}

double FFTChargeDipole_LS::polarization(double k2) {
		return 1.0;
	}






