// pb_FFTInteractionTypeCodes.h

#include <new>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <set>
#include "../fit/PositionUtils.h"
#include "../utils/AtomSpecifier.h"
#include "../gmath/Physics.h"
#include "../gcore/System.h"




#include "FFTInteractionTypeCodes.h"
using pb::FFTInteractionTypeCodes;

 FFTInteractionTypeCodes::FFTInteractionTypeCodes(){
    this->lsType=0;
    this->rfType=1;
    this->scType=2;
   
     
    }

/*
 int FFTInteractionTypeCodes::get_lsType(){return lsType;}
 int FFTInteractionTypeCodes::get_scType(){return scType;}
 int FFTInteractionTypeCodes::get_rfType(){return rfType;}
*/
