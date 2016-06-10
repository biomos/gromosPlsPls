// pb_FFTChargeShapingFunction.cc


#include <new>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <set>
#include "../fit/PositionUtils.h"
#include "../utils/AtomSpecifier.h"
#include "../gmath/Physics.h"
#include "../gcore/System.h"


#include "FFTChargeShapingFunction.h"
#include "PB_Parameters.h"
#include "FFTChargeShapingFunction_Hat.h"
#include "FFTChargeShapingFunction_Parabola.h"



using pb::FFTChargeShapingFunction;

FFTChargeShapingFunction::FFTChargeShapingFunction(){

	this->tinynum = ppp.getTiny_real();
      
        this->HatType=0;
        this->ParabolaType=1;
    //    this->OtherParabolaType=2;

 

}

     /*  FFTChargeShapingFunction FFTChargeShapingFunction::getCSF(int type){
       try{
       if (type == HatType){return FFTChargeShapingFunction_Hat();}
       else if (type == ParabolaType){
           return FFTChargeShapingFunction_Parabola();
       }
       else{
         throw gromos::Exception("FFTChargeShapingFunction","Invalid charge shaping function code. Exiting.");
       }
       }// end of try
        catch (const gromos::Exception &e){
                                         cerr << e.what() << endl;
                                         exit(1);
                                }
       }
	*/

	
	