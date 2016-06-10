// pb_FFTDipoleDipole_LS.h

#ifndef INCLUDED_PB_FFTDipoleDipole_LS
#define INCLUDED_PB_FFTDipoleDipole_LS
#ifndef INCLUDED_PB_FFTDipoleDipole
#include "FFTDipoleDipole.h"
#endif


namespace pb{



class FFTDipoleDipole_LS: public FFTDipoleDipole{

	
	double k0EpsilonFactor;
	double es_1_es;
        double epsB;
        
public:
        //constructor

	FFTDipoleDipole_LS(double epsB,double epsS);

        //deconstructor

        ~FFTDipoleDipole_LS(){};


        // methods
        
	void updateTensorK0(double (& tensor)[3][3]);
	
	double computeAFactor(double k2);
	
	double computeBFactor(double k2);




}; // class
} // namespace


#endif
