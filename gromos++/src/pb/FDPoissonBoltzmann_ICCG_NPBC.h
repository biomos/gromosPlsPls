// pb_FDPoissonBoltzmann_ICCG_NPBC.h

#ifndef INCLUDED_PB_FDPoissonBoltzmann_ICCG_NPBC
#define INCLUDED_PB_FDPoissonBoltzmann_ICCG_NPBC

namespace pb{




class FDPoissonBoltzmann_ICCG_NPBC{

	
	
	int GPX;
	int GPY;
	int GPZ;
	int GPXGPY;
	int GPXGPYGPZ;
	
	int index;
	
public:
           //constructor
           FDPoissonBoltzmann_ICCG_NPBC(int GPX, int GPY, int GPZ);
	   // deconstructor
           ~FDPoissonBoltzmann_ICCG_NPBC(){}

           //methods
	
	
            void initpciccg(std::vector<double> &ldiag,
            std::vector<double> &epsCgrid,
            std::vector<double> &epsIgrid,
            std::vector<double> &epsJgrid,
            std::vector<double> &epsKgrid);

	
	
	
	    void pciccg(std::vector<double> &zvec, std::vector<double> &ldiag,
            std::vector<double> &rhogrid,
	    std::vector<double> &epsIgrid,
            std::vector<double> &epsJgrid,
            std::vector<double> &epsKgrid);

	
            void gqact(std::vector<double> &zvec, std::vector<double> &pvec,
            std::vector<double> &epsCgrid,
            std::vector<double> &epsIgrid,
            std::vector<double> &epsJgrid,
            std::vector<double> &epsKgrid);




}; //class
} //namespace
#endif