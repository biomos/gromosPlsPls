// pb_FDPoissonBoltzmann.h

#ifndef INCLUDED_PB_FDPoissonBoltzmann
#define INCLUDED_PB_FDPoissonBoltzmann

#ifndef INCLUDED_PB_FDPoissonBoltzmann_ICCG_PBC
#include "FDPoissonBoltzmann_ICCG_PBC.h"
#endif
#ifndef INCLUDED_PB_FDPoissonBoltzmann_ICCG_NPBC
#include "FDPoissonBoltzmann_ICCG_NPBC.h"
#endif
#ifndef INCLUDED_PB_PB_Parameters
#include "PB_Parameters.h"
#endif

namespace pb{


    
class FDPoissonBoltzmann{
        utils::AtomSpecifier atoms;
      //  gcore::System sys;
        utils::AtomSpecifier atoms_to_charge;

	// vector for radii to store and manipulate
	std::vector<double> radii;
        
        // whether we are in pbc or npbc
        bool pbc;

        // number of grid points along X
         int GPX;

        
      // number of grid points along Y
         int GPY;

       // number of grid points along Z
         int GPZ;


        // multiples
         int GPXGPY;
         int GPXGPYGPZ;

         // density and potential grids
         std::vector<double> rhogrid;
         std::vector<double> phigrid; // (will be the solution when solving for the electrostatic potential)

         // epsIgrid[GPXGPYGPZ] stores the permittivity surrounding each face
         // of a gridpoint[I][J][K] shifted in direction I
          std::vector<double> epsIgrid;

         // epsJgrid[GPXGPYGPZ] stores the permittivity surrounding each face
         // of a gridpoint[I][J][K] shifted in direction J
          std::vector<double> epsJgrid;

         // epsKgrid[GPXGPYGPZ] stores the permittivity surrounding each face
         // of a gridpoint[I][J][K] shifted in direction K
          std::vector<double> epsKgrid;

       // epsCgrid[GPXGPYGPZ] stores the sum of the permittivities at the six faces
       // (diagonal of coefficient matrix A)
          std::vector<double> epsCgrid;




	  // gridspacing
         double gridspacing;

	 // epsilon of solute and solvent
         double epssolute;
         double epssolvent;

	 // grid origin
         double gridstart[3];
	 // grid center
	 double gridcenter[3];

         pb::PB_Parameters ppp;
        

         
       // pb::FDPoissonBoltzmann_ChargeGrid PBchgGrid;
    //    pb::FDPoissonBoltzmann_EpsGrid PBeps;

      //  pb::FDPoissonBoltzmann_ICCG iccg;
  




 public:
  // constructor

     FDPoissonBoltzmann(utils::AtomSpecifier atoms,utils::AtomSpecifier atoms_to_charge,
             int gridpointsX, int gridpointsY, int gridpointsZ,
             double gridspace, bool pbc,
			double epssolvent, ofstream &os);
 
  
   // deconstructor
  ~FDPoissonBoltzmann(){}


  // methods

  void setupGrid(bool newphi, ofstream &os, double gridstartX=0, double gridstartY=0, double gridstartZ=0, double gridcenterX=0, double gridcenterY=0, double gridcenterZ=0);
  //return functions
  void getgridcenter(double& X, double& Y, double& Z);
  void getgridstart(double& X, double& Y, double& Z);
  //void getgridcenter(double gc[]); //double& X, double& Y, double& Z);
  //void getgridstart(double gs[]); //double& X, double& Y, double& Z);
  bool solveforpotential_pbc(int maxits, double acceptance,FDPoissonBoltzmann_ICCG_PBC iccg, ofstream &os);
  bool solveforpotential_npbc(int maxits, double acceptance,FDPoissonBoltzmann_ICCG_NPBC iccg, ofstream &os);
  double dGelec(ofstream &os, vector<double> *potentials=NULL);
  double getdG();
  double getdG_restricted(ofstream &os, vector<double> *potentials=NULL);
  void increasebox(ofstream &os);
  void increasegrid(ofstream &os);
  void atomshift(ofstream &os);
  void gridcheck(ofstream &os);
  int index(int x, int y, int z);


  
 // void setpermittivity();
  void setboundarySolvent();


   void radiusboundaryEPSI( std::vector<double>& epsgrid) ;
   void radiusboundaryEPSJ( std::vector<double>& epsgrid) ;
   void radiusboundaryEPSK( std::vector<double>& epsgrid) ;


    void chargeGridtrilinear(ofstream &os);

    void DebyeHueckel(double gridspacing, double gridstart[3],
		      double kappa,  std::vector<double> & rhogrid, ofstream &os);

 /* double getpotential();
  double getgridstart();*/


}; // class
} // namespace


#endif
