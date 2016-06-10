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

         pb::PB_Parameters ppp;
        

         
       // pb::FDPoissonBoltzmann_ChargeGrid PBchgGrid;
    //    pb::FDPoissonBoltzmann_EpsGrid PBeps;

      //  pb::FDPoissonBoltzmann_ICCG iccg;
  




 public:
  // constructor

     FDPoissonBoltzmann(utils::AtomSpecifier atoms,utils::AtomSpecifier atoms_to_charge,
             int gridpointsX, int gridpointsY, int gridpointsZ,
             double gridspace, bool pbc,
              double epssolvent);
 
  
   // deconstructor
  ~FDPoissonBoltzmann(){}


  // methods

  void setupGrid(bool newphi);
  bool solveforpotential_pbc(int maxits, double acceptance,FDPoissonBoltzmann_ICCG_PBC iccg);
  bool solveforpotential_npbc(int maxits, double acceptance,FDPoissonBoltzmann_ICCG_NPBC iccg);
  double dGelec();
  double getdG();
  double getdG_restricted();
  void gridcheck();
  int index(int x, int y, int z);
 // void setpermittivity();
  void setboundarySolvent();


   void radiusboundaryEPSI( std::vector<double>& epsgrid) ;
   void radiusboundaryEPSJ( std::vector<double>& epsgrid) ;
   void radiusboundaryEPSK( std::vector<double>& epsgrid) ;


    void chargeGridtrilinear();

    void DebyeHueckel(double gridspacing, double gridstart[3],
              double kappa,  std::vector<double> & rhogrid);

 /* double getpotential();
  double getgridstart();*/


}; // class
} // namespace


#endif