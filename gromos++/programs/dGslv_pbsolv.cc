/**
 * @file dGslv_pbsolv.cc
 * calculates the correction to the
 * electrostatic component of the
 * solvation free energy for either
 * going from LS/PBC to CB/NPBC or for
 * going from RF/PBC to CB/NPBC
 */

/**
 * @page programs Program Documentation
 *
 * @anchor dGslv_pbsolv
 * @section dGslv_pbsolv calculates a correction to the electrostatic component of the solvation free energy
 * @author @ref Maria Reif
 * @date 25-10-10
 *
 * Progam dGslv_pbsolv will compute two electrostatic components of the solvation free energy,
 * namely one for CB/NPBC and one for the user-specified electrostatics scheme (LS or RF).
 *
 * In the following, the abbreviation dGslv will be used to denote an electrostatic component of
 * the solvation free energy.
 *
 * dGslv will be computed for a user-specified group of atoms.
 * Note that all atoms of the solute topology block will be nonpolarizable,
 * i.e. they will be assigned a relative dielectric permittivity of one.
 *
 *
 * The CB/NPBC dGslv will be computed from a FD algorithm.
 *
 * 
 * The LS/PBC dGslv will be computed from both a FD and a FFT algorithm (for comparison).
 * But the user should use the FD value to cancel possible grid discretization errors.
 * That is, the resulting correction is: CB/NPBC[FD] - LS/PBC[FD]
 *
 *
 * The RF/PBC dGslv will be computed from a FFT algorithm.
 * But the user should use the corresponding LS calculation to compute a correction
 * to cancel possible grid discretization errors.
 * That is, the resulting correction is: CB/NPBC[FD] - LS/PBC[FD] + LS/PBC[FFT] - RF/PBC[FFT]
 *
 *
 * In the LS-scheme, tinfoil boundary conditions are used and a hat charge shaping function
 * will be used.
 *
 * In the RF-scheme, a user-specified relative dielectric permittivity is used.
 * Note that a relative dielectric permittivity of one implies no application of a reaction-fied correction.
 *
 *
 *
 * The solute will be centered in the computational box, with its center of geometry.
 *
 * The algorithms employed are from these papers:
 * FD: Comput. Phys. Commun. 62, 187-197 (1991)
 * FFT: J. Chem. Phys. 116, 7434-7451 (2002), J. Chem. Phys. 119, 12205-12223 (2003)
 *
 *
 *  *
 * Example:
 * @verbatim
 dGslv_pbsolv 
 @topo topo.top
 @atoms 1:a
 @atoms_to_charge 1:a
 @frame_g96 coor.dat
 @scheme RF
 @eps  66.6
 @epsRF 66.6
 @epsNPBC 78.4
 @rcut 1.4
 @gridspacing 0.02
 @gridpointsXYZ 158 158 158
 @maxiter 600
 @FFT_Cubes 4
 @probeIAC 5
 @probeRAD 0.14
 @HRAD 0.05
 @pbc r gbond
 @radscal 1.0
 @rminORsigma 0
 @endverbatim
 *
 *
 *
 */

#include <cassert>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <map>
#include <cmath>


#include <fftw3.h>


#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/AtomicRadii.h"
#include "../src/fit/PositionUtils.h"

#include "../src/pb/FDPoissonBoltzmann.h"
#include "../src/pb/FDPoissonBoltzmann_ICCG_NPBC.h"
#include "../src/pb/FDPoissonBoltzmann_ICCG_PBC.h"
#include "../src/pb/PB_Parameters.h"
#include "../src/pb/FFTPoisson.h"
#include "../src/pb/FFTGridType.h"
#include "../src/pb/FFTBoundaryCondition.h"



using namespace std;
using namespace gcore;
using namespace gio;
using namespace gmath;
using namespace bound;
using namespace args;
using namespace utils;
using namespace pb;

int main(int argc, char **argv){

  Argument_List knowns;
  knowns << "topo" 
         << "pbc"
         << "atoms" <<  "atoms_to_charge" << "frame_g96" << "scheme" << "eps"
         << "epsRF" << "rcut"
         << "gridspacing" << "gridpointsXYZ" << "maxiter" // << "convergence"
         << "FFT_Cubes" << "probeIAC" << "probeRAD" << "HRAD" <<  "epsNPBC" <<  "radscal" << "rminORsigma"; // << "NWAT" << "boxc1"; // << "sphereRAD" << "sphereatoms";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo    <molecular topology file>\n";
  usage += "\t@pbc     <boundary type> [<gather method>]\n";
  usage += "\t@atoms   <atoms to include>\n";
  usage += "\t@atoms_to_charge   <atoms to charge>\n";
  usage += "\t@frame_g96    <g96 coordinates>\n";
  usage += "\t@scheme   <electrostatics scheme: LS or RF>\n";
  usage += "\t@eps    <solvent relative dielectric permittivity >\n";
  usage += "\t@epsRF    <reaction field relative dielectric permittivity >\n";
  usage += "\t@rcut    <cutoff distance (ONLY USED IF scheme==RF; default 1.4)>\n";
  usage += "\t@gridspacing    <grid spacing>\n";
  usage += "\t@gridpointsXYZ   <number of gridpoints along the X,Y,Z direction>\n";
  usage += "\t@maxiter   <maximum number of iteration steps>\n";
  //usage += "\t@convergence   <convergence threshold criterion (kJ/mol)>\n";
  usage += "\t@FFT_Cubes   <number of cubes in FFT for boundary smoothing>\n";
  usage += "\t@probeIAC   <integer atom code to take for radius calculation (for water, it would be 4 or 5 depending on the ff)>\n";
  usage += "\t@probeRAD   <probe radius (for water, it would be 0.14 [nm]; default 0.0 nm)>\n";
  usage += "\t@HRAD   <your desired hydrogen radius; default 0.05 nm>\n";
 // usage += "\t@NWAT   <number of water molecules (for C1-correction estimate)>\n";
  usage += "\t@epsNPBC   <relative dielectric permittivity for NPBC calculation (reasonable: 78.4)>\n";
  //usage += "\t@chargefac <modify charges given in topology for atoms to be charged by this multiplicative factor>   >\n";
  usage += "\t@radscal <scale non-H radii with this factor (in case you want to play with radii; default 1.0)>\n";
 // usage += "\t@boxc1 <box edge for c1 calculation>\n";
  usage += "\t@rminORsigma <which radii to use: rmin (0) or sigma (1); default 0>\n";
  
 try{
  Arguments args(argc, argv, knowns, usage);
  

  // read topology
  args.check("topo",1);
  gio::InTopology it(args["topo"]);
  System sys(it.system());

  System refSys(it.system());

  // set atoms for which we want to compute the solvation free energy
    utils::AtomSpecifier atoms(sys);
    utils::AtomSpecifier atoms_to_charge(sys);
     utils::AtomSpecifier sphereatoms(sys);
  
    Arguments::const_iterator iter1=args.lower_bound("atoms");
    Arguments::const_iterator to=args.upper_bound("atoms");
    for(;iter1!=to;iter1++){
      atoms.addSpecifier(iter1->second.c_str());
    } 
     
    if (atoms.size()==0)  throw gromos::Exception("dGslv_pbsolv","No atoms specified. Exiting ...");

     iter1=args.lower_bound("atoms_to_charge");
     to=args.upper_bound("atoms_to_charge");
    for(;iter1!=to;iter1++){
      atoms_to_charge.addSpecifier(iter1->second.c_str());
    }
  
    // iter1=args.lower_bound("sphereatoms");
    // to=args.upper_bound("sphereatoms");
   // for(;iter1!=to;iter1++){
   //   sphereatoms.addSpecifier(iter1->second.c_str());
   // }



    if (atoms_to_charge.size()==0)  throw gromos::Exception("dGslv_pbsolv","No atoms to charge specified. Exiting ...");

   std::cout << "# READ: atoms_to_charge " << endl;
   for (int i=0;i< atoms_to_charge.size();i++){
       std::cout << "# READ: mol " <<  atoms_to_charge.mol(i) << " " << atoms_to_charge.name(i) << endl;
   }


  // if (sphereatoms.size()==0)  throw gromos::Exception("dGslv_pbsolv","No sphereatoms specified. Exiting ...");

  // std::cout << "# READ: sphereatoms " << endl;
  // for (int i=0;i< sphereatoms.size();i++){
  //     std::cout << "# READ: mol " <<  sphereatoms.mol(i) << " " << sphereatoms.name(i) << endl;
  // }


  // read scheme
    string scheme = "";
    
    if(args.count("scheme")>0){
      scheme = args["scheme"];
      transform(scheme.begin(), scheme.end(), scheme.begin(), static_cast<int (*)(int)>(std::toupper));
    }
      if(scheme!="LS" && scheme !="RF")
	throw gromos::Exception("dGslv_pbsolv","scheme format "+scheme+" unknown. Exiting ...");
      
  std::cout << "# READ: scheme " << scheme << endl;

  // read cutoff distance
  double rcut=1.4;
  if(args.count("rcut")>0) rcut=atof(args["rcut"].c_str());
  if (rcut<0)  throw gromos::Exception("dGslv_pbsolv","The cutoff must not be negative. Exiting ...");

  std::cout << "# READ: rcut " << rcut << endl;

   // read probe IAC
  int probe_iac=0;
  if(args.count("probeIAC")>0) probe_iac=atoi(args["probeIAC"].c_str());
  if (probe_iac<=0)  throw gromos::Exception("dGslv_pbsolv","The probe IAC must not be negative or 0. Exiting ...");

  std::cout << "# READ: probe_iac " << probe_iac << endl;

   // read probe radius
  double probe_rad=0;
  if(args.count("probeRAD")>0) probe_rad=atof(args["probeRAD"].c_str());
  if (probe_rad<0)  throw gromos::Exception("dGslv_pbsolv","The probe radius must not be negative. Exiting ...");

  std::cout << "# READ: probe_rad " << probe_rad << endl;


   // read hydrogen radius
  double hydrogen_rad=0.05;
  if(args.count("HRAD")>0) hydrogen_rad=atof(args["HRAD"].c_str());
  if (hydrogen_rad<0)  throw gromos::Exception("dGslv_pbsolv","The hydrogen radius must not be negative. Exiting ...");

  std::cout << "# READ: hydrogen_rad " << hydrogen_rad << endl;

  // read NWAT
 // int nwat=0;
 // if(args.count("NWAT")>0) nwat=atoi(args["NWAT"].c_str());
 // if (nwat<=0)  throw gromos::Exception("dGslv_pbsolv","The number of water molecules must be positive. Exiting ...");

//std::cout << "# READ: NWAT (only used for C1-correction estimate)" << nwat << endl;




   // read chargefac
  //double chargefac=0;
  //if(args.count("chargefac")>0) chargefac=atof(args["chargefac"].c_str());


  //std::cout << "# READ: chargefac " << chargefac << endl;


     // read radscal
  double radscal=1.0;
  if(args.count("radscal")>0) radscal=atof(args["radscal"].c_str());
  std::cout << "# READ: radscal " << radscal << endl;


       // read boxc1
 //double boxc1=0.0;
 // if(args.count("boxc1")>0) boxc1=atof(args["boxc1"].c_str());
 // std::cout << "# READ: boxc1 " << boxc1 << endl;

     // read sphere
//  double sphererad=0;
//  if(args.count("sphereRAD")>0) sphererad=atof(args["sphereRAD"].c_str());
//  if (sphererad<0)  throw gromos::Exception("dGslv_pbsolv","The sphere radius must not be negative. Exiting ...");

 // std::cout << "# READ: sphererad " << sphererad << endl;



    // read epsilon
  double epssolvent=0.0;
  if(args.count("eps")>0) epssolvent=atof(args["eps"].c_str());
  if (epssolvent<1)  throw gromos::Exception("dGslv_pbsolv","The solvent permittivity must not be smaller than 1. Exiting ...");


  std::cout << "# READ: epssolvent " << epssolvent << endl;


    // read RF epsilon
  double epsRF=0.0;
  if(args.count("epsRF")>0) epsRF=atof(args["epsRF"].c_str());
    if (epsRF<1 && epsRF != 0.0)  throw gromos::Exception("dGslv_pbsolv","The reaction field permittivity must not be smaller than 1. Exiting ...");

  std::cout << "# READ: epsRF " << epsRF << endl;


      // read NPBC epsilon
  double epsNPBC=0.0;
  if(args.count("epsNPBC")>0) epsNPBC=atof(args["epsNPBC"].c_str());
    if (epsNPBC<1)  throw gromos::Exception("dGslv_pbsolv","The NPBC permittivity must not be smaller than 1. Exiting ...");

  std::cout << "# READ: epsNPBC " << epsNPBC << endl;

 // read gridspacing
  double gridspacing=0.0;
  if(args.count("gridspacing")>0) gridspacing=atof(args["gridspacing"].c_str());
    if (gridspacing<0)  throw gromos::Exception("dGslv_pbsolv","The grid spacing must not be negative. Exiting ...");

  std::cout << "# READ: gridspacing " << gridspacing << endl;


  // read number of gridpoints
  int ngrid_x=0;
  int ngrid_y=0;
  int ngrid_z=0;
  Arguments::const_iterator iter2=args.lower_bound("gridpointsXYZ");

  if(iter2!=args.upper_bound("gridpointsXYZ")){
     // std::sstream s(iter->second);
     // s >> ngrid_x;
	ngrid_x=atoi(iter2->second.c_str());
	++iter2;
      }

      if(iter2!=args.upper_bound("gridpointsXYZ")){
        ngrid_y=atoi(iter2->second.c_str());
        ++iter2;
      }
       if(iter2!=args.upper_bound("gridpointsXYZ")){
        ngrid_z=atoi(iter2->second.c_str());
        }
        if (ngrid_x<=0 || ngrid_y <=0 || ngrid_z<=0)  throw gromos::Exception("dGslv_pbsolv","The number of grid points must be positive. Exiting ...");

      std::cout << "# READ: gridX,Y,Z " << ngrid_x << " " << ngrid_y << " " << ngrid_z << endl;

   // read maxiter
  int maxiter=0;
  if(args.count("maxiter")>0) maxiter=atoi(args["maxiter"].c_str());
  if (maxiter<=0)  throw gromos::Exception("dGslv_pbsolv","The maximum number of iterations must be positive. Exiting ...");


  std::cout << "# READ: maxiter " << maxiter << endl;


  // read radius definition
  int rminorsigma=0;
  if(args.count("rminORsigma")>0) rminorsigma=atoi(args["rminORsigma"].c_str());
   if (rminorsigma != 0 && rminorsigma != 1)  throw gromos::Exception("dGslv_pbsolv","The rminorsigma flag should be 0 (rmin) or 1 (sigma). Exiting ...");

  std::cout << "# READ: rminORsigma " << rminorsigma << endl;


     // read convergence
 // double convergence=0.0;
//  if(args.count("convergence")>0) convergence=atof(args["convergence"].c_str());
//  if (convergence <0)  throw gromos::Exception("dGslv_pbsolv","The convergence criterion must not be negative. Exiting ...");

 // std::cout << "# READ: convergence " << convergence << endl;


   // read FFT_Cubes
  int fftcub=0;
  if(args.count("FFT_Cubes")>0) fftcub=atoi(args["FFT_Cubes"].c_str());
   if (fftcub<=0)  throw gromos::Exception("dGslv_pbsolv","The FFT cubelet number must be positive. Exiting ...");

  std::cout << "# READ: fftcub " << fftcub << endl;


//Boundary *pbc = BoundaryParser::boundary(sys, args);
//Boundary::MemPtr gathmethod = args::GatherParser::parse(args);


    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);






  //cout << "# read in arguments " << endl;










      // read  coordinates
      
      InG96 ic(args["frame_g96"]);
        if(args.count("frame_g96")>0){
      ic.select("ALL");
      ic >> sys;
      (*pbc.*gathmethod)();
      ic.close();
    }

     //   cout << "# read in coordinates " << endl;





      if (!sys.hasBox)throw gromos::Exception("dGslv_pbsolv","No box block in coordinate file. Exiting ...");
      if (!sys.hasPos)throw gromos::Exception("dGslv_pbsolv","No position block in coordinate file. Exiting ...");

 





      PB_Parameters ppp(epssolvent);

// get convergence
      double convergence_fd=ppp.get_convergence_fd();
      double convergence_fft=ppp.get_convergence_fft();


      // get atomic radii (these are rmin)
      std::cout << "# RADII: based on probe_iac " << probe_iac << endl;
      if (rminorsigma == 1 ){
          std::cout << "# RADII: based on sigma" << endl;
      }
      else
      {
          std::cout << "# RADII: based on rmin" << endl;
      }
      utils::compute_atomic_radii_vdw(probe_iac-1,probe_rad,sys, it.forceField());



      // loop over atoms and set H radii to hradius_param
      // also, if sigma is used, adapt appropriately
      for (int i=0;i<atoms.size();i++){
  
          if (fabs(atoms.radius(i))< ppp.tiny_real){
     //     sys.mol(atoms.mol(i)).topology().atom(atoms.atom(i)).setradius( ppp.get_hradius() );
            sys.mol(atoms.mol(i)).topology().atom(atoms.atom(i)).setradius( hydrogen_rad );
          std::cout << "# !!!! H atom number (0 radius) " << i << " : assign rad(i) = " << atoms.radius(i) << endl;
      }
          else{
// scale radius
              if ( rminorsigma == 1){
              // use sigma instead (i.e. add probe_rad again, then divide by 2^(1/6) then subtract probe_rad)
              double tmprad=(atoms.radius(i) + probe_rad)/( exp((1.0/6.0)  * log(2.0) ) )-probe_rad;
              // now scale
              tmprad=tmprad*radscal;
              sys.mol(atoms.mol(i)).topology().atom(atoms.atom(i)).setradius( tmprad);
              }
              else{
                  // rmin ok, just scale
              double tmprad=atoms.radius(i)*radscal;
              sys.mol(atoms.mol(i)).topology().atom(atoms.atom(i)).setradius( tmprad);
              }
          }
          std::cout << "# radius of atom " << i << " : rad(i) = " << atoms.radius(i)  << endl;
      }
       for (int i=0;i<atoms_to_charge.size();i++){
          std::cout << "# radius of atom_to_charge " << i << " : rad(i) = " << atoms_to_charge.radius(i)  << endl;
      }


      // modify the charges of the atoms to be charged by multiplication with chargefac

      // for (int i=0;i<atoms_to_charge.size();i++){
      //     sys.mol(atoms.mol(i)).topology().atom(atoms.atom(i)).setCharge( atoms_to_charge.charge(i)*chargefac ) ;
//
 //          std::cout << "# modified with chargefac " << chargefac << " : charge of atom " << i << " is " << atoms_to_charge.charge(i) << endl;
  // }



      // compute a suggestion for the LS and RF C1 correction for water
     // double volbox = boxc1*boxc1*boxc1;
     // double volcutoffsphere=4.0/3.0*ppp.getPI()*rcut*rcut*rcut;
     // double volion=0.0;
     // double chargeion=0.0;
     //  for (int i=0;i<atoms_to_charge.size();i++){
     //       volion+=atoms_to_charge.radius(i) * atoms_to_charge.radius(i) * atoms_to_charge.radius(i) * 4.0/3.0*ppp.getPI();
     //      chargeion+=atoms_to_charge.charge(i);
     // }
     // double volwater=volbox-volion;
      //double numden=nwat/volwater;
      //double exclpot=ppp.get_quadr()*ppp.getFPEPSI()*4*ppp.getPI()/6.0*numden;
      //double C1_ls=exclpot*(-1.0)*chargeion*(1.0-volion/volbox);
      //double C1_rf=exclpot*(-1.0)*chargeion*(1.0-volion/volcutoffsphere)*2.0*(epsRF-1)/(2*epsRF+1);


      //std::cout << "# volbox " << volbox << " volion " << volion << " exclpot " << exclpot << endl;
      //std::cout << "# Estimate for C1_LS (using spc water as solvent): " << C1_ls << endl;
      //std::cout << "# Estimate for C1_RF (using spc water as solvent): " << C1_rf << endl;








     // sys.mol(atoms.mol(atoms.size()-1)).topology().atom(atoms.atom(atoms.size()-1)).setradius(sphererad);
     // atoms.pos(atoms.size()-1)[0]=fit::PositionUtils::cog(*atoms.sys(), sphereatoms)[0];
    //  atoms.pos(atoms.size()-1)[1]=fit::PositionUtils::cog(*atoms.sys(), sphereatoms)[1];
    //  atoms.pos(atoms.size()-1)[2]=fit::PositionUtils::cog(*atoms.sys(), sphereatoms)[2];
     

      // std::cout << "# sphereatoms cog " << atoms.pos(atoms.size()-1)[0] << " " << atoms.pos(atoms.size()-1)[1] << " " << atoms.pos(atoms.size()-1)[2]  << endl;
     //  sys.mol(atoms.mol(atoms.size()-1)).topology().atom(atoms.atom(atoms.size()-1)).setCharge(0.0);






      // FOR DEBUGGING
    //  std::cout << "!!!!!!!!!!!!!! WARNING SET RADIUS FOR DEBUGGING: 0.309640222378209" << endl;
    //        sys.mol(atoms.mol(0)).topology().atom(atoms.atom(0)).setradius( 0.309640222378209 );
      
 //cout << "# got radii" << endl;

     //   double   tmp1=ppp.get_hradius();
     //   std::cout << "# tmp1 " <<  tmp1 << endl;
     //   double   tmp2=ppp.get_alpha1();
     //   std::cout << "# tmp2 " <<  tmp2 << endl;
        
    //    std::cout << "# alpha1 " << ppp.get_alpha1() << endl;
    //    std::cout << "# alpha2 " << ppp.get_alpha2() << endl;
   //     std::cout << "# nalias1 " << ppp.get_nalias1() << endl;
   //     std::cout << "# nalias2 " << ppp.get_nalias2() << endl;

  //    std::cout << "#   LLL " << ppp.get_FFTlambda() << endl;



     
      // do two FD calculations and LS FFT in any case

          //*****************************************************************
          //DO THE FD ONES NOW


          // first, npbc
    /*
     OK */

      



    
/* MARIA COMMENT OUT FOR RF FFT DEBUG  ONLY */
      {
        FDPoissonBoltzmann_ICCG_NPBC iccg_npbc(ngrid_x,ngrid_y,ngrid_z);
        
        cout << "# ************************************************** " <<   endl;  
        cout << "# *** FD LS NPBC EPSSOLV *** " <<   endl;
        cout << "# ************************************************** " <<   endl;
        FDPoissonBoltzmann pbsolv_NPBC_epssolvent(atoms,atoms_to_charge,ngrid_x,ngrid_y,ngrid_z,gridspacing, false,epsNPBC);
        pbsolv_NPBC_epssolvent.setupGrid(true);
        pbsolv_NPBC_epssolvent.solveforpotential_npbc(maxiter, convergence_fd,iccg_npbc);
        double dgresult_LS_FD_NPBC_epssolvent = pbsolv_NPBC_epssolvent.dGelec();
       
        cout << "# ************************************************** " <<   endl;  
        cout << "# *** FD LS NPBC VAC *** " <<   endl;
        cout << "# ************************************************** " <<   endl;
        FDPoissonBoltzmann pbsolv_NPBC_vac(atoms,atoms_to_charge,ngrid_x,ngrid_y,ngrid_z,gridspacing, false,1.0); // vacuum calc
        pbsolv_NPBC_vac.setupGrid(true);
        pbsolv_NPBC_vac.solveforpotential_npbc(maxiter, convergence_fd,iccg_npbc);
        double dgresult_LS_FD_NPBC_vac = pbsolv_NPBC_vac.dGelec();
        double result1=dgresult_LS_FD_NPBC_epssolvent-dgresult_LS_FD_NPBC_vac;
        cout << "# DGRESULT NPBC " <<  result1   << endl;


         // kill objects to free memory
         cout << "# done with FD-NPBC ... now deconstructing iccg_npbc, pbsolv_NPBC_epssolvent, pbsolv_NPBC_vac"  << endl;
         //iccg_npbc.~FDPoissonBoltzmann_ICCG_NPBC();
         //pbsolv_NPBC_epssolvent.~FDPoissonBoltzmann();
       //  ~pbsolv_NPBC_epssolvent();
       //  pbsolv_NPBC_vac.~FDPoissonBoltzmann();

      }
      {
          // then pbc
         FDPoissonBoltzmann_ICCG_PBC iccg_pbc(ngrid_x,ngrid_y,ngrid_z);
     
         cout << "# ************************************************** " <<   endl;
         cout << "# *** FD LS PBC EPSSOLV *** " <<   endl;
         cout << "# ************************************************** " <<   endl;
         FDPoissonBoltzmann pbsolv_PBC_epssolvent(atoms,atoms_to_charge,ngrid_x,ngrid_y,ngrid_z,gridspacing, true,epssolvent);
         pbsolv_PBC_epssolvent.setupGrid(true);
         pbsolv_PBC_epssolvent.solveforpotential_pbc(maxiter, convergence_fd,iccg_pbc);
         double dgresult_LS_FD_PBC_epssolvent = pbsolv_PBC_epssolvent.dGelec();
         
         cout << "# ************************************************** " <<   endl;
         cout << "# *** FD LS PBC VAC *** " <<   endl;
         cout << "# ************************************************** " <<   endl;
         FDPoissonBoltzmann pbsolv_PBC_vac(atoms,atoms_to_charge,ngrid_x,ngrid_y,ngrid_z,gridspacing, true,1.0); //vacuum calc
         pbsolv_PBC_vac.setupGrid(true);
         pbsolv_PBC_vac.solveforpotential_pbc(maxiter, convergence_fd,iccg_pbc);
         double dgresult_LS_FD_PBC_vac = pbsolv_PBC_vac.dGelec();
         double result2=dgresult_LS_FD_PBC_epssolvent-dgresult_LS_FD_PBC_vac;
         cout << "# DGRESULT PBC " << result2   << endl;


         // kill objects to free memory
         cout << "# done with FD-PBC ... now deconstructing iccg_pbc, pbsolv_PBC_epssolvent, pbsolv_PBC_vac"  << endl;
      //   iccg_pbc.~FDPoissonBoltzmann_ICCG_PBC();
      //   pbsolv_PBC_epssolvent.~FDPoissonBoltzmann();
      //   pbsolv_PBC_vac.~FDPoissonBoltzmann();
      }

 /* END MARIA COMMENT OUT */
         //*****************************************************************
         // DO THE FFT ONE NOW

      FFTGridType gt(ngrid_x, ngrid_y, ngrid_z, ngrid_x*gridspacing, ngrid_y*gridspacing, ngrid_z*gridspacing, fftcub);

      {

        if (scheme == "RF"){



        // DO AN ADDITIONAL LS FFT


        cout << "# ************************************************** " <<   endl;
        cout << "# *** FFT LS (PBC) *** " <<   endl;
        cout << "# ************************************************** " <<   endl;
       // make the grid
        


        // make the boundary conditions for LS
        FFTBoundaryCondition bc_LS(0, "LS",
        ppp.get_alpha1(), ppp.get_alpha2(), ppp.get_nalias1(), ppp.get_nalias2(), rcut, epsRF);
         // (rcut and epsRF are basically dummies)
        // and print them
        bc_LS.dumpparameters();

        // setup the main object
        FFTPoisson fftp_LS(atoms, atoms_to_charge, gt, bc_LS, maxiter, convergence_fft, ppp.get_FFTlambda(),
        epssolvent, false);

       // and now we iterate
        cout << "# call solve_poisson ..." << endl;
       fftp_LS.solve_poisson();


          // kill objects to free memory
         cout << "# done with FFT-LS ... now deconstructing bc_LS, fftp_LS"  << endl;
     //    fftp_LS.~FFTPoisson();
      //   bc_LS.~FFTBoundaryCondition();

         
                      }// end of RF
      }
         //*****************************************************************

      {
       if (scheme == "RF"){



        // DO AN ADDITIONAL RF FFT



        cout << "# ************************************************** " <<   endl;
        cout << "# *** FFT RF (PBC) *** " <<   endl;
        cout << "# ************************************************** " <<   endl;

        

        if (fabs(epsRF- 1.0)> ppp.tiny_real ){
        cout << "# WITH REACTION FIELD CORRECTION : RF ... "  << endl;
        // make the boundary conditions for RF
        FFTBoundaryCondition bc_RF(1, "RF",
        ppp.get_alpha1(), ppp.get_alpha2(), ppp.get_nalias1(), ppp.get_nalias2(), rcut, epsRF);
       // setup the main objects for RF
        FFTPoisson fftp_RF(atoms, atoms_to_charge,  gt, bc_RF, maxiter, convergence_fft, ppp.get_FFTlambda(),
        epssolvent, false);
          // print params and iterate : RF
        bc_RF.dumpparameters();
        cout << "# call solve_poisson ..." << endl;
        fftp_RF.solve_poisson();
        // kill objects to free memory
         cout << "# done with FFT-RF ... now deconstructing bc_RF, fftp_RF"  << endl;
      //   fftp_RF.~FFTPoisson();
      //   bc_RF.~FFTBoundaryCondition();
        }
        else{
        cout << "# WITHOUT REACTION FIELD CORRECTION : SC ... "  << endl;
         // make the boundary conditions for SC
        FFTBoundaryCondition bc_SC(2, "SC",
        ppp.get_alpha1(), ppp.get_alpha2(), ppp.get_nalias1(), ppp.get_nalias2(), rcut, epsRF);
        // setup the main objects for RF
        FFTPoisson fftp_SC(atoms,  atoms_to_charge,  gt, bc_SC, maxiter, convergence_fft, ppp.get_FFTlambda(),
        epssolvent, false);
        // print params and iterate : SC
        bc_SC.dumpparameters();
        cout << "# call solve_poisson ..." << endl;
        fftp_SC.solve_poisson();
        // kill objects to free memory
         cout << "# done with FFT-RF ... now deconstructing bc_SC, fftp_SC"  << endl;
      //   fftp_SC.~FFTPoisson();
   //      bc_SC.~FFTBoundaryCondition();
        }


         //*****************************************************************
          } // end of RF
      }



 } // end of try for argument reading

    catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }






}
