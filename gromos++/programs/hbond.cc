//Hydrogen bond analysis Hbond
//dedicated to wilfred: "you know, mika: proahb IS my favorite program"
//--mika

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/Neighbours.h"
#include "../src/utils/Hbond.h"
#include "../src/utils/Hbondcalc.h"


#include <iostream>
#include <iomanip>



using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace gmath;
using namespace utils;
using namespace std;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "type",  "ref", "SoluteDonorAtoms", "SoluteAcceptorAtoms", "SolventDonorAtoms", "SolventAcceptorAtoms","Hbparas", "time", "massfile", "traj", "molrange"};
  int nknowns = 13;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t@type   <1 2 3 4 5 6>\n"
    "\t\tChoices are:\n"
    "\t\t1: Intramolecular\n"
    "\t\t2: Intermolecular\n"
    "\t\t3: native Intramolecular\n"
    "\t\t4: native Intermolecular\n"
    "\t\t5: Solute <-> Solvent\n"
    "\t\t6: Solvent <-> Solvent (SLOW...)\n";
  usage += "\t@SoluteDonorAtoms <atomspecifier>\n";
  usage += "\t@SoluteAcceptorAtoms <atomspecifier>\n";
  usage += "\t@SolventDonorAtoms <atomspecifier>\n";
  usage += "\t@SolventAcceptorAtoms <atomspecifier>\n";
  usage += "\t@ref <reference coordinates for native H-bonds>\n";
  usage += "\t@Hbparas <distance [nm] and angle; default: 0.25, 135>\n";
  usage += "\t@time   <time and dt>\n";
  usage += "\t[@massfile   <massfile>]\n";
  usage += "\t[@molrange   <number of molecules and atom range, e.g. 2 1-10 11-20. >]\n";
  usage += "\t@traj   <trajectory files>\n";
  
 
    try{
    Arguments args(argc, argv, nknowns, knowns, usage);
   

  args.check("topo",1);
  InTopology it(args["topo"]);
  System sys(it.system());


  Hbondcalc HB(sys,args);

  //get time

  double time=0, dt=1; 
  {
    Arguments::const_iterator iter=args.lower_bound("time");
    if(iter!=args.upper_bound("time")){
      time=atof(iter->second.c_str());
      ++iter;
    }
    if(iter!=args.upper_bound("time"))
        dt=atof(iter->second.c_str());
  }

  HB.settime(time, dt);

  // get the paras
  double maxdist=0.25, minangle = 135;
  {
    Arguments::const_iterator iter=args.lower_bound("Hbparas");
    if(iter!=args.upper_bound("Hbparas")){
      maxdist=atof(iter->second.c_str());
      ++iter;
    }
    if(iter!=args.upper_bound("Hbparas"))
        minangle=atof(iter->second.c_str());
  }

  HB.setmaxdist(maxdist);
  HB.setminangle(minangle);

  HB.determineAtoms();

 
  //check for massfile
  string mfile;
  try{
   args.check("massfile");
   {
   Arguments::const_iterator iter=args.lower_bound("massfile");
   if(iter!=args.upper_bound("massfile")){
     mfile=(iter->second.c_str());
   }
   }
   HB.readinmasses(mfile);
   HB.determineAtomsbymass();
  }
   catch(Arguments::Exception e){
   }


  //get type
   int type=0;
     {
   Arguments::const_iterator iter=args.lower_bound("type");
   if(iter!=args.upper_bound("type")){
     type=atoi(iter->second.c_str());
   }
    }

   
  Hbondcalc::MemPtr calcmethod;
  

  //switch statement

  switch(type){
  case 1: { HB.calcHintra_init(); calcmethod = &Hbondcalc::calcHintra; }
  break;
  case 2: { HB.calcHinter_init(); calcmethod = &Hbondcalc::calcHinter; }
  break;
  case 3: { HB.calcHintra_native_init(); calcmethod = &Hbondcalc::calcHintra; }
  break;
  case 4: { HB.calcHinter_native_init(); calcmethod = &Hbondcalc::calcHinter; }
  break;
  case 5: { HB.calcHsolusolv_init(); calcmethod = &Hbondcalc::calcHsolusolv; }
  break;
  case 6: { HB.calcHsolvsolv_init(); calcmethod = &Hbondcalc::calcHsolvsolv; }
  break;
  default:
  throw gromos::Exception("Hond", args["type"] + 
                                " unknown. Known types are 1, 2, 3, 4, 5 and 6");
  }


  InG96 ic;

  // loop over all trajectories
      for(Arguments::const_iterator 
      iter=args.lower_bound("traj"),
      to=args.upper_bound("traj");
      iter!=to; ++iter){

     // open file
       ic.open((iter->second).c_str());

     // loop over single trajectory
   
       while(!ic.eof()){

       ic >> sys;
 
       (HB.*calcmethod)();

       }
       ic.close();

      }

      std::cout << "Statistics of the run:" << endl;
      std::cout << "#"  
             << setw(8)  <<  "#HB"   
             << setw(8)  <<   "Donor" << setw(8) << "Acceptor" 
             << setw(10) <<   "D-"    << setw(5) << "H..."    <<  setw(5) << "A" 
             << setw(8)  << "DIST"    << setw(8) << "ANGLE"   <<  setw(8) << "OCCUR"
             << setw(8)  << "%"       << endl;
 

     if (type == 5) HB.printstatistics_solusolv();
     else if (type == 6) {}
     else HB.printstatistics();
  
    }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
