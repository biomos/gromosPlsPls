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

  char *knowns[] = {"topo", "pbc", "ref", "DonorAtomsA", "AcceptorAtomsA", "DonorAtomsB", "AcceptorAtomsB","Hbparas", "time", "massfile", "traj"};
  int nknowns = 11;
  
  string usage = argv[0];
  usage += "\n\t@topo           <topology>\n";
  usage += "\t@pbc            <boundary type>\n";
  usage += "\t@DonorAtomsA    <atomspecifier>\n";
  usage += "\t@AcceptorAtomsA <atomspecifier>\n";
  usage += "\t@DonorAtomsB    <atomspecifier>\n";
  usage += "\t@AcceptorAtomsB <atomspecifier>\n";
  usage += "\t[@ref           <reference coordinates for native H-bonds>]\n";
  usage += "\t@Hbparas        <distance [nm] and angle; default: 0.25, 135>\n";
  usage += "\t@time           <time and dt>\n";
  usage += "\t[@massfile      <massfile>]\n";
  usage += "\t@traj           <trajectory files>\n";
  
 
  try{
    Arguments args(argc, argv, nknowns, knowns, usage);
    
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
    if(args.count("massfile")>0){
      Arguments::const_iterator iter=args.lower_bound("massfile");
      if(iter!=args.upper_bound("massfile")){
	mfile=(iter->second.c_str());
      }
      HB.readinmasses(mfile);
      HB.determineAtomsbymass();
    }
    

    // initialize the calculation
    HB.init();
    
    // do native?
    bool do_native=false;
    if(args.count("ref")>0){
      InG96 ic(args["ref"]);
      ic.select("ALL");
      ic >> sys;
      ic.close();
      
      // calculate the hb.
      HB.calc();
      // and clear the statistics and reset the time
      HB.clear();
      HB.settime(time, dt);

      do_native=true;
    }
    
    InG96 ic;
    
    // loop over all trajectories
    for(Arguments::const_iterator 
	  iter=args.lower_bound("traj"),
	  to=args.upper_bound("traj");
	iter!=to; ++iter){
      
      // open file
      ic.open((iter->second).c_str());
      ic.select("ALL");
      
      // loop over single trajectory
      
      while(!ic.eof()){
	
	ic >> sys;
	if(do_native)
	  HB.calc_native();
	else
	  HB.calc();

	HB.writets();
      }
      ic.close();
      
    }
    
    std::cout << "Statistics of the run:" << endl;
    std::cout << "#" 
	      << setw(3) << "HB"  
	      << setw(16)<< "Donor" 
	      << setw(13)<< "Acceptor" 
	      << setw(12)<<   "D -"    
	      << setw(12)<< "H ..."
	      << setw(10) << "A" 
	      << setw(8) << "DIST"    
	      << setw(8) << "ANGLE"   
	      << setw(8) << "OCCUR"
	      << setw(8) << "%"       << endl;
    
    
    HB.printstatistics();
    
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
