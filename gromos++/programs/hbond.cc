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
#include "../src/utils/Hbond3c.h"
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

  char *knowns[] = {"topo", "pbc", "ref", "DonorAtomsA", "AcceptorAtomsA", "DonorAtomsB", "AcceptorAtomsB","Hbparas", "threecenter", "time", "massfile", "traj"};
  int nknowns = 12;
  
  string usage = argv[0];
  usage += "\n\t@topo           <topology>\n";
  usage += "\t@pbc            <boundary type>\n";
  usage += "\t@DonorAtomsA    <atomspecifier>\n";
  usage += "\t@AcceptorAtomsA <atomspecifier>\n";
  usage += "\t@DonorAtomsB    <atomspecifier>\n";
  usage += "\t@AcceptorAtomsB <atomspecifier>\n";
  usage += "\t[@ref           <reference coordinates for native H-bonds>]\n";
  usage += "\t@Hbparas        <distance [nm] and angle; default: 0.25, 135>\n";
  usage += "\t@threecenter    <distances [nm]> <angles> <sum> <dihedral>\n";
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
    // and for three-center hydrogen bonds
    double maxdist3c=0.27, minangle3c = 90;
    double minanglesum3c = 340, maxdihedral3c = 15;
    bool hbond3c=false;
    if(args.count("threecenter")>=0){
      hbond3c=true;
      Arguments::const_iterator iter=args.lower_bound("threecenter");
      if(args.count("threecenter")>0){
	maxdist3c = atof(iter->second.c_str());
	++iter;
      }
      if(args.count("threecenter")>1){
	minangle3c = atof(iter->second.c_str());
	++iter;
      }
      if(args.count("threecenter")>2){
	minanglesum3c = atof(iter->second.c_str());
	++iter;
      }
      if(args.count("threecenter")>3){
	maxdihedral3c = atof(iter->second.c_str());
      }
    }
    
    HB.setmaxdist(maxdist);
    HB.setminangle(minangle);
    HB.setmaxdist3c(maxdist3c);
    HB.setminangle3c(minangle3c);
    HB.setminanglesum3c(minanglesum3c);
    HB.setmaxdihedral3c(maxdihedral3c);
    HB.determineAtoms();
    
    // print out the information
    if(hbond3c){
      std::cout << "#\n"
		<< "# 2-Centered hydrogen bond D-H..A counted if:\n"
		<< "#     Distance H..A is at most " << maxdist << "\n"
		<< "#     Angle    D-H..A is at least " << minangle << "\n";
      std::cout << "#\n"
		<< "# 3-Centered hydrogen bond D-H..A1\n"
		<< "#                              \\A2 counted if:\n"
		<< "#     Distance H..A1 is at most " << maxdist3c << "\n"
		<< "#     Distance H..A2 is at most " << maxdist3c << "\n"
		<< "#     Angle    D-H..A1 is at least " << minangle3c << "\n"
		<< "#     Angle    D-H..A2 is at least " << minangle3c << "\n"
		<< "#     Sum of angles D-H..A1, D-H..A2, A1..H..A2 is at least "
		<< minanglesum3c << "\n"
		<< "#     Dihedral angle D..A1..A2..H is at most " << maxdihedral3c << "\n";
      std::cout << "#\n#\n";
    }
      
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
	else if(hbond3c)
	  HB.calc3c();
	else
	  HB.calc();
	
	HB.writets();
      }
      ic.close();
      
    }
    
    std::cout << "# Statistics of the run:" << endl;
    if(hbond3c)
      std::cout << endl << "# Two-centered hydrogen bonds:\n";
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

    if(hbond3c){
      std::cout << endl << endl
		<< "# Three-centered hydrogen bonds:" << endl
		<< "#" 
		<< setw(82) << " "
		<< setw(12) << "2-CENTER"
		<< setw(30) << "3-CENTER" << endl
		<< "#"
		<< setw(3) << "HB"  
		<< setw(16)<< "Donor" 
		<< setw(13)<< "Acceptor" 
		<< setw(12)<<   "D -"    
		<< setw(12)<< "H ..."
		<< setw(10) << "A" 
		<< setw(8) << "DIST"    
		<< setw(8) << "ANGLE"   
		<< setw(8) << "OCCUR"
		<< setw(8) << "%"       
		<< setw(10) << "SUM"
		<< setw(8) << "DIHED."
		<< setw(8) << "OCCUR"
		<< setw(8) << "%" << endl;
      
      HB.printstatistics3c();
    }
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
