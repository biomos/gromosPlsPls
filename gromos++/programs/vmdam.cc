// convert.cc

#include "../src/args/Arguments.h"
#include "../src/utils/Rmsd.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gio/OutPdb.h"
#include "../src/gio/Outvmdam.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/TruncOct.h"
#include "../src/bound/Vacuum.h"
#include "../src/bound/RectBox.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/args/ReferenceParser.h"
#include <vector>
#include <iomanip>
#include <fstream>
#include <iostream>

using namespace std;
using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "traj", "class", "atoms", "time", "nframes", 
		    "pbc", "ref", "mol", "outformat"};
  int nknowns = 10;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@mol <molecules to be considered for fit>\n";
  usage += "\t@class <classes of atoms to consider for fit>\n";
  usage += "\t@atoms <atoms to consider to consider for fit>\n";
  usage += "\t@time <time and dt>\n";
  usage += "\t@nframes <total number of frames in trajectory>\n";
  usage += "\t@ref <reference coordinates>\n";
  usage += "\t@traj <trajectory files>\n";
  usage += "\t@outformat <output format>\n";


  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // get simulation time
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

    double nf=0;
    {
      Arguments::const_iterator iter=args.lower_bound("nframes");
            if(iter!=args.upper_bound("nframes"))
        nf=atof(iter->second.c_str());
    }

    // read topology
    InTopology it(args["topo"]);
    System refSys(it.system());
    
    // read reference coordinates...
    InG96 ic;

    try{
      // Adding references
      args.check("ref",1);
      ic.open(args["ref"]);
    }
    catch(const Arguments::Exception &){
      args.check("traj",1);
      ic.open(args["traj"]);
    }
      ic >> refSys;
      ic.close();

    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(refSys, args);
    // GatherParser
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

    // gather reference system
    (*pbc.*gathmethod)();
 
    delete pbc;
  
    // parse the reference system
    Reference ref(&refSys);
    ReferenceParser refP(refSys, args, ref);
    refP.add_ref();

    // System for calculation
    System sys(refSys);
    // Parse boundary conditions for sys
    pbc = BoundaryParser::boundary(sys, args);

    OutCoordinates *oc;
    
    try{
      string format = args["outformat"];
      if(format == "pdb")
	oc = new OutPdb();
      else if(format == "g96")
	oc = new OutG96S();
      else if(format == "vmdam")
        oc = new Outvmdam();
      else 
	throw gromos::Exception("Convert","output format "+format+" unknown.\n");
    }
    catch(Arguments::Exception &){
      oc = new OutG96S();
    }

  
    RotationalFit rf(&ref);
    // loop over all trajectories
cout << "Number of config.:   " << nf << "," << "Initial time:   " << time << "," << "Time between config.:  " << dt << "\n";

    int numFrames = 0;
    Vec bla(0.0,0.0,0.0);
    for(Arguments::const_iterator iter=args.lower_bound("traj");
	iter!=args.upper_bound("traj"); ++iter){
     ic.open(iter->second);
      // loop over all frames
      while(!ic.eof()){
       numFrames++;

	oc->open(cout);

	ic >> sys;
	(*pbc.*gathmethod)();
		try{
	  args.check("ref");

	    rf.fit(&sys);
		}
		catch (Arguments::Exception &){}

    
	if (numFrames == 1){
         ofstream os("ref.pdb");
         OutPdb opdb(os);
         opdb << sys;
         os.close();
	}
	*oc << sys;
	oc->close();
     
      }
      ic.close();

    }

    
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

