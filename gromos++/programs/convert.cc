// convert.cc

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/args/ReferenceParser.h"
#include "../src/utils/Rmsd.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gio/OutPdb.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/TruncOct.h"
#include "../src/bound/RectBox.h"
#include "../src/bound/Vacuum.h"
#include <vector>
#include <iomanip>
#include <fstream>
#include <iostream>

using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;
using namespace std;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "traj", "class", "atoms", 
		    "pbc", "ref", "mol", "out", "outformat"};
  int nknowns = 9;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@mol <molecules to be considered for fit>\n";
  usage += "\t@class <classes of atoms to consider for fit>\n";
  usage += "\t@atoms <atoms to consider to consider for fit>\n";
  usage += "\t@ref <reference coordinates>\n";
  usage += "\t@traj <trajectory files>\n";
  usage += "\t@outformat <output format>\n";
  usage += "\t@out <output file names>\n";


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

    // read topology
    InTopology it(args["topo"]);
    System refSys(it.system());
    
    // read reference coordinates...
    InG96 ic;

    Reference ref(&refSys);
    ReferenceParser refP(refSys, args, ref);
    refP.add_ref();
    try{
      // Adding references
      args.check("ref",1);
      ic.open(args["ref"]);
      ic >> refSys;
      ic.close();
      
      ofstream os("ref.pdb");
      OutPdb oc(os);
      oc << refSys;

    }
    catch(Arguments::Exception &){}

    // System for calculation
    System sys(refSys);

    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);;
    // Parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

    OutCoordinates *oc;
    
    try{
      string format = args["outformat"];
      if(format == "pdb")
	oc = new OutPdb();
      else if(format == "g96")
	oc = new OutG96S();
      else 
	throw gromos::Exception("Convert","output format "+format+" unknown.\n");
    }
    catch(Arguments::Exception &){
      oc = new OutG96S();
    }

    RotationalFit rf(&ref);
    
    
    // loop over all trajectories
    Arguments::const_iterator oit=args.lower_bound("out");

    for(Arguments::const_iterator iter=args.lower_bound("traj");
	iter!=args.upper_bound("traj"); ++iter){
      ic.open(iter->second);
      
      // loop over all frames
      while(!ic.eof()){
	if(oit==args.upper_bound("out"))
	  throw gromos::Exception("convert","not all frames considered;\ntoo few output file names.");
	cerr << "writing to " << oit->second.c_str() << endl;
	
	ofstream os(oit->second.c_str());
	oc->open(os);
	oc->writeTitle("Converted from " + iter->second);

	ic >> sys;
	(*pbc.*gathmethod)();
	
	try{
	  args.check("ref");
	  //	  tf.fitToCog(&sys);
       	  rf.fit(&sys);
	}
	catch (Arguments::Exception &){}
	
	*oc << sys;
	oc->close();
	os.close();
	++oit;
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

