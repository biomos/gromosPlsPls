// rmsd.cc

#include <cassert>
#include <vector>
#include <iostream>
#include <iomanip>

#include "../src/args/Arguments.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/args/ReferenceParser.h"
#include "../src/utils/Rmsd.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gio/OutPdb.h"
#include "../src/gmath/Vec.h"


using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;
using namespace std;
using namespace gmath;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "traj", "atomsfit", "atomsrmsd", "pbc", "ref", "time"};
  int nknowns = 7;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@time <time and dt>\n";
  usage += "\t@atomsrmsd <atomspecifier: atoms to consider for rmsd>\n";
  usage += "\t[@atomsfit  <atomspecifier: atoms to consider for fit>]\n"; 
  usage += "\t@ref <reference coordinates>\n";
  usage += "\t@traj <trajectory files>\n";


  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);

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
    try{
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
    
    Reference reffit(&refSys);
    // System for calculation
    System sys(refSys);
    Reference refrmsd(&refSys);
    AtomSpecifier fitatoms(refSys);
    AtomSpecifier rmsdatoms(sys);

    //get rmsd atoms
    {
       Arguments::const_iterator iter = args.lower_bound("atomsrmsd");
       Arguments::const_iterator to = args.upper_bound("atomsrmsd");

       for(;iter!=to;iter++){
        string spec=iter->second.c_str();
        rmsdatoms.addSpecifier(spec);
       }
    }  

    refrmsd.addAtomSpecifier(rmsdatoms);
    
    //try for fit atoms
    try{
      args.check("fitatoms",1);
     
      {
       Arguments::const_iterator iter = args.lower_bound("atomsfit");
       Arguments::const_iterator to = args.upper_bound("atomsfit");

       for(;iter!=to;iter++){
        string spec=iter->second.c_str();
        fitatoms.addSpecifier(spec);
       }
      }

    }
    catch(const Arguments::Exception &){
      fitatoms = rmsdatoms;
    }
    
    reffit.addAtomSpecifier(fitatoms);

    // Parse boundary conditions for sys
    pbc = BoundaryParser::boundary(sys, args);

    RotationalFit rf(&reffit);
   
    Rmsd rmsd(&refrmsd);


    // loop over all trajectories
    for(Arguments::const_iterator iter=args.lower_bound("traj");
	iter!=args.upper_bound("traj"); ++iter){
      ic.open(iter->second);
      
      // loop over all frames
      while(!ic.eof()){
	ic >> sys;
	
	(*pbc.*gathmethod)();

	rf.fit(&sys);

	double r = rmsd.rmsd(sys);
	cout.precision(2);
	cout << setw(10) << time;
	cout.precision(5);
	cout << setw(10) << r << endl;
	time+=dt;
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

