//time series tser

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/fit/Reference.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/PropertyContainer.h"
#include <vector>
#include <iomanip>
#include <math.h>
#include <iostream>

using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "time", "prop", "traj"};
  int nknowns = 5;

  string usage = argv[0];
  usage += "\n\t@topo   <topology>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t@time   <time and dt>\n";  
  usage += "\t@prop   <property specifier>\n";
  usage += "\t@traj   <trajectory files>\n";
  
 
    try{
    Arguments args(argc, argv, nknowns, knowns, usage);
 

    //   get simulation time
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
  
    //  read topology
    args.check("topo",1);
    InTopology it(args["topo"]);
    System sys(it.system());
    
    // get atoms into AtomSpecifier
    PropertyContainer props(sys);
    {
      Arguments::const_iterator iter=args.lower_bound("prop");
      Arguments::const_iterator to=args.upper_bound("prop");
      for(; iter!=to; iter++)
	{
	  string spec=iter->second.c_str();
	  props.addSpecifier(spec);
	}    
    }
  
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);

    // define input coordinate
    InG96 ic;
    //define pi
    //const double pi=3.1415926535898;
   
    // title

    cout << "#" << endl;
  
    cout << "#" << setw(9) << "time";
    cout << "\t\t" << props.toTitle();
    cout << endl;
    
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
      pbc->gather();
   
      // calculate the props
      props.calc();

      // print the properties      
      cout << setw(10) << time << "\t\t";
      cout << props;
      
      // check the boundaries
      cout << props.checkBounds();
      
      time += dt;
    }

  }
  ic.close();
  cout << "# Averages over run: (<average> <rmsd from z-value>)\n#\t" 
       << props.averageOverRun() << endl;
  }
  catch (const gromos::Exception &e){
    cerr << "EXCEPTION:\t";
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
