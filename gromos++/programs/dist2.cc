//distributions dist

#include "../src/args/Arguments.h"
#include "../src/bound/TruncOct.h"
#include "../src/bound/Vacuum.h"
#include "../src/bound/RectBox.h"
#include "../src/args/BoundaryParser.h"
#include "../src/fit/Reference.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Distribution.h"
#include "../src/utils/PropertyContainer.h"
#include <vector>
#include <iomanip>
#include <math.h>
#include <iostream>

using namespace fit;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "prop", "dist", "traj"};
  int nknowns = 5;

  string usage = argv[0];
  usage += "\n\t@topo   <topology>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t@dist   <lower and upper boundary and number of steps>\n";
  usage += "\t@prop  <propertyspecifier>\n";
  usage += "\t@traj   <trajectory files>\n";
  
 
try{
  Arguments args(argc, argv, nknowns, knowns, usage);

  //   get distribution parameters
  double begin=0, end=0;
  int nsteps=0; 
  {
    Arguments::const_iterator iter=args.lower_bound("dist");
    if(iter!=args.upper_bound("dist")){
      begin=atof(iter->second.c_str());
      ++iter;
    }
    if(iter!=args.upper_bound("dist")){
      end=atof(iter->second.c_str());
      ++iter;
    }
    if(iter!=args.upper_bound("dist")){
      nsteps=atoi(iter->second.c_str());
    }     
  }
  
  //  read topology
  args.check("topo",1);
  InTopology it(args["topo"]);
  System sys(it.system());

  // get properties into PropertySpecifier
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
  
  // Parse boundary conditions
  Boundary *pbc;
  try{
    char b=args["pbc"].c_str()[0];
    switch(b){
      case 't':
        pbc=new TruncOct(&sys);
        break;
      case 'v':
        pbc=new Vacuum(&sys);
        break;
      case 'r':
        pbc=new RectBox(&sys);
        break;
      default:
        throw gromos::Exception("Boundary", args["pbc"] + 
				" unknown. Known boundaries are t, r and v");
	
    }
  }
  catch(Arguments::Exception &e){
    pbc = new Vacuum(&sys);
  }

  // define input coordinate
  InG96 ic;

  // set up distribution arrays
  gmath::Distribution dist(begin, end, nsteps);
  props.addDistribution(dist);

  // set pi
  // const double pi = 3.1415926535898;
  
  // loop over all trajectories
  for(Arguments::const_iterator 
        iter=args.lower_bound("traj"),
        to=args.upper_bound("traj");
      iter!=to; ++iter){

    // open file
    ic.open((iter->second).c_str());
       // ic.select("ALL");
    
      // loop over single trajectory
    while(!ic.eof()){
      ic >> sys;
      pbc->gather();
      props.calc();
      cout << props.checkBounds();
      
    }
  }
  ic.close();
  // print out the distribution, calculate the average and rmsd
  cout << "#" << endl;  
  cout << "# number of values calculated: "   
       << props.getDistribution().nVal() << endl;
  cout << "# average value:               "   
       << props.getDistribution().ave() << endl;
  cout << "# RMSD (from distribution):    "   
       << props.getDistribution().rmsd() << endl;

  cout << "# time\t\t" <<  props.toTitle() << endl;

  props.getDistribution().write(cout);
    
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}





