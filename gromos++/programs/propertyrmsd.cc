//distance deviations: lower-bound upper-bound rmsd

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

  char *knowns[] = {"topo", "pbc", "time", "prop", "traj"};
  int nknowns = 5;

  string usage = argv[0];
  usage += "\n\t@topo   <topology>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t@time   <start time and timestep>\n";
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
				" unknown. Known pbc are t, r and v");
	
    }
  }
  catch(Arguments::Exception &e){
    pbc = new Vacuum(&sys);
  }

  // define input coordinate
  InG96 ic;

  cout << endl;
  cout << "# Average over properties (per timestep)\n";
  cout << "# " << props.toTitle() << endl;
  
  cout << "#" << setw(9) << "time" << setw(10) << "average" << setw(10)
       << "rmsd" << setw(10) << "rmsd-zval" << setw(10) << "lowest" << setw(10) 
       << "highest" << endl;
  
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
      
      // calculate the properties
      props.calc();
      // print out boundary violations
      cout << props.checkBounds();
      // and print out the averageOverProperties
      cout << setw(10) << time << "\t" 
	   << props.averageOverProperties() << endl;
      
      time += dt;
    }
  }
  
  ic.close();

  cout << "# avg,rmsd,... has been calculated of " 
       << props.size() << " values (per timestep).\n";
  
  } 
  
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}





