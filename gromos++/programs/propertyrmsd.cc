//distance deviations: lower-bound upper-bound rmsd

// this is the third of a triade of programs
// tser, dist and propertyrmsd

// this program obviously gives rmsd, upper and lower bound and average of a property

// read tser.cc for a simple introduction of properties
// read dist.cc if you want to add your own properties

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/bound/TruncOct.h"
#include "../src/bound/Vacuum.h"
#include "../src/bound/RectBox.h"
#include "../src/fit/Reference.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Distribution.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/PropertyContainer.h"
#include <vector>
#include <iomanip>
#include <math.h>
#include <iostream>
#include <sstream>

using namespace std;
using namespace fit;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "time", "prop", "traj", "skip", "stride"};
  int nknowns = 7;

  string usage = argv[0];
  usage += "\n\t@topo   <topology>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t@time   <start time and timestep>\n";
  usage += "\t@prop   <property specifier>\n";
  usage += "\t@traj   <trajectory files>\n";
  usage += "\t[@skip     <skip n first frames>\n";
  usage += "\t[@stride   <take every n-th frame>\n";
  
 
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

  // Parse boundary conditions
  Boundary *pbc = BoundaryParser::boundary(sys, args);
  //parse gather method
  Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

  // get properties into PropertySpecifier
  // we should at least get one property
  if (args.count("prop") <= 0)
    throw Arguments::Exception("specify at least one property");
  
  // declare a property container and read the given properties into it
  PropertyContainer props(sys, pbc);
  {
    Arguments::const_iterator iter=args.lower_bound("prop");
    Arguments::const_iterator to=args.upper_bound("prop");
    for(; iter!=to; iter++)
      {
	string spec=iter->second.c_str();
	props.addSpecifier(spec);
      }    
  }

  int skip = 0;
  if (args.count("skip") > 0){
    std::istringstream is(args["skip"]);
    if (!(is >> skip))
      throw Arguments::Exception("could not read skip");
  }
  
  int stride = 1;
  if (args.count("stride") > 0){
    std::istringstream is(args["stride"]);
    if (!(is >> stride))
      throw Arguments::Exception("could not read stride");      
  }

  // define input coordinate
  InG96 ic(skip, stride);

  // put out some titles
  cout << endl;
  cout << "# Average over properties (per timestep)\n";
  // this will print which properties are calculated
  // note that all will be added to the same distribution!
  // cout << "# " << props.toTitle() << endl;

  cout << "# PROPERTIES\n#\n";
  for(size_t i=0; i<props.size(); ++i){
    cout << "#\tProperty " << i+1 << ": "
	 << props[i]->toTitle() << "\n";
  }
  cout << "#\n# END\n\n";
  
  // this are the properties of the distribution that are calculated
  cout << "#" << setw(9) << "time"
       << setw(19) << "average " << setw(13) << "rmsd "
       << setw(13) << "rmsd-zval " << setw(13) << "lowest "
       << setw(13) << "highest " 
       << setw(7) << "lowid " << setw(7) << "highid " << "\n";
  
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
      if (ic.stride_eof()) break;

      // take care with gathering -> ask mika,chris
      (*pbc.*gathmethod)();
      
      // calculate the properties
      // the container will loop through the properties and add the calculated value
      // to the distribution
      props.calc();
      // print out boundary violations
      // this checks whether any of the calculated values lies outside the range, which
      // the user can specify when giving the property-specifier
      cout << props.checkBounds();
      // and print out the averageOverProperties
      // this is of all properties per time step
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





