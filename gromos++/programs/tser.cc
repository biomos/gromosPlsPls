/**
 * @file tser.cc
 * time series
 */

/**
 * @page programs Program Documentation
 *
 * @anchor tser
 * @section tser time series of properties
 * @author @ref mc
 * @date 22. 11. 2004
 *
 * calculate the time series of properties.
 * 
 * arguments:
 * - topo topology
 * - pbc [v,r,t,c] [gathermethod]
 * - time t0 dt
 * - prop [@ref PropertySpecifier "property specifier"]
 * - traj trajectory
 * 
 * <b>See also</b> @ref PropertySpecifier "property specifier"
 *
 * Example:
 * @verbatim
  tser
    @topo ex.top
    @pbc  r
    @time 0 0.1
    @prop t%1:1,3,5,6
    @traj ex.tr

    @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <sstream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/fit/Reference.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/PropertyContainer.h"
#include <vector>
#include <iomanip>
#include <math.h>
#include <iostream>

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "time", "prop", "traj", "skip", "stride"};
  int nknowns = 7;

  string usage = argv[0];
  usage += "\n\t@topo      <topology>\n";
  usage += "\t@pbc       <boundary type>\n";
  usage += "\t@time      <time and dt>\n";  
  usage += "\t@prop      <property specifier>\n";
  usage += "\t@traj      <trajectory files>\n";
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
  
    // read in a property
    
    // it's nice to actually store the properties read in somewhere
    // let's use a PropertyContainer, as it can handle some of the
    // work later on
    // the PropertyContainer should know about the system, so it can
    // check the input whether ie the atoms exist in the topology
    PropertyContainer props(sys);
    {
      Arguments::const_iterator iter=args.lower_bound("prop");
      Arguments::const_iterator to=args.upper_bound("prop");
      // we read in all properties specified by the user
      if(iter==to)
	throw Arguments::Exception("no property given");
      for(; iter!=to; iter++)
	{
	  string spec=iter->second.c_str();
	  // and that's how easy it is to add a standard property
	  // like distance, angle, torsional angle
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

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

    // define input coordinate
    InG96 ic(skip, stride);

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
	ic.select("ALL");
      
    // loop over single trajectory
    while(!ic.eof()){
      ic >> sys;
      if (ic.stride_eof()) break;
      
      (*pbc.*gathmethod)();
      
      // calculate the props
      // this is now the place, where a property-container is very handy
      // it knows that we want to loop over all properties and calculate
      // their 'value'. This works best if the property can be represented
      // by one number (which is often the case)
      props.calc();

      // print the properties
      // this is a time series, so let's just print out the properties
      // the << operator is overloaded for the property container as well
      // as for the single properties
      cout << setw(10) << time << "\t\t";
      cout << props;

      // check the boundaries
      // for every property, there is the (standard) possibility, to add
      // a range within that the property is allowed to be.
      // if there is a violation, a messge will be printed out
      cout << props.checkBounds();
      
      time += dt;
    }

  }

    // and that's all! because the properties have to worry themselves,
    // how they are calculated, it does not matter, what properties you
    // add to the container, they just have to be known (= be implemented)

    // if you need to add non-standard properties, have a look at the dist
    // program!

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



