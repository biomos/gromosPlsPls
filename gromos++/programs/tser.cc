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
// i will use properties
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
  
    // read in a property
    // we should get at least one property
    if (args.check("prop") < 1)
      throw Arguments::Exception("no property given");
    
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
      for(; iter!=to; iter++)
	{
	  string spec=iter->second.c_str();
	  // and that's how easy it is to add a standard property
	  // like distance, angle, torsional angle
	  props.addSpecifier(spec);
	}    
    }
  
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);

    // define input coordinate
    InG96 ic;
    //define pi
    // pi is no longer needed here, 'cause the properties know themselves
    // how to be calculated
    //const double pi=3.1415926535898;
   
    // title

    cout << "#" << endl;
  
    cout << "#" << setw(9) << "time";
    // this will write out a title line, properties are abbreviated in
    // the same style as the user input...
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



