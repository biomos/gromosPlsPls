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
 * - skip <nr> (skip nr initial structures)
 * - stride <nr> (stride through structures)
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
 @skip 0
 @stride 1

 @endverbatim
 *
 * @bug Mar 22 2005: nearestImage calls in properties were missing
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

  char *knowns[] = {"topo", "pbc", "time", "prop",
		    "traj", "skip", "stride", "nots",
		    "dist", "norm", "solv"};
  int nknowns = 11;

  string usage = argv[0];
  usage += "\n\t@topo      <topology>\n";
  usage += "\t@pbc       <boundary type>\n";
  usage += "\t@time      <time and dt>\n";  
  usage += "\t@prop      <property specifier>\n";
  usage += "\t[@nots     do not write time series]\n";
  usage += "\t[@dist     <steps [min max]>]\n";
  usage += "\t[@norm     normalise distribution\n";
  usage += "\t[@solv     read in solvent\n";
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

    bool do_tser = true;
    if (args.count("nots") >= 0)
      do_tser = false;

    // if dist_min and dist_max are
    // not supplied, the maximum value
    // (of Stat) lies outside (!) the
    // distribution and is missed.
    // Still, the behaviour seems to be
    // correct.
    // probably it would be best to add
    // a tiny number to max in this case,
    // but then: what is a tiny number?
    bool do_dist = false;
    bool dist_boundaries = false;
    double dist_min=0, dist_max=0;
    int dist_steps=0; 
    if (args.count("dist") > 0)
    {
      do_dist = true;
      Arguments::const_iterator iter=args.lower_bound("dist");
      if(iter!=args.upper_bound("dist")){
	std::istringstream is(iter->second);
	is >> dist_steps;
	++iter;
      }
      if(iter!=args.upper_bound("dist")){
	dist_boundaries = true;
	std::istringstream is(iter->second);
	is >> dist_min;
	++iter;
      }
      if(iter!=args.upper_bound("dist")){
	std::istringstream is(iter->second);
	is >> dist_max;
      }     
    }

    bool normalize = false;
    if (args.count("norm") != -1)
      normalize = true;

    //  read topology
    args.check("topo",1);
    InTopology it(args["topo"]);
    
    System sys(it.system());

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

    // read in a property
    
    // it's nice to actually store the properties read in somewhere
    // let's use a PropertyContainer, as it can handle some of the
    // work later on
    // the PropertyContainer should know about the system, so it can
    // check the input whether ie the atoms exist in the topology
    PropertyContainer props(sys, pbc);
    {
      std::string prop;
      Arguments::const_iterator iter=args.lower_bound("prop");
      Arguments::const_iterator to=args.upper_bound("prop");
      // we read in all properties specified by the user
      if(iter==to)
	throw Arguments::Exception("no property given");
      for(; iter!=to; iter++){
	string spec=iter->second.c_str();
	prop += " " + spec;
	// props.addSpecifier(spec);
      }    
      // and that's how easy it is to add a standard property
      // like distance, angle, torsional angle
      props.addSpecifier(prop);
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

    bool solvent = false;
    if (args.count("solv") != -1)
      solvent = true;

    // define input coordinate
    InG96 ic(skip, stride);

    // title
    if (do_tser){
      cout << "#" << endl;
      cout << "#" << setw(9) << "time";
      cout << "\t\t" << props.toTitle();
      cout << endl;
      std::cout.precision(6);
      std::cout.setf(std::ios::fixed, std::ios::floatfield);
    }
    
    // loop over all trajectories
    for(Arguments::const_iterator 
	  iter=args.lower_bound("traj"),
	  to=args.upper_bound("traj");
	iter!=to; ++iter){

      // open file
      ic.open((iter->second).c_str());
      if (solvent) ic.select("ALL");
      
      // loop over single trajectory
      while(!ic.eof()){
	ic >> sys;
	if (ic.stride_eof()) break;
      
	(*pbc.*gathmethod)();
      
	// calculate the props
	// this is now the place, where a property-container is very handy
	// it knows that we want to loop over all properties and calculate
	// their 'value'.
	props.calc();

	// print the properties
	// this is a time series, so let's just print out the properties
	// the << operator is overloaded for the property container as well
	// as for the single properties
	if (do_tser){
	  cout << setw(10) << time << "\t\t";
	  cout << props;
	}
	
	time += dt;
      }
    }

    ic.close();
    if (do_tser){
      cout << "# Averages over run: (<average> <rmsd> <error estimate>)\n" ;
      for(unsigned int i=0; i<props.size(); ++i){
	std::cout << "# " << props[i]->toTitle() << "\t" << props[i]->average() << "\n";
      }
    }
    
    // do we write distributions? (only scalar ones!)
    if (do_dist){
      std::cout.precision(6);
      std::cout.setf(std::ios::fixed, std::ios::floatfield);
      
      for(unsigned int i=0; i<props.size(); ++i){
	
	gmath::Stat<double> & stat = props[i]->getScalarStat();
	if (dist_boundaries)
	  stat.dist_init(dist_min, dist_max, dist_steps);
	else
	  stat.dist_init(dist_steps);

	cout << "#" << endl;  
	cout << "# Distribution of     " << props[i]->toTitle() << endl;
	cout << "# values:             " << stat.n() << endl;
	cout << "# average value:      " << stat.ave() << endl;
	cout << "# rmsd:               " << stat.rmsd() << endl;
	cout << "# error estimate:     " << stat.ee() << endl;

	if (normalize)
	  stat.distribution().write_normalized(cout);
	else
	  stat.distribution().write(cout);
      }
    }
  }
  catch (const gromos::Exception &e){
    cerr << "EXCEPTION:\t";
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}



