/**
 * @file dist.cc
 * @page programs Program Documentation
 *
 * @anchor dist
 * @section dist distribution of properties
 * @author @ref mc
 * @date 22. 11. 2004
 *
 * <h3>calculate the distribution of properties</h3>
 * 
 * <h4>arguments:</h4>
 * - topo <topology>
 * - pbc <v,r,t,c> [<gathermethod>]
 * - dist <lower bound  upper bound  steps>
 * - prop <property specifier> [<property specifier>] [...]
 * - traj <trajectory>
 * - norm normalise distribution
 * - solv read in solvent as well
 * 
 * <b>See also</b> @ref PropertySpecifier "property specifier"
 *
 * <h4>Example:</h4>
 * @verbatim
  dist
    @topo ex.top
    @pbc  r
    @dist 0 360 361
    @prop t%1:1,3,5,6 t%2:1,2,6,8
    @traj ex.tr
    @norm

    @endverbatim
 *
 * <h4>test reports</h4>
 * - 8/2/2002: Uses the same scheme to calculate the properties as tser.
 * Has been tested against protcf and with specific test cases.
 * Versions of before 8/2/2002 (cvs versions 1.1 and 1.2) contain errors
 * in the calculation of dihedral angles.
 *
 * <hr>
 */

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/bound/TruncOct.h"
#include "../src/bound/Vacuum.h"
#include "../src/bound/RectBox.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
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
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "prop", "dist", "traj", "norm", "solv", "skip", "stride"};
  int nknowns = 9;

  string usage = argv[0];
  usage += "\n\t@topo   <topology>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t@dist   <lower and upper boundary and number of steps>\n";
  usage += "\t@prop   <propertyspecifier>\n";
  usage += "\t@traj   <trajectory files>\n";
  usage += "\t[@norm   normalize the distribution\n";
  usage += "\t[@solv   read in solvent as well\n";
  usage += "\t[@skip     <skip n first frames>\n";
  usage += "\t[@stride   <take every n-th frame>\n";
 
try{
  Arguments args(argc, argv, nknowns, knowns, usage);

  //   get distribution parameters
  //   maybe change to a statistics class, get distribution from there...
  //   CHRIS???
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

  // Parse boundary conditions and get gather method
  Boundary *pbc = BoundaryParser::boundary(sys, args); 
  Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

  // get properties into PropertySpecifier
  // these are the standard properties we want to calculate
  // at every timestep
  // these will be added to the standard PropertyContainer distribution
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

  // set up distribution arrays
  gmath::Distribution dist(begin, end, nsteps);
  props.addDistribution(dist);

  // define input coordinate
  InG96 ic(skip, stride);

  bool normalize = false;
  if (args.count("norm") != -1)
    normalize = true;

  bool solvent = false;
  if (args.count("solv") != -1)
		 solvent = true;

  // the "real" average (distribution takes the middle of the bins)
  // CHRIS: this would also be better with a statistics class
  double average = 0.0;
  int steps = 0;
  
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
      props.calc();
      cout << props.checkBounds();
      double av, rmsd, zrmsd, lb, ub;
      int lp, up;
      
      props.averageOverProperties(av, rmsd, zrmsd, lb, ub, lp, up);
      average += av;
      ++steps;
    }
  }
  ic.close();
  // this already has been the main program!!!

  // print out the distribution, calculate the average and rmsd
  cout << "#" << endl;  
  cout << "# values per frame: " << props.size() << endl
       << "# frames: " << props.getDistribution().nVal() / props.size() << endl;
  cout << "# total number of values calculated: "   
       << props.getDistribution().nVal() << endl;
  cout << "# average value:               "   
       << props.getDistribution().ave() << endl;
  cout << "# RMSD (from distribution):    "   
       << props.getDistribution().rmsd() << endl;

  cout << "# real average\t\t" << average / steps << endl;

  //cout << "# time\t\t" <<  props.toTitle() << endl;

  if (normalize)
    props.getDistribution().write_normalized(cout);
  else
    props.getDistribution().write(cout);

  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

