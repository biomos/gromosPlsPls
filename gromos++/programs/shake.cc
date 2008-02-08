// try using SHAKE from GromosXX
// markus christen

#include "../config.h"

#ifdef HAVE_GROMOSXX

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/fit/Reference.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Vec.h"

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <math.h>
#include <sstream>

// the gromosXX stuff

#include <gromosXX/math/gmath.h>
#include <gromosXX/util/debug.h>
#include <gromosXX/util/timing.h>
#include <gromosXX/util/error.h>

#include <gromosXX/io/message.h>
#include <gromosXX/io/argument.h>

#include <gromosXX/algorithm/algorithm.h>
#include <gromosXX/topology/topology.h>
#include <gromosXX/simulation/simulation.h>
#include <gromosXX/configuration/configuration.h>

#include <gromosXX/algorithm/algorithm/algorithm_sequence.h>

#include <gromosXX/interaction/interaction.h>
#include <gromosXX/interaction/interaction_types.h>
#include <gromosXX/interaction/forcefield/forcefield.h>

#include <gromosXX/util/parse_verbosity.h>

#include <gromosXX/io/topology/in_topology.h>

#include <gromosXX/algorithm/constraints/create_constraints.h>
#include <gromosXX/util/create_simulation.h>

using namespace std;
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

    //==================================================
    //=== XX INITIALIZATION                           ==
    //==================================================

    // debug_level = 100;

    util::simulation_struct a_xx_sim;

    a_xx_sim.sim.param().system.npm = 1;
    a_xx_sim.sim.param().system.nsm = 0;
    a_xx_sim.sim.param().constraint.ntc = 3;
    a_xx_sim.sim.param().constraint.solute.algorithm = simulation::constr_shake;
    a_xx_sim.sim.param().constraint.solvent.algorithm = simulation::constr_off;
    a_xx_sim.sim.param().constraint.solute.shake_tolerance = 0.0001;
    a_xx_sim.sim.param().constraint.solvent.shake_tolerance = 0.0001;

    
    for(int m=0, t=0; m < sys.numMolecules(); t+=sys.mol(m++).numAtoms()){
      a_xx_sim.sim.param().submolecules.submolecules.push_back(t);
    }
    
    {
      // create a XX In_Topology
      io::In_Topology in_topo;
      
      // no output...
      in_topo.quiet = true;
    
      if (util::create_simulation(args["topo"],
				  "",
				  "",
				  "",
				  a_xx_sim,
				  in_topo)){
	std::cerr << "creating the XX system failed!" << std::endl;
	return 1;
      }
    
      if (algorithm::create_constraints(a_xx_sim.md,
					a_xx_sim.topo,
					a_xx_sim.sim,
					in_topo,
					true)){
	std::cerr << "creating the constraints algorithm failed!" << std::endl;
	return 1;
      }
      
    } // don't need the In_Topology any more...
    
    a_xx_sim.conf.resize(a_xx_sim.topo.num_atoms());

    if (args["pbc"] == "r")
      a_xx_sim.conf.boundary_type = math::rectangular;
    else if (args["pbc"] == "v")
      a_xx_sim.conf.boundary_type = math::vacuum;
    else
      throw std::string("wrong boundary condition");

    // initialise arrays but do not gather
    a_xx_sim.conf.init(a_xx_sim.topo,
                       a_xx_sim.sim.param(),
                       false);

    //==================================================
    //=== XX INITIALIZED                              ==
    //==================================================
      
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);
      
    // define input coordinate
    InG96 ic;
    OutG96 oc(cout);
    
    // loop over all trajectories
    for(Arguments::const_iterator 
	  iter=args.lower_bound("traj"),
	  to=args.upper_bound("traj");
	iter!=to; ++iter){

      // open file
      ic.open((iter->second).c_str());

      // title
      std::ostringstream is;
      is << ic.title() << "\n\tshaken by GromosXX";
      oc.writeTitle(is.str());
      
      // loop over single trajectory
      while(!ic.eof()){
	ic >> sys;
	(*pbc.*gathmethod)();

	// parse the positions over
	for(int m=0, t=0; m < sys.numMolecules(); ++m){
	  for(int a=0; a < sys.mol(m).numAtoms(); ++a, ++t){
	    for(int d=0; d<3; ++d){
	      a_xx_sim.conf.current().pos(t)(d) = sys.mol(m).pos(a)[d];
	    }
	  }
	}
	// exchange them (for shake, we need new and old)
	a_xx_sim.conf.exchange_state();
	a_xx_sim.conf.current().pos = a_xx_sim.conf.old().pos;

	// set the box block
	a_xx_sim.conf.current().box(0) = math::Vec(sys.box()[0], 0.0, 0.0);
	a_xx_sim.conf.current().box(1) = math::Vec(0.0, sys.box()[1], 0.0);
	a_xx_sim.conf.current().box(2) = math::Vec(0.0, 0.0, sys.box()[2]);
	
	// shake
	a_xx_sim.md.run(a_xx_sim.topo, a_xx_sim.conf, a_xx_sim.sim);

	// recover positions
	// parse the positions over
	for(int m=0, t=0; m < sys.numMolecules(); ++m){
	  for(int a=0; a < sys.mol(m).numAtoms(); ++a, ++t){
	    for(int d=0; d<3; ++d){
	      sys.mol(m).pos(a)[d] = a_xx_sim.conf.current().pos(t)(d);
	    }
	  }
	}
      
	time += dt;
	oc << sys;
	
      }
	
    }

    ic.close();
    oc.close();
    
  }
  catch (const gromos::Exception &e){
    cerr << "EXCEPTION:\t";
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

// no GromosXX
#else

#include <iostream>

int main()
{
  std::cout << "\nconfigure could not find the GromosXX libraries" << std::endl
	    << "needed to run this program." << std::endl << std::endl
	    << "You need to add them to your CPPFLAGS, CXXFLAGS, LDFLAGS" << std::endl
            << "or run ./configure --with-gromosxx=<path>" << std::endl << std::endl
	    << "Reconfigure and recompile to use this program" << std::endl;
  return 1;
}

#endif
