/* gathtraj.cc
 * write a gathered trajectory to stdout.
 * Vincent Kraeutler, March 2002.
 */

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include <vector>
#include <iomanip>
#include <iostream>

using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "traj"};
  int nknowns = 3;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@traj <trajectory files>\n";
 
  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    //  read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);

    // define input coordinate
    InG96 ic;
    OutCoordinates *oc;
    oc = new OutG96();
    oc->open(cout);  
    
    ic.open(args["traj"].c_str());
    oc->writeTitle(ic.title());
      
    // loop over single trajectory
    while(!ic.eof()){
      ic >> sys;
      pbc->gather();
      *oc << sys;
    }

    ic.close();
    oc->close();
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
