// bound_RectBox.t.cc

#include <cassert>
#include <cstdlib>
#include "../gio/InG96.h"
#include "../gcore/System.h"
#include "../gio/InTopology.h"
#include "../gio/OutG96.h"
#include "RectBox.h"
#include "../gcore/System.h"
#include "../gcore/Solvent.h"
#include <string>
#include <iostream>
#include "../gmath/Vec.h"


using namespace gcore;
using namespace gio;
using bound::RectBox;

using namespace std;

int main(int argc, char *argv[]) {
  if (argc != 3) {
    cerr << "Usage: " + string(argv[0]) + " <Topology> <Filename>\n";
    exit(1);
  }
  string top = argv[1];
  string file = argv[2];

  // read Simulation data

  InTopology it(top);
  System sys(it.system());
  Solvent sol = sys.sol(0);

  InG96 ic;


  ic.open(file);
  ic.select("ALL");
  cout << "sys.numSolvents: " << sys.numSolvents() << endl;
  cout << "sol.numCoords: " << sol.numPos() << endl;
  OutG96 oc;

  ic >> sys;
  cout << "sys.numSolvents: " << sys.numSolvents() << endl;
  RectBox pbc(&sys);
  cout << "sys.numSolvents after RectBox pbc(&sys): " << sys.numSolvents() << endl;

  oc.open(cout);
  oc.select("ALL");
  oc.writeTitle(ic.title());

  while (!ic.eof()) {
    pbc.coggather();
    oc << sys;
    ic >> sys;
  }
  pbc.coggather();
  oc << sys;
  return 0;
}
