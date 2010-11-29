// bound_TruncOct.t.cc

#include <cassert>
#include <cstdlib>
#include "../gio/InG96.h"
#include "../gcore/System.h"
#include "../gio/InTopology.h"
#include "../gio/OutG96.h"
#include "TruncOct.h"
#include <string>
#include <iostream>

using namespace gcore;
using namespace gio;
using bound::TruncOct;

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

  TruncOct pbc(&sys);

  InG96 ic;
  ic.open(file);
  OutG96 oc;

  oc.open(cout);
  oc.writeTitle(ic.title());

  while (!ic.eof()) {
    ic >> sys;
    pbc.gather();
    oc << sys;
  }
  return 0;
}
